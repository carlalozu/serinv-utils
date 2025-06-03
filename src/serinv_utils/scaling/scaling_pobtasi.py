import time
import argparse

import numpy as np
import scipy.linalg as np_la

from serinv.algs import pobtaf, pobtasi
from serinv_utils.storage.utils_bba import dd_bba, fill_bba
from serinv_utils.storage.parameters import calculate_parameters_tri_diagonal
from serinv_utils.flops.flops import T_flops_POBTAF, T_flops_POBTASI

try:
    import cupy as cp
    import cupyx.scipy.linalg as cu_la

    CUPY_AVAIL = True
    xp = cp
    la = cu_la

except ImportError:

    CUPY_AVAIL = False
    xp = np
    la = np_la


def run_pobtasi(
    M_diagonal_blocks,
    M_lower_diagonal_blocks,
    M_arrow_bottom_blocks,
    M_arrow_tip_block,
):

    # Do inversion on compressed format
    # pobtaf time
    start_time = time.time()
    (
        L_diagonal_blocks,
        L_lower_diagonal_blocks,
        L_arrow_bottom_blocks,
        L_arrow_tip_block
    ) = pobtaf(
        M_diagonal_blocks,
        M_lower_diagonal_blocks,
        M_arrow_bottom_blocks,
        M_arrow_tip_block,
        False,
    )
    end_time = time.time()
    time_c = end_time - start_time

    # pobtasi time
    start_time = time.time()
    (
        I_diagonal_blocks,
        I_lower_diagonal_blocks,
        I_arrow_bottom_blocks,
        I_arrow_tip_block
    ) = pobtasi(
        L_diagonal_blocks,
        L_lower_diagonal_blocks,
        L_arrow_bottom_blocks,
        L_arrow_tip_block,
        False,
    )
    end_time = time.time()
    time_in = end_time - start_time

    return time_c, time_in


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Configure parameters.")
    parser.add_argument('--n', type=int, required=True,
                        help="Matrix size.")
    parser.add_argument('--bandwidth', type=int, required=True,
                        help="Bandwidth.")
    parser.add_argument('--arrowhead_blocksize', type=int, required=True,
                        help="Arrowhead block width.")
    parser.add_argument('--n_runs', type=int, required=False, default=1,
                        help="Number of runs.")

    # Parse arguments
    args = parser.parse_args()
    parameters = calculate_parameters_tri_diagonal(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
    )

    # Set values
    dtype = np.float64


    print(args.n_runs, end=',')
    print(parameters['parameters']['matrix_size'], end=',')
    print(parameters['parameters']['bandwidth'], end=',')
    print(parameters['parameters']['arrowhead_blocksize'], end=',')

    if not parameters['flag']:
        print('NA,NA,NA,NA,NA,NA,NA')

    else:
        print(parameters['parameters']['effective_bandwidth'], end=',')
        print(parameters['parameters']['diagonal_blocksize'], end=',')
        print(parameters['parameters']['n_offdiags'], end=',')
        print(parameters['parameters']['n_t'], end=',')

        (
            M_diagonal_blocks,
            M_lower_diagonal_blocks,
            M_arrow_bottom_blocks,
            M_arrow_tip_block,
        ) = dd_bba(
            1,
            parameters['parameters']['diagonal_blocksize'],
            parameters['parameters']['arrowhead_blocksize'],
            n_t=parameters['parameters']['n_t'],
            dtype=dtype,
        )


        # Warm up runs
        warmup_runs = 3
        if args.n_runs < warmup_runs:
            warmup_runs = 1
        for _ in range(warmup_runs):
            times = run_pobtasi(
                M_diagonal_blocks,
                M_lower_diagonal_blocks,
                M_arrow_bottom_blocks,
                M_arrow_tip_block,
            )
            fill_bba(
                M_diagonal_blocks,
                M_lower_diagonal_blocks,
                M_arrow_bottom_blocks,
                M_arrow_tip_block,
                factor=int(xp.sqrt(parameters['parameters']['n_t']))
            )

        time_c = []
        time_in = []
        for _ in range(args.n_runs):

            times = run_pobtasi(
                M_diagonal_blocks,
                M_lower_diagonal_blocks,
                M_arrow_bottom_blocks,
                M_arrow_tip_block,
            )
            time_c.append(times[0])
            time_in.append(times[1])

            fill_bba(
                M_diagonal_blocks,
                M_lower_diagonal_blocks,
                M_arrow_bottom_blocks,
                M_arrow_tip_block,
                factor=int(xp.sqrt(parameters['parameters']['n_t']))
            )

        ## Print results
        print(np.median(time_c), end=',')
        print(np.std(time_c), end=',')

        print(np.median(time_in), end=',')
        print(np.std(time_in), end=',')

        # GET FLOPS
        flops_c = T_flops_POBTAF(
            nt=parameters['parameters']['n_t'],
            ns=parameters['parameters']['diagonal_blocksize'],
            nb=parameters['parameters']['arrowhead_blocksize']
        )
        print(flops_c, end=',')

        flops_si = T_flops_POBTASI(
            nt=parameters['parameters']['n_t'],
            ns=parameters['parameters']['diagonal_blocksize'],
            nb=parameters['parameters']['arrowhead_blocksize']
        )
        print(flops_si)


if __name__ == "__main__":
    main()
