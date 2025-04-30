import numpy as np
import time

from serinv.algs.work_in_progress.scpobbaf import scpobbaf_c
from serinv.algs.work_in_progress.scpobbasi import scpobbasi_c
from serinv_utils.storage.utils_bba import dd_bba, fill_bba
from serinv_utils.storage.parameters import calculate_parameters_n_diagonal
from serinv_utils.flops.scpobbaf_flops import scpobbaf_flops
from serinv_utils.flops.scpobbasi_flops import scpobbasi_flops

import argparse

try:
    import cupy as cp
    import cupyx.scipy.linalg as cu_la

    CUPY_AVAIL = True
    xp = cp
    la = cu_la

except ImportError:

    import scipy.linalg as np_la

    CUPY_AVAIL = False
    xp = np
    la = np_la


def run_scpobbasi(
    M_diagonal_blocks,
    M_lower_diagonal_blocks,
    M_arrow_bottom_blocks,
    M_arrow_tip_block,
    overwrite=True
):

    start_time = time.time()
    (
        L_diagonal_blocks,
        L_lower_diagonal_blocks,
        L_arrow_bottom_blocks,
        L_arrow_tip_block
    ) = scpobbaf_c(
        M_diagonal_blocks,
        M_lower_diagonal_blocks,
        M_arrow_bottom_blocks,
        M_arrow_tip_block,
        overwrite,
    )
    end_time = time.time()
    time_c = end_time - start_time

    # Do inversion on compressed format
    start_time = time.time()
    (
        I_diagonal_blocks,
        I_lower_diagonal_blocks,
        I_arrow_bottom_blocks,
        I_arrow_tip_block
    ) = scpobbasi_c(
        L_diagonal_blocks,
        L_lower_diagonal_blocks,
        L_arrow_bottom_blocks,
        L_arrow_tip_block,
        overwrite,
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
    parser.add_argument('--n_offdiags_blk', type=int, required=True,
                        help="Number of lower off-diagonal blocks.")
    parser.add_argument('--overwrite', type=int, required=False, default=1,
                        help="Overwrite the original arrays.")
    parser.add_argument('--n_runs', type=int, required=False, default=1,
                        help="Number of runs.")

    # Parse arguments
    args = parser.parse_args()

    # Set values
    dtype = np.float64
    overwrite = bool(args.overwrite)

    # Calculate parameters
    parameters = calculate_parameters_n_diagonal(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
        n_offdiags_=args.n_offdiags_blk,
    )

    # print("n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error")

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
            parameters['parameters']['n_offdiags'],
            parameters['parameters']['diagonal_blocksize'],
            parameters['parameters']['arrowhead_blocksize'],
            n_t=parameters['parameters']['n_t'],
            dtype=dtype,
        )

        time_c = 0
        time_in = 0
        for _ in range(args.n_runs):

            times = run_scpobbasi(
                M_diagonal_blocks,
                M_lower_diagonal_blocks,
                M_arrow_bottom_blocks,
                M_arrow_tip_block,
                overwrite
            )
            time_c += times[0]
            time_in += times[1]

            fill_bba(
                M_diagonal_blocks,
                M_lower_diagonal_blocks,
                M_arrow_bottom_blocks,
                M_arrow_tip_block,
                factor=int(xp.sqrt(parameters['parameters']['n_t']))
            )

        print(time_c, end=',')
        print(time_in, end=',')

        flops_c, _ = scpobbaf_flops(
            n_diag_blocks=parameters['parameters']['n_t'],
            diagonal_blocksize=parameters['parameters']['diagonal_blocksize'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
            n_offdiags_blk=parameters['parameters']['n_offdiags'],
        )
        print(flops_c, end=',')

        flops_si, _ = scpobbasi_flops(
            n_diag_blocks=parameters['parameters']['n_t'],
            diagonal_blocksize=parameters['parameters']['diagonal_blocksize'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
            n_offdiags_blk=parameters['parameters']['n_offdiags'],
        )
        print(flops_si)


if __name__ == "__main__":
    main()
