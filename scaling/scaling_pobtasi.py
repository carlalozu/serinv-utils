import time
import argparse

import numpy as np
import scipy.linalg as np_la

from serinv.algs import pobtaf, pobtasi
from storage.utils_bba import dd_bba, bba_arrays_to_dense, bba_dense_to_arrays
from storage.parameters import calculate_parameters_tri_diagonal
from flops.flops import T_flops_POBTAF, T_flops_POBTASI

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


def run_pobtasi(diagonal_blocksize, arrowhead_blocksize, n_t, dtype):

    out = ""

    # Create matrix in compressed format
    (
        M_diagonal_blocks,
        M_lower_diagonal_blocks,
        M_arrow_bottom_blocks,
        M_arrow_tip_block
    ) = dd_bba(
        1,
        diagonal_blocksize,
        arrowhead_blocksize,
        n_t,
        dtype,
    )

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
    out += f"{end_time - start_time},"

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
    out += f"{end_time - start_time}"

    return out


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Configure parameters.")
    parser.add_argument('--n', type=int, required=True,
                        help="Matrix size.")
    parser.add_argument('--bandwidth', type=int, required=True,
                        help="Bandwidth.")
    parser.add_argument('--arrowhead_blocksize', type=int, required=True,
                        help="Arrowhead block width.")

    # Parse arguments
    args = parser.parse_args()
    parameters = calculate_parameters_tri_diagonal(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
    )

    # Set values
    dtype = np.float64

    # print("n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,pobtaf_time,pobtasi_time,pobtaf_FLOPS,pobtasi_FLOPS")

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

        out = run_pobtasi(
            diagonal_blocksize=parameters['parameters']['diagonal_blocksize'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
            n_t=parameters['parameters']['n_t'],
            dtype=dtype,
        )
        print(out, end=',')

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
