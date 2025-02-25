import numpy as np
import time

from serinv.algs.work_in_progress.scpobbaf import scpobbaf_c
from ..storage.utils_bba import dd_bba, bba_arrays_to_dense, bba_dense_to_arrays
from ..storage.parameters import calculate_parameters_n_diagonal
from ..flops.scpobbaf_flops import scpobbaf_flops

import argparse

try:
    import cupy as cp
    from serinv.cupyfix.cholesky_lowerfill import cholesky_lowerfill

    CUPY_AVAIL = True

    xp = cp
    cholesky = cholesky_lowerfill

except ImportError:

    CUPY_AVAIL = False

    xp = np
    cholesky = np.linalg.cholesky


def run_pobbaf(
        n_offdiags_blk, diagonal_blocksize, arrowhead_blocksize, n_diag_blocks,
        dtype, fits_memory=False, overwrite=False):
    out = ""

    # Create matrix in compressed format
    (
        M_diagonal_blocks,
        M_lower_diagonal_blocks,
        M_arrow_bottom_blocks,
        M_arrow_tip_block
    ) = dd_bba(
        n_offdiags_blk,
        diagonal_blocksize,
        arrowhead_blocksize,
        n_diag_blocks,
        dtype,
    )

    # Do inversion on compressed format
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
    out += f"{end_time - start_time},"

    if fits_memory:
        # Get dense matrix
        M = bba_arrays_to_dense(
            M_diagonal_blocks,
            M_lower_diagonal_blocks,
            M_arrow_bottom_blocks,
            M_arrow_tip_block,
            symmetric=True
        )
        start_time = time.time()
        L_ref = cholesky(M)
        end_time = time.time()
        out += f"{end_time - start_time},"

        # From compressed to dense
        (
            L_ref_diagonal_blocks,
            L_ref_lower_diagonal_blocks,
            L_ref_arrow_bottom_blocks,
            L_ref_arrow_tip_block
        ) = bba_dense_to_arrays(
            L_ref,
            n_offdiags_blk,
            diagonal_blocksize,
            arrowhead_blocksize,
            lower=False
        )

        # Error between dense decomposition and numpy
        error = (
            xp.mean(L_ref_diagonal_blocks - L_diagonal_blocks) +
            xp.mean(L_ref_lower_diagonal_blocks - L_lower_diagonal_blocks) +
            xp.mean(L_ref_arrow_bottom_blocks - L_arrow_bottom_blocks) +
            xp.mean(L_ref_arrow_tip_block - L_arrow_tip_block)
        )/4
        out += f"{error.item()}"

    else:
        out += "NA,NA"

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
    parser.add_argument('--n_offdiags_blk', type=int, required=True,
                        help="Number of lower off-diagonal blocks.")
    parser.add_argument('--numpy_compare', type=int, required=False, default=0,
                        help="Fits memory and compare against numpy.")
    parser.add_argument('--overwrite', type=int, required=False, default=1,
                        help="Overwrite the original arrays.")

    # Parse arguments
    args = parser.parse_args()

    # Set values
    dtype = np.float64
    numpy_compare = bool(args.numpy_compare)
    overwrite = bool(args.overwrite)

    # Calculate parameters
    parameters = calculate_parameters_n_diagonal(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
        n_offdiags_=args.n_offdiags_blk,
    )

    # print("n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error")

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

        out = run_pobbaf(
            n_offdiags_blk=parameters['parameters']['n_offdiags'],
            diagonal_blocksize=parameters['parameters']['diagonal_blocksize'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
            n_diag_blocks=parameters['parameters']['n_t'],
            dtype=dtype,
            fits_memory=numpy_compare,
            overwrite=overwrite,
        )
        print(out, end=',')

        # GET FLOPS
        flops = scpobbaf_flops(
            n_diag_blocks=parameters['parameters']['n_t'],
            diagonal_blocksize=parameters['parameters']['diagonal_blocksize'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
            n_offdiags_blk=parameters['parameters']['n_offdiags'],
        )
        print(flops)


if __name__ == "__main__":
    main()
