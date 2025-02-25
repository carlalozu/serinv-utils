import time
import argparse
import numpy as np

from serinv.algs import pobtaf
from storage.utils_bba import dd_bba, bba_arrays_to_dense, bba_dense_to_arrays
from storage.parameters import calculate_parameters_tri_diagonal
from flops.flops import T_flops_POBTAF

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


def run_pobtaf(
        diagonal_blocksize, arrowhead_blocksize, n_t,
        dtype, fits_memory=False, streaming=False):

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

    if fits_memory:
        # Get dense matrix
        M = bba_arrays_to_dense(
            xp.copy(M_diagonal_blocks),
            xp.copy(M_lower_diagonal_blocks),
            xp.copy(M_arrow_bottom_blocks),
            xp.copy(M_arrow_tip_block),
            symmetric=True
        )
        start_time = time.time()
        I_ref = cholesky(M)
        end_time = time.time()
        out += f"{end_time - start_time},"  # numpy time

    # Do inversion on compressed format
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
        streaming,
    )
    end_time = time.time()
    out += f"{end_time - start_time},"  # time

    if fits_memory:
        # From compressed to dense
        (
            L_ref_diagonal_blocks,
            L_ref_lower_diagonal_blocks,
            L_ref_arrow_bottom_blocks,
            L_ref_arrow_tip_block
        ) = bba_dense_to_arrays(
            I_ref,
            1,
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
        out += f"{error.item()}"  # error

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
    parser.add_argument('--streaming', type=int, required=False, default=0,
                        help="Device streaming.")
    parser.add_argument('--numpy_compare', type=int, required=False, default=0,
                        help="Fits memory and compare against numpy.")

    # Parse arguments
    args = parser.parse_args()
    parameters = calculate_parameters_tri_diagonal(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
    )

    # Set values
    dtype = np.float64
    numpy_compare = bool(args.numpy_compare)
    streaming = bool(args.streaming)

    # print("n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error")

    print(parameters['parameters']['matrix_size'], end=',')
    print(parameters['parameters']['bandwidth'], end=',')
    print(parameters['parameters']['arrowhead_blocksize'], end=',')

    if not parameters['flag']:
        print('NA,NA,NA,NA,NA,NA,NA')
        print(parameters['error'])

    else:
        print(parameters['parameters']['effective_bandwidth'], end=',')
        print(parameters['parameters']['diagonal_blocksize'], end=',')
        print(parameters['parameters']['n_offdiags'], end=',')
        print(parameters['parameters']['n_t'], end=',')

        out = run_pobtaf(
            diagonal_blocksize=parameters['parameters']['diagonal_blocksize'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
            n_t=parameters['parameters']['n_t'],
            dtype=dtype,
            fits_memory=numpy_compare,
            streaming=streaming,
        )
        print(out, end=',')

        # GET FLOPS
        flops = T_flops_POBTAF(
            nt=parameters['parameters']['n_t'],
            ns=parameters['parameters']['diagonal_blocksize'],
            nb=parameters['parameters']['arrowhead_blocksize']
        )
        print(flops)


if __name__ == "__main__":
    main()
