import numpy as np
import time

from serinv.algs.work_in_progress.scpobaf import scpobaf
from storage.utils_ba import dd_ba, ba_arrays_to_dense, ba_dense_to_arrays
from storage.parameters import calculate_parameters_banded
from flops.flops import scpobaf_flops

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


def run_scpobaf(n, n_offdiags, arrowhead_size, overwrite, dtype, fits_memory=False):

    out = ""

    (
        M_diagonal,
        M_lower_diagonals,
        M_arrow_bottom,
        M_arrow_tip
    ) = dd_ba(
        n_offdiags,
        arrowhead_size,
        n,
        dtype,
        factor=int(np.sqrt(n))
    )

    if fits_memory:
        M = ba_arrays_to_dense(
            xp.copy(M_diagonal),
            xp.copy(M_lower_diagonals),
            xp.copy(M_arrow_bottom),
            xp.copy(M_arrow_tip),
            symmetric=True
        )
        # Cholesky numpy
        start_time = time.time()
        L_ref = cholesky(M)
        end_time = time.time()
        out += f"{end_time - start_time},"

    # Cholesky serinv
    start_time = time.time()
    (
        L_diagonal,
        L_lower_diagonals,
        L_arrow_bottom,
        L_arrow_tip
    ) = scpobaf(
        M_diagonal,
        M_lower_diagonals,
        M_arrow_bottom,
        M_arrow_tip,
        overwrite
    )
    end_time = time.time()
    out += f"{end_time - start_time},"

    if fits_memory:

        (
            L_ref_diagonal,
            L_ref_lower_diagonals,
            L_ref_arrow_bottom,
            L_ref_arrow_tip
        ) = ba_dense_to_arrays(L_ref, n_offdiags, arrowhead_size)

        # Error
        error = (
            xp.mean(L_ref_diagonal - L_diagonal) +
            xp.mean(L_ref_lower_diagonals - L_lower_diagonals) +
            xp.mean(L_ref_arrow_bottom - L_arrow_bottom) +
            xp.mean(L_ref_arrow_tip - L_arrow_tip)
        )/4

        out += f"{error.item()}"
    else:
        out += "NA,NA"

    return out


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Configure matrix parameters.")
    parser.add_argument('--n', type=int, required=True, help="Matrix size.")
    parser.add_argument('--bandwidth', type=int,
                        required=True, help="Bandwidth.")
    parser.add_argument('--arrowhead_blocksize', type=int,
                        required=True, help="Arrowhead width.")
    parser.add_argument('--numpy_compare', type=int, required=False, default=0,
                        help="Fits memory and compare against numpy.")
    parser.add_argument('--overwrite', type=int, required=False, default=1,
                        help="Overwrite the original arrays.")

    # Parse arguments
    args = parser.parse_args()
    parameters = calculate_parameters_banded(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
    )

    overwrite = True
    dtype = np.float64
    numpy_compare = bool(args.numpy_compare)
    overwrite = bool(args.overwrite)

    # print("n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,time,numpy_time,error,FLOPS")

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

        out = run_scpobaf(
            n=parameters['parameters']['matrix_size'],
            n_offdiags=parameters['parameters']['n_offdiags'],
            arrowhead_size=parameters['parameters']['arrowhead_blocksize'],
            overwrite=overwrite,
            dtype=dtype,
            fits_memory=numpy_compare
        )
        print(out, end=',')

        # GET FLOPS
        flops = scpobaf_flops(
            n_diag_blocks=parameters['parameters']['n_t'],
            n_offdiags=parameters['parameters']['n_offdiags'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
        )
        print(flops)


if __name__ == "__main__":
    main()
