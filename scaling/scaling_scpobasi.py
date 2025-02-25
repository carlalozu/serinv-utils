import numpy as np
import time
import argparse

from serinv.algs.work_in_progress.scpobaf import scpobaf
from serinv.algs.work_in_progress.scpobasi import scpobasi
from storage.utils_ba import dd_ba
from storage.parameters import calculate_parameters_banded
from flops.scpobaf_flops import scpobaf_flops
from flops.scpobasi_flops import scpobasi_flops

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


def run_scpobasi(n, n_offdiags, arrowhead_size, overwrite, dtype):
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

    # Cholesky decomposition serinv
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

    # Inversion serinv
    start_time = time.time()
    (
        I_diagonal,
        I_lower_diagonals,
        I_arrow_bottom,
        I_arrow_tip
    ) = scpobasi(
        L_diagonal,
        L_lower_diagonals,
        L_arrow_bottom,
        L_arrow_tip,
        overwrite
    )
    end_time = time.time()
    out += f"{end_time - start_time}"

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
    parser.add_argument('--overwrite', type=int, required=False, default=1,
                        help="Overwrite the original arrays.")

    # Parse arguments
    args = parser.parse_args()
    parameters = calculate_parameters_banded(
        matrix_size=args.n,
        arrowhead_width=args.arrowhead_blocksize,
        bandwidth=args.bandwidth,
    )

    dtype = np.float64
    overwrite = bool(args.overwrite)

    # print("n,bandwidth,arrowhead_blocksize,effective_bandwidth,diagonal_blocksize,n_offdiags,n_t,scpobaf_time,scpobasi_time,scpobaf_FLOPS,scpobasi_FLOPS")

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

        out = run_scpobasi(
            n=parameters['parameters']['matrix_size'],
            n_offdiags=parameters['parameters']['n_offdiags'],
            arrowhead_size=parameters['parameters']['arrowhead_blocksize'],
            overwrite=overwrite,
            dtype=dtype,
        )
        print(out, end=',')

        # GET FLOPS
        flops_c = scpobaf_flops(
            n_diagonals=parameters['parameters']['n_t'],
            n_offdiags=parameters['parameters']['n_offdiags'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
        )
        print(flops_c, end=',')

        flops_si = scpobasi_flops(
            n_diagonals=parameters['parameters']['n_t'],
            n_offdiags=parameters['parameters']['n_offdiags'],
            arrowhead_blocksize=parameters['parameters']['arrowhead_blocksize'],
        )
        print(flops_si)


if __name__ == "__main__":
    main()
