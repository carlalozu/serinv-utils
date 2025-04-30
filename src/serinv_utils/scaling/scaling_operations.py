from timeit import Timer
import argparse
import numpy as np
import scipy.linalg as np_la

try:
    import cupy as cp
    from serinv.cupyfix.cholesky_lowerfill import cholesky_lowerfill
    import cupyx.scipy.linalg as cu_la

    CUPY_AVAIL = True

    xp = cp
    la = cu_la
    cholesky = cholesky_lowerfill

except ImportError:

    CUPY_AVAIL = False
    xp = np
    la = np_la
    cholesky = np.linalg.cholesky


def spd(m, factor):
    # Make matrix positive def
    for i in range(m.shape[0]):
        m[i, i] = (1 + xp.sum(m[i, :]))*factor
    m[:, :] = (m[:, :] + m[:, :].conj())/2


def run_operations(diag_blocksize, arrowhead_blocksize, n_runs, dtype):

    print(n_runs, end=',')

    block_ns_ns = xp.zeros(
        (diag_blocksize, diag_blocksize),
        dtype=dtype
    )
    block_ns_ns[:, :] = xp.random.rand(
        *(diag_blocksize, diag_blocksize),
    )
    # Make matrix positive def
    spd(block_ns_ns, int(xp.sqrt(diag_blocksize)))
    block_ns_ns_out = xp.zeros(
        (diag_blocksize, diag_blocksize),
        dtype=dtype
    )
    block_nb_ns = xp.zeros(
        (arrowhead_blocksize, diag_blocksize),
        dtype=dtype
    )
    block_nb_ns[:, :] = xp.random.rand(
        *(arrowhead_blocksize, diag_blocksize),
    )
    block_nb_ns_out = xp.zeros(
        (arrowhead_blocksize, diag_blocksize),
        dtype=dtype
    )

    block_nb_nb = xp.zeros(
        (arrowhead_blocksize, arrowhead_blocksize),
        dtype=dtype
    )
    block_nb_nb[:, :] = xp.random.rand(
        *(arrowhead_blocksize, arrowhead_blocksize),
    )
    block_nb_nb_out = xp.zeros(
        (arrowhead_blocksize, arrowhead_blocksize),
        dtype=dtype
    )

    vector_ns = xp.random.rand(diag_blocksize)
    vector_ns_out = xp.zeros(diag_blocksize,
                             dtype=dtype)

    vector_nb = xp.random.rand(arrowhead_blocksize)
    vector_nb_out = xp.zeros(arrowhead_blocksize,
                             dtype=dtype)

    element = xp.random.rand(1)
    element_out = xp.zeros(1, dtype=dtype)
    print(f"{diag_blocksize},{arrowhead_blocksize}", end=',')

    # Do Cholesky
    def cholesky_op():
        block_ns_ns_out[:, :] = cholesky(block_ns_ns)
    time_taken = Timer(cholesky_op).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do triangular solve ns^3
    def triangular_solve_ns3():
        block_ns_ns_out[:, :] = la.solve_triangular(
            block_ns_ns, block_ns_ns.T, lower=True
        ).T
    time_taken = Timer(triangular_solve_ns3).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do triangular solve ns^2nb
    def triangular_solve_ns2nb():
        block_nb_ns_out[:, :] = la.solve_triangular(
            block_ns_ns, block_nb_ns.T, lower=True
        ).T
    time_taken = Timer(triangular_solve_ns2nb).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMM_ns^3
    def dgemm_ns3():
        block_ns_ns_out[:, :] = block_ns_ns.T @ block_ns_ns
    time_taken = Timer(dgemm_ns3).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMM_ns^2nb
    def dgemm_ns2nb():
        block_nb_ns_out[:, :] = block_nb_ns @ block_ns_ns.T
    time_taken = Timer(dgemm_ns2nb).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMM_nb^2ns
    def dgemm_nb2ns():
        block_nb_nb_out[:, :] = block_nb_ns @ block_nb_ns.T
    time_taken = Timer(dgemm_nb2ns).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMM_nsnb^2
    def dgemm_nsnb2():
        block_nb_ns_out[:, :] = block_nb_nb @ block_nb_ns
    time_taken = Timer(dgemm_nsnb2).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMV_ns^2
    def dgemv_ns2():
        vector_ns_out[:] = block_ns_ns @ vector_ns
    time_taken = Timer(dgemv_ns2).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMV_nsnb
    def dgemv_nsnb():
        vector_ns_out[:] = block_nb_ns.T @ vector_nb
    time_taken = Timer(dgemv_nsnb).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do DGEMV_nbns
    def dgemv_nbns():
        vector_nb_out[:] = block_nb_ns @ vector_ns
    time_taken = Timer(dgemv_nbns).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do ns \cdot ns
    def ns_dot_ns():
        element_out[0] = vector_ns.T @ vector_ns
    time_taken = Timer(ns_dot_ns).timeit(number=n_runs)
    print(time_taken, end=',')

    # Do scale ns
    def scale_ns():
        vector_ns_out[:] = vector_ns.T * element[0]
    time_taken = Timer(scale_ns).timeit(number=n_runs)
    print(time_taken)


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Configure parameters.")
    parser.add_argument('--diag_blocksize', type=int, required=True,
                        help="Diagonal block width.")
    parser.add_argument('--arrowhead_blocksize', type=int, required=True,
                        help="Arrowhead block width.")
    parser.add_argument('--n_runs', type=int, required=True,
                        help="n_runs for the algorithms.")

    # Parse arguments
    args = parser.parse_args()

    # Set values
    dtype = np.float64

    # print("diag_blocksize,arrowhead_blocksize,cholesky_ns3,triang_solve_ns3,triang_solve_ns2nb,dgemm_ns3,dgemm_ns2nb,dgemm_nb2ns,dgemm_nsnb2,matrix_vector_nsns,matrix_vector_nsnb,matrix_vector_nbns,dot_prod_ns,scale_ns")
    run_operations(
        diag_blocksize=args.diag_blocksize,
        arrowhead_blocksize=args.arrowhead_blocksize,
        n_runs=args.n_runs,
        dtype=dtype,
    )


if __name__ == "__main__":
    main()
