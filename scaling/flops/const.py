
OPERATIONS_FLOPS = {
    'cholesky_ns3': {
        'name': 'POTRF$n_s^3$',
        'flops': lambda ns, nb: 1/3*ns**3+1/2*ns**2+1/6*ns,
    },
    'cholesky_nb3': {
        'name': 'POTRF$n_b^3$',
        'flops': lambda ns, nb: 1/3*nb**3+1/2*nb**2+1/6*nb,
    },
    'triang_solve_ns3': {
        'name': 'TRSM$n_s^3$',
        'flops': lambda ns, nb: ns**3,
    },
    'triang_solve_nb3': {
        'name': 'TRSM$n_b^3$',
        'flops': lambda ns, nb: nb**3,
    },
    'triang_solve_ns2nb': {
        'name': 'TRSM$n_s^2n_b$',
        'flops': lambda ns, nb: ns**2*nb,
    },
    'dgemm_ns3': {
        'name': 'GEMM$n_s^3$',
        'flops': lambda ns, nb: 2*ns**3,
    },
    'dgemm_nb3': {
        'name': 'GEMM$n_b^3$',
        'flops': lambda ns, nb: 2*nb**3,
    },
    'dgemm_ns2nb': {
        'name': 'GEMM$n_s^2n_b$',
        'flops': lambda ns, nb: 2*ns**2*nb,
    },
    'dgemm_nsnb2': {
        'name': 'GEMM$n_sn_b^2$',
        'flops': lambda ns, nb: 2*ns*nb**2,
    },
    'dgemm_nb2ns': {
        'name': 'GEMM$n_b^2n_s$',
        'flops': lambda ns, nb: 2*ns*nb**2,
    },
    'matrix_vector_nsns': {
        'name': 'GEMV$n_s^2$',
        'flops': lambda ns, nb: 2*ns**2,
    },
    'matrix_vector_nsns_1': {
        'name': 'GEMV$n_sn_{s-1}$',
        'flops': lambda ns, nb: 2*ns*(ns-1),
    },
    'matrix_vector_nsnb': {
        'name': 'GEMV$n_sn_b$',
        'flops': lambda ns, nb: 2*ns*nb,
    },
    'matrix_vector_nbnb': {
        'name': 'GEMV$n_b^2$',
        'flops': lambda ns, nb: 2*nb**2,
    },
    'dot_prod_ns': {
        'name': 'DOT$n_s$',
        'flops': lambda ns, nb: 2*ns,
    },
    'dot_prod_nb': {
        'name': 'DOT$n_b$',
        'flops': lambda ns, nb: 2*nb,
    },
    'scale_ns': {
        'name': 'SCAL$n_s$',
        'flops': lambda ns, nb: ns,
    },
    'scale_nb': {
        'name': 'SCAL$n_b$',
        'flops': lambda ns, nb: nb,
    },
    'div': {
        'name': 'div',
        'flops': lambda ns, nb: 1,
    },
    'sqrt': {
        'name': 'sqrt',
        'flops': lambda ns, nb: 1,
    },
    'mult': {
        'name': 'mult',
        'flops': lambda ns, nb: 1,
    },
    'ger_nb2': {
        'name': 'GER$n_b^2$',
        'flops': lambda ns, nb: 2*nb**2,
    },
}

ALG_OPERATION_COUNT = {
    'POBAF': {
        'scale_ns': lambda nt, n: nt-1,
        'scale_nb': lambda nt, n: nt,
        'div': lambda nt, n: nt,
        'sqrt': lambda nt, n: nt,
        'dot_prod_ns': lambda nt, n: nt-1,
        'matrix_vector_nsns_1': lambda nt, n: nt-1,
        'matrix_vector_nsnb': lambda nt, n: nt-1,
        'ger_nb2': lambda nt, n: nt,
        'cholesky_nb3': lambda nt, n: 1,
    },
    'POBASI': {
        'triang_solve_nb3': lambda nt, n: 1,
        'scale_ns': lambda nt, n: nt-1,
        'scale_nb': lambda nt, n: nt,
        'mult': lambda nt, n: nt,
        'div': lambda nt, n: nt,
        'dot_prod_ns': lambda nt, n: nt-1,
        'dot_prod_nb': lambda nt, n: nt,
        'matrix_vector_nsns': lambda nt, n: nt-1,
        'matrix_vector_nsnb': lambda nt, n: nt-1,
        'matrix_vector_nbnb': lambda nt, n: nt,
        'dgemm_nb3': lambda nt, n: 1,
    },
    'POBTAF': {
        'dgemm_ns3': lambda nt, n: nt-1,
        'dgemm_ns2nb': lambda nt, n: nt-1,
        'dgemm_nb2ns': lambda nt, n: nt,
        'triang_solve_ns3': lambda nt, n: nt-1,
        'triang_solve_ns2nb': lambda nt, n: nt,
        'cholesky_ns3': lambda nt, n: nt,
        'cholesky_nb3': lambda nt, n: 1,
    },
    'POBTASI': {
        'triang_solve_ns3': lambda nt, n: nt,
        'triang_solve_nb3': lambda nt, n: 1,
        'dgemm_ns3': lambda nt, n: 4*(nt-1)+1,
        'dgemm_ns2nb': lambda nt, n: 4*(nt-1)+2,
        'dgemm_nsnb2': lambda nt, n: nt,
        'dgemm_nb3': lambda nt, n: 1,
    },
    'POBBAF': {
        'dgemm_ns3': lambda nt, n: (nt-1)*n**2,
        'dgemm_ns2nb': lambda nt, n: (nt-1)*n,
        'dgemm_nb2ns': lambda nt, n: nt,
        'triang_solve_ns3': lambda nt, n: (nt-1)*n,
        'triang_solve_ns2nb': lambda nt, n: nt,
        'cholesky_ns3': lambda nt, n: nt,
        'cholesky_nb3': lambda nt, n: 1,
    },
    'POBBASI': {
        'triang_solve_ns3': lambda nt, n: nt,
        'triang_solve_nb3': lambda nt, n: 1,
        'dgemm_ns3': lambda nt, n: (nt-1)*(n+1)**2+1,
        'dgemm_ns2nb': lambda nt, n: 2*(nt-1)*(n+1)+2,
        'dgemm_nsnb2': lambda nt, n: nt,
        'dgemm_nb3': lambda nt, n: 1,
    },
}
