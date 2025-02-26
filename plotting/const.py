
PLT_PARAMS = {
    'axes.labelsize': 16,  # X and Y label font size
    'axes.titlesize': 16,  # Title font size
    'xtick.labelsize': 14,  # X-tick labels font size
    'ytick.labelsize': 14,  # Y-tick labels font size
    'legend.fontsize': 14,   # Legend font size
    'lines.linewidth': 3,
}

    
PEAK_PERFORMANCE = {
    'fritz': 2765, # GFLOPS
    'alex':  37400,  # 37.4 TFLOPS
}

OEPRATIONS_FLOPS = {
    'cholesky_ns3': 
        {'name': f'potrf $n_s^3$',
        'flops': lambda ns, nb: 1/3*ns**3+1/2*ns**2+1/6*ns,
        },
    'triang_solve_ns3': 
        {'name': f'trsm $n_s^3$',
        'flops': lambda ns, nb: ns**3,
        },
    'triang_solve_ns2nb': 
        {'name': f'trsm $n_s^2n_b$',
        'flops': lambda ns, nb: ns**2*nb,
        },
    'dgemm_ns3': 
        {'name': f'gemm $n_s^3$',
        'flops': lambda ns, nb: 2*ns**3,
        },
    'dgemm_ns2nb': 
        {'name': f'gemm $n_s^2n_b$',
        'flops': lambda ns, nb: 2*ns**2*nb,
        },
    'dgemm_nsnb2': 
        {'name': f'gemm $n_sn_b^2$',
        'flops': lambda ns, nb: 2*ns*nb**2,
        },
    'dgemm_nb2ns': 
        {'name': f'gemm $n_b^2n_s$',
        'flops': lambda ns, nb: 2*ns*nb**2,
        },
    'matrix_vector_nsns': 
        {'name': f'gemv $n_sn_s$',
        'flops': lambda ns, nb: 2*ns**2,
        },
    'matrix_vector_nsnb': 
        {'name': f'gemv $n_sn_b$',
        'flops': lambda ns, nb: 2*ns*nb,
        },
    'matrix_vector_nbns': 
        {'name': f'gemv $n_bn_s$',
        'flops': lambda ns, nb: 2*ns*nb,
        },
    'dot_prod_ns': 
        {'name': f'dot $n_s$',
        'flops': lambda ns, nb: 2*ns,
        },
    'scale_ns': 
        {'name': f'scal $n_s$',
        'flops': lambda ns, nb: ns,
        },
}