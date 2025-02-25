import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

# plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update({
    'axes.labelsize': 16,  # X and Y label font size
    'axes.titlesize': 16,  # Title font size
    'xtick.labelsize': 14,  # X-tick labels font size
    'ytick.labelsize': 14,  # Y-tick labels font size
    'legend.fontsize': 14,   # Legend font size
    'lines.linewidth': 3,
})

operations = {
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

operations_block_chol = [
    'cholesky_ns3',
    'triang_solve_ns3',
    'triang_solve_ns2nb',
    'dgemm_ns3',
    'dgemm_ns2nb',
    'dgemm_nb2ns',
]

operations_block_chol_inv = [
    'triang_solve_ns3',
    'dgemm_ns3',
    'dgemm_ns2nb',
    'dgemm_nb2ns',
]

operations_banded_chol = [
    'scale_ns',
    'dot_prod_ns',
    'matrix_vector_nsns',
    'matrix_vector_nsnb',
    'cholesky_ns3',
    'dgemm_ns3'
]

operations_list = operations_block_chol_inv
def main(filename, imgname, type, cluster):
    data = pd.read_csv(filename)

    data_g = data.groupby('diag_blocksize').mean()

    plt.figure(figsize=(8, 6))

    if type=='runtime':
        for alg_ in operations_list:
            plt.plot(
                data_g.index,
                data_g[alg_]/data_g['repetitions'],
                label=operations[alg_]['name'],
                marker='o',
                linestyle='-',
            )
        plt.ylabel('Runtime (secs)')

    if type=='performance':
        for alg_ in operations_list:
            flops = np.array([operations[alg_]['flops'](ns, 64) for ns in data_g.index])
            plt.plot(
                data_g.index,
                flops*(1e-9)/np.array(data_g[alg_]/data_g['repetitions']),
                label=operations[alg_]['name'],
                marker='o',
                linestyle='-',
            )
        peak_perf = 2764.8 # GFLOPS
        if cluster=='alex':
            peak_perf =  37400  # 37.4 TFLOPS
        plt.hlines(
            peak_perf, # peak node performance
            xmin=data_g.index[0], xmax=2048,
            color='red', alpha=0.5,
            label="peak",
            linestyle = '-',
        )
        plt.ylabel('Performance (GFLOPS/sec)')

    # Customize the plot
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)
    plt.title(f'$n_b$=64')

    plt.xscale('log', base=2)
    plt.xticks(data_g.index, [str(int(tick)) for tick in data_g.index])
    plt.xlim(data_g.index[0], 2048)
    plt.xlabel(f'$n_s$')

    plt.yscale('log', base=10)
    # plt.ylim(bottom=10e-2)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)

if __name__ == "__main__":

    for cluster in ['alex', 'fritz']:
        for type in ['runtime', 'performance']:
        
            filename = f"../jobs/{cluster}/results/operations.txt"
            imgname = f"../jobs/{cluster}/images/operations_{type}_blocked_inv.pdf"
            main(filename, imgname, type, cluster)
