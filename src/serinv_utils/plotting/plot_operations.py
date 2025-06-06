import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from const import PEAK_PERFORMANCE, PLT_PARAMS, FIG_SIZE
from serinv_utils.flops.const import OPERATIONS_FLOPS
PATH = ".."

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)

operations_list = {
    "banded": [
        'dgemm_ns3',
        'scale_ns',
        'dot_prod_ns',
        'matrix_vector_nsns',
        'matrix_vector_nsnb',
        'cholesky_ns3',
        'triang_solve_ns3',
    ],
    "block_chol": [
        'dgemm_ns3',
        'dgemm_ns2nb',
        'dgemm_nb2ns',
        'triang_solve_ns3',
        'triang_solve_ns2nb',
        'cholesky_ns3',
    ],
    "block_inv": [
        'dgemm_ns3',
        'dgemm_ns2nb',
        'dgemm_nb2ns',
        'triang_solve_ns3',
    ],
}

lims_x = {
    "alex": {
        'block_chol': (16, 4096),
        'block_inv': (16, 4096),
        'banded': (64, 2048),
    },
    "fritz": {
        'block_chol': (16, 4096),
        'block_inv': (16, 4096),
        'banded': (64, 2048),
    }
}

lims_y = {
    'alex' : {
        'block_chol': (10e-6, 10e-3),
        'block_inv': {"auto":True},
        'banded': {"auto":True},
    },
    'fritz': {
        'block_chol': {"bottom":10e-7, "top":10e-1},
        'block_inv': {"auto":True},
        'banded': {"auto":True},
    }
}


def main(filename, imgname, type_, routine, cluster):
    
    data = pd.read_csv(filename)
    data_g = data.groupby('diag_blocksize').sum()

    plt.figure(figsize=FIG_SIZE)

    if type_ == 'runtime':
        for alg_ in operations_list[routine]:
            plt.plot(
                data_g.index,
                data_g[alg_]/data_g['n_runs'],
                label=OPERATIONS_FLOPS[alg_]['name'],
                marker='o',
                linestyle='-',
            )
        plt.ylabel('Runtime (secs)')

    if type_ == 'performance':
        for alg_ in operations_list[routine]:
            flops = np.array([OPERATIONS_FLOPS[alg_]['flops'](ns, 64)
                             for ns in data_g.index])
            plt.plot(
                data_g.index,
                flops*(1e-9)/np.array(data_g[alg_]/data_g['n_runs']),
                label=OPERATIONS_FLOPS[alg_]['name'],
                marker='o',
                linestyle='-',
            )

        plt.hlines(
            PEAK_PERFORMANCE[cluster],  # peak node performance
            xmin=data_g.index[0], xmax=data_g.index[-1],
            color='red', alpha=0.5,
            label="peak",
            linestyle='-',
        )
        plt.ylabel('Performance (GFLOPS/sec)')

    # Customize the plot
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    if 'block' in routine:
        plt.xscale('log', base=2)
    plt.xticks(data_g.index, [str(int(tick)) for tick in data_g.index])
    plt.xlim(lims_x[cluster][routine])
    plt.ylim(**lims_y[cluster][routine])
    plt.xlabel('$n_s$')

    plt.yscale('log', base=10)
    # plt.ylim(bottom=10e-2)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":

    for cluster in ['alex', 'fritz']:
        for type_ in ['runtime', 'performance']:
            for routine in operations_list:
                filename = f"{PATH}/jobs/{cluster}/results/operations.csv"
                imgname = f"{PATH}/jobs/{cluster}/images/operations_{routine}_{type_}.pdf"
                main(filename, imgname, type_, routine, cluster)
