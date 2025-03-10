import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
from const import PEAK_PERFORMANCE, PLT_PARAMS
from serinv_utils.scaling.flops.const import OPERATIONS_FLOPS

# plt.style.use("seaborn-v0_8-colorblind")
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

lims = {
    'block_chol_alex': (16, 4096),
    'block_inv_alex': (16, 4096),
    'block_chol_fritz': (16, 2048),
    'block_inv_fritz': (16, 2048),
    'banded_alex': (64, 2048),
    'banded_fritz': (64, 2048),
}


def main(filename, imgname, type_, routine, cluster):
    data = pd.read_csv(filename)

    data_g = data.groupby('diag_blocksize').sum()

    plt.figure(figsize=(8, 6))

    if type_ == 'runtime':
        for alg_ in operations_list[routine]:
            plt.plot(
                data_g.index,
                data_g[alg_]/data_g['repetitions'],
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
                flops*(1e-9)/np.array(data_g[alg_]/data_g['repetitions']),
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
    plt.xlim(*lims[f'{routine}_{cluster}'])
    plt.xlabel('$n_s$')

    plt.yscale('log', base=10)
    # plt.ylim(bottom=10e-2)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":

    for cluster in ['alex', 'fritz']:
        for type_ in ['runtime', 'performance']:
            for routine in operations_list:
                filename = f"../jobs/{cluster}/results/operations.txt"
                imgname = f"../jobs/{cluster}/images/operations_{routine}_{type_}.pdf"
                main(filename, imgname, type_, routine, cluster)
