import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from const import PEAK_PERFORMANCE, PLT_PARAMS, FIG_SIZE
from serinv_utils.config import PATH

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)

lims_x = {
    'fritz': (16, 4096),
    'alex': (16, 4096),
}

lims_y = {
    'fritz': (10e-2, 10e2),
    'alex': (10e-2, 10e2),
}

def main(filename, filename_ref, imgname, alg, type_, cluster):
    data_ = pd.read_csv(filename)
    data_pobtaf = pd.read_csv(filename_ref)

    plt.figure(figsize=FIG_SIZE)
    for b in sorted([513, 2049, 8193])[::-1]:
        data = data_[data_['bandwidth'] == b]

        bandwidth = int(list(data['bandwidth'])[0])
        label = f'$b$: {int(bandwidth)}'
        # label = f'Bandwidth: {int(bandwidth)}'

        grouped_data = data.groupby('diagonal_blocksize')
        data_pobtaf_ = data_pobtaf[data_pobtaf['bandwidth'] == b]
        if type_ == 'runtime':
            grouped_data = grouped_data.sum()
            data_pobtaf_ = data_pobtaf_.sum()
        else:
            grouped_data = grouped_data.mean()
            data_pobtaf_ = data_pobtaf_.mean()

        grouped_data['time'] = grouped_data[f'scpobba{alg}_time'] / \
            grouped_data['n_runs']
        
        grouped_data[f'scpobba{alg}_FLOPS'] = data.groupby('diagonal_blocksize').mean()[f'scpobba{alg}_FLOPS']

        # Plot mean times with error bars
        if type_ == 'runtime':
            print(grouped_data['n_runs'])
            time_ref = (data_pobtaf_[f'pobta{alg}_time']/data_pobtaf_['n_runs'])
            p = plt.plot(
                grouped_data.index,
                grouped_data['time'],
                # yerr=data.groupby('diagonal_blocksize')['time'].std(),
                label=label,
                marker='o',
                linestyle='-',
            )
            ylabel = 'Runtime (secs)'
            plt.hlines(time_ref, 16, 4096, color=p[0].get_color(), alpha=0.5, linestyles='dotted')
            plt.ylim(lims_y[cluster])
            plt.xlim(lims_x[cluster])
            plt.yscale('log', base=10)

        else: # performance
            flops_ref = data_pobtaf_[f'scpoba{alg}_FLOPS']
            print(grouped_data[f'scpobba{alg}_FLOPS']*10e-9)
            print("flops_ref", flops_ref*10e-9)
            plt.plot(
                grouped_data.index,
                flops_ref/grouped_data[f'scpobba{alg}_FLOPS']*100,
                label=label,
                marker='o',
                linestyle='-',
            )
            ylabel = 'Effective FLOPS (%)'
            plt.title(f'SCPOBBA{alg}')

    # Plot customization
    plt.xscale('log', base=2)
    nss = sorted(data_['diagonal_blocksize'].unique())
    plt.xticks(nss, [str(tick) for tick in nss])
    plt.xlim(nss[1], nss[-1])

    # if type_ == 'performance':
    #     plt.hlines(
    #         PEAK_PERFORMANCE[cluster],  # peak node performance
    #         xmin=nss[0], xmax=nss[-1],
    #         color='red', alpha=0.5,
    #         label="peak",
    #         linestyle='-',
    #     )

    plt.ylabel(ylabel)

    plt.xlabel('$n_s$')

    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":

    arrow = 64
    for cluster in ['fritz', 'alex']:
        for alg in ['f', 'si']:
            # f=cholesky, si=selected inversion

            for type_ in ['runtime', 'performance']:
                filename = f"{PATH}/jobs/{cluster}/results/scpobbasi_blocks_{arrow}.txt"
                filename_ref = f'{PATH}/jobs/{cluster}/results/pobtasi_{arrow}_runs.txt'
                if type_ == 'performance':
                    filename_ref = f'{PATH}/jobs/{cluster}/results/scpoba{alg}_{arrow}_df.txt'
                imgname = f"{PATH}/jobs/{cluster}/images/scpobba{alg}_blocks_{arrow}_{type_}.pdf"

                main(filename, filename_ref, imgname, alg, type_, cluster)
