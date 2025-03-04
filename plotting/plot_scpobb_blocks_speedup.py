import pandas as pd
import matplotlib.pyplot as plt
from const import PLT_PARAMS

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)

lims = {
    'alex': (16, 4096),
    'fritz': (16, 2048),
}


def main(filename, filename_ref, imgname, cluster, alg):
    data_scpobbaf = pd.read_csv(filename)
    data_pobtaf = pd.read_csv(filename_ref)

    plt.figure(figsize=(8, 6))
    for i in sorted(data_scpobbaf['id'].unique())[::-1]:
        data = data_scpobbaf[data_scpobbaf['id'] == i]
        data_pobtaf_ = data_pobtaf[data_pobtaf['id'] == i].sum()

        time_ref = (data_pobtaf_[f'pobta{alg}_time']/data_pobtaf_['n_runs'])
        print(time_ref)

        bandwidth = list(data['bandwidth'])[0]

        grouped_data = data.groupby(
            'diagonal_blocksize').sum().reset_index()
        grouped_data['time'] = grouped_data[f'scpobba{alg}_time'] / \
            grouped_data['n_runs']
        

        # Scaling
        grouped_data[f'scpobba{alg}_speedup'] = time_ref / \
            grouped_data['time']

        # Plot mean times with error bars
        plt.errorbar(
            grouped_data['diagonal_blocksize'],
            grouped_data[f'scpobba{alg}_speedup'],
            label=f'$b$={int(bandwidth)}',
            marker='o',
            linestyle='-',
        )

    plt.axhline(1, color='black', linestyle='-')

    plt.xlabel('$n_s$')
    plt.xscale('log', base=2)

    nss = sorted(data_scpobbaf['diagonal_blocksize'].unique())
    plt.xticks(nss, [str(tick) for tick in nss])
    plt.xlim(*lims[cluster])

    # plt.yscale('log', base=10)
    plt.ylabel('Speedup')

    # Customize the plot
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    # Save the plot
    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":

    arrow = 64
    for cluster in ['alex']:  # , 'alex']:  # , 'alex']:
        for alg in ['f', 'si']:
            # f=cholesky, si=selected inversion

            filename = f"../jobs/{cluster}/results/scpobbasi_blocks_{arrow}.txt"
            filename_ref = f'../jobs/{cluster}/results/pobtasi_{arrow}_runs.txt'
            imgname = f"../jobs/{cluster}/images/scpobba{alg}_blocks_{arrow}_speedup.pdf"
            main(filename, filename_ref, imgname, cluster, alg)
