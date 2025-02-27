import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from const import PLT_PARAMS
import sys

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)

lims = {
    'alex':(16,4096),
    'fritz':(16,2048),
}

def main(filename, filename_ref, imgname, cluster):
    data_scpobbaf = pd.read_csv(filename)
    data_pobtaf = pd.read_csv(filename_ref)

    plt.figure(figsize=(8, 6))
    for i in sorted(data_scpobbaf['id'].unique())[::-1]:
        data = data_scpobbaf[data_scpobbaf['id']==i]
        time_ref = data_pobtaf[data_pobtaf['id']==i][f'pobta{alg}_time'].mean()

        bandwidth = list(data['bandwidth'])[0]
        arrowhead = list(data['arrowhead_blocksize'])[0]
        matrix_size = list(data['n'])[0]

        grouped_data = data.groupby('diagonal_blocksize').mean().reset_index()

        # Scaling
        grouped_data[f'scpobba{alg}_speedup'] = time_ref/grouped_data[f'scpobba{alg}_time']

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
    for cluster in ['alex']: #, 'alex']:
        for alg in ['f', 'si']:
            # f=cholesky, si=selected inversion
    
            filename = f"../jobs/{cluster}/results/scpobbasi_blocks_{arrow}.txt"
            filename_ref = f'../jobs/{cluster}/results/pobtasi_{arrow}.txt'
            imgname = f"../jobs/{cluster}/images/scpobba{alg}_blocks_{arrow}_speedup.pdf"
            main(filename, filename_ref, imgname, cluster)
