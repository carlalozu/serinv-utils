import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update({
    'axes.labelsize': 16,  # X and Y label font size
    'axes.titlesize': 16,  # Title font size
    'xtick.labelsize': 14,  # X-tick labels font size
    'ytick.labelsize': 14,  # Y-tick labels font size
    'legend.fontsize': 14,   # Legend font size
    'lines.linewidth': 3,
})


def main(filename, filename_ref, imgname):
    data_scpobbaf = pd.read_csv(filename)
    data_pobtaf = pd.read_csv(filename_ref)

    plt.figure(figsize=(8, 6))
    for i in data_scpobbaf['id'].unique()[::-1]:
        data = data_scpobbaf[data_scpobbaf['id']==i]
        time_ref = data_pobtaf[data_pobtaf['id']==i]['time'].mean()

        bandwidth = list(data['bandwidth'])[0]
        arrowhead = list(data['arrowhead_blocksize'])[0]
        matrix_size = list(data['n'])[0] - arrowhead

        grouped_data = data.groupby('diagonal_blocksize').agg({
            'time': 'mean',
            'numpy_time': 'mean',
        }).reset_index()
        # Scaling
        grouped_data['time'] = time_ref/grouped_data['time']

        # Plot mean times with error bars
        plt.errorbar(
            grouped_data['diagonal_blocksize'], 
            grouped_data['time'], 
            label=f'$b$={int(bandwidth)}',
            marker='o',
            linestyle='-',
        )

    plt.axhline(1, color='black', linestyle='-')

    plt.xlabel('$n_s$')
    plt.xscale('log', base=2)

    nss = sorted(data_scpobbaf['diagonal_blocksize'].unique())
    plt.xticks(nss, [str(tick) for tick in nss])
    plt.xlim(nss[1], nss[-1])

    # plt.yscale('log', base=10)
    plt.ylabel('Speedup')

    plt.title(f'matrix size={int(matrix_size)}, $n_b$={int(arrowhead)}')

    # Customize the plot
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    # Save the plot
    plt.tight_layout()
    plt.savefig(imgname, dpi=300)

if __name__ == "__main__":
    for cluster in ['fritz', 'alex']:
        for alg in ['f', 'si']:
            # f=cholesky, si=selected inversion
            for arrow in [64]: #, 128]: #, 512]:
        
                filename = f"../jobs/{cluster}/results/scpobba{alg}_blocks_{arrow}.txt"
                filename_ref = f'../jobs/{cluster}/results/pobta{alg}_{arrow}.txt'
                imgname = f"../jobs/{cluster}/images/scpobba{alg}_blocks_{arrow}_speedup.pdf"
                main(filename, filename_ref, imgname)
