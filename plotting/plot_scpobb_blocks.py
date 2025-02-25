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

def main(filename, imgname, alg, type, cluster):
    data_ = pd.read_csv(filename)

    plt.figure(figsize=(8, 6))
    for i in sorted(data_['id'].unique())[::-1]:
        data = data_[data_['id']==i]

        bandwidth = int(list(data['bandwidth'])[0])
        arrowhead = int(list(data['arrowhead_blocksize'])[0])
        matrix_size = int(list(data['n'])[0])

        grouped_data = data.groupby('diagonal_blocksize').median()


        label = f'$b$: {int(bandwidth)}'
        # label = f'Bandwidth: {int(bandwidth)}'
        # Plot mean times with error bars
        if type=='runtime':
            plt.plot(
                grouped_data.index, 
                grouped_data[f'scpobba{alg}_time'], 
                # yerr=data.groupby('diagonal_blocksize')['time'].std(),
                label=label,
                marker='o',
                linestyle='-',
            )
            ylabel = 'Runtime (secs)'
        else:
            plt.plot(
                grouped_data.index, 
                grouped_data[f'scpobba{alg}_FLOPS']/grouped_data[f'scpobba{alg}_time']*1e-9, 
                label=label,
                marker='o',
                linestyle='-',
            )
            ylabel = 'Performance (GFLOPS/sec)'


    # Plot customization
    plt.xscale('log', base=2)
    nss = sorted(data_['diagonal_blocksize'].unique())
    plt.xticks(nss, [str(tick) for tick in nss])
    plt.xlim(nss[1], nss[-1])

    if type=='performance':
        peak_perf = 2764.8 # GFLOPS
        if cluster=='alex':
            peak_perf =  37400  # 37.4 TFLOPS
        plt.hlines(
            peak_perf, # peak node performance
            xmin=nss[0], xmax=nss[-1],
            color='red', alpha=0.5,
            label="peak",
            linestyle = '-',
        )

    plt.yscale('log', base=10)
    plt.ylabel(ylabel)

    plt.xlabel(f'$n_s$')
    
    title = f'matrix size={matrix_size}, $n_b$={int(arrowhead)}'
    plt.title(title)

    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)

if __name__ == "__main__":

    for cluster in ['fritz', 'alex']:
        for alg in ['f', 'si']:
            # f=cholesky, si=selected inversion
            for arrow in [64]: #, 128]:
                for type in ['runtime', 'performance']:
                    filename = f"../../jobs/{cluster}/results/scpobbasi_blocks_{arrow}.txt"
                    imgname = f"../../jobs/{cluster}/images/scpobba{alg}_blocks_{arrow}_{type}.pdf"

                    main(filename, imgname, alg, type, cluster)