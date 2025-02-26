import pandas as pd
import matplotlib.pyplot as plt
from const import PEAK_PERFORMANCE, PLT_PARAMS

plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)


def main(filename, filename_ref, imgname, alg='f', cluster='fritz'):
    data_ = pd.read_csv(filename)
    data_ref_ = pd.read_csv(filename_ref)

    fig, ax1 = plt.subplots(figsize=(8, 6))

    arrowhead = list(data_['arrowhead_blocksize'])[0]

    grouped_data = data_.groupby('n_offdiags').median()
    grouped_data_ref = data_ref_.groupby('diagonal_blocksize').median()

    # Plot mean times with error bars
    ax1.errorbar(
        grouped_data.index, 
        grouped_data[f'scpoba{alg}_time'],
        yerr=data_.groupby('n_offdiags')[f'scpoba{alg}_time'].std(),
        label=f'poba{alg}',
    )

    # Plot ref times with error bars
    ax1.errorbar(
        grouped_data_ref.index, 
        grouped_data_ref[f'pobta{alg}_time'],
        yerr=data_ref_.groupby('diagonal_blocksize')[f'pobta{alg}_time'].std(),
        label=f'pobta{alg}',
    )

    ax2 = ax1.twinx()
    # Plot mean times with error bars
    ax2.plot(
        grouped_data.index, 
        grouped_data[f'scpoba{alg}_FLOPS']*(1e-9)/grouped_data[f'scpoba{alg}_time'],
        label=f'poba{alg}',
        linestyle = '--',
    )

    # Plot ref times with error bars
    ax2.plot(
        grouped_data_ref.index, 
        grouped_data_ref[f'pobta{alg}_FLOPS']*(1e-9)/grouped_data_ref[f'pobta{alg}_time'],
        label=f'pobta{alg}',
        linestyle = '--',
    )

    ax2.hlines(
        PEAK_PERFORMANCE[cluster], # peak node performance
        xmin=grouped_data.index[0],
        xmax=grouped_data.index[-1],
        color='red',
        alpha=0.5,
        label="peak",
        linestyle = '--',
    )

    plt.xticks(grouped_data.index, [str(tick) for tick in grouped_data.index])
    plt.xlim(grouped_data.index[0], 2048)
    ax1.set_xlabel(f'$n_s$')

    # ax1.set_ylim(10e-2, 10e1)
    ax1.set_yscale('log', base=10)
    ax1.set_ylabel('Runtime (sec)')
    
    # ax2.set_ylim(10e-2, 10e4)
    ax2.set_yscale('log', base=10)
    ax2.set_ylabel('Performance (GFLOPS/sec)')

    plt.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend(loc='upper left', frameon=True, shadow=False)
    ax2.legend(loc='lower right', frameon=True, shadow=False)


    # Save as pdf
    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":
    for cluster in ['fritz', 'alex']:
        arrowhead = '64'
        for alg in ['f','si']: # goes with diagonal_blocksize
            # alg = 'scpobasi' # goes with n_offdiags

            filename = f"../jobs/{cluster}/results/scpobasi_{arrowhead}.txt"
            filename_ref = f"../jobs/{cluster}/results/pobtasi_{arrowhead}.txt"
            imgname = f"../jobs/{cluster}/images/scpoba{alg}_{arrowhead}.pdf"
            main(filename, filename_ref, imgname, alg, cluster)
