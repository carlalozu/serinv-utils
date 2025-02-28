from serinv_utils.scaling.flops.const import ALG_OPERATION_COUNT, BENCHMARKED_OPS
from serinv_utils.scaling.flops.scpobbaf_flops import scpobbaf_flops
from serinv_utils.scaling.flops.scpobbasi_flops import scpobbasi_flops

from serinv_utils.scaling.storage.parameters import calculate_parameters_n_diagonal
from const import PLT_PARAMS
import pandas as pd
import matplotlib.pyplot as plt


plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)


lims = {
    'POBBAF_alex': (16, 4096),
    'POBBAF_fritz': (16, 2048),
    'POBBASI_alex': (16, 4096),
    'POBBASI_fritz': (16, 2048),
}

NB: int = 64
M: int = 2**16+NB


def main(filename, filename_alg, imgname, routine, cluster):
    data_ops = pd.read_csv(filename)
    data_rou = pd.read_csv(filename_alg)

    operations_df = data_ops.groupby('diag_blocksize').sum().reset_index()

    plt.figure(figsize=(8, 6))

    for b in sorted(data_rou['bandwidth'].unique())[::-1]:
        data_rou_ = data_rou[data_rou['bandwidth'] == b]

        times = []

        data_rou_b = data_rou_.groupby('n_offdiags').mean()
        for n in sorted(data_rou_['n_offdiags'].unique()):
            p = calculate_parameters_n_diagonal(M, b, NB, n)['parameters']

            print(
                'b', b,
                'n_t', p['n_t'],
                'diagonal_blocksize', p['diagonal_blocksize'],
                'arrowhead_blocksize', p['arrowhead_blocksize'],
                'n_offdiags', int(p['n_offdiags']),
            )
            if routine == 'POBBASI':
                count_alg = scpobbasi_flops
            else:
                count_alg = scpobbaf_flops

            _, counts = count_alg(
                p['n_t'],
                p['diagonal_blocksize'],
                p['arrowhead_blocksize'],
                int(p['n_offdiags']),
            )

            if p['diagonal_blocksize'] < 16:
                times.append(0)
                continue

            routine_time = 0
            for alg_ in ALG_OPERATION_COUNT[routine]:

                counts_ = counts[alg_]
                # counts_ = ALG_OPERATION_COUNT[routine][alg_](p['n_t'], p['n_offdiags'])

                if alg_ in BENCHMARKED_OPS:

                    alg_df = operations_df.loc[operations_df['diag_blocksize']
                                               == p['diagonal_blocksize']].copy()

                    alg_df['time'] = alg_df[alg_]/alg_df['repetitions']

                    print('adding', alg_, 'with', counts_)
                    routine_time += (alg_df['time']*counts_).item()

                else:
                    alg_ns_ = alg_
                    if 'b' in alg_:
                        alg_ns_ = alg_ns_.replace('b', 's')
                        print('replace', alg_, alg_ns_)

                    if alg_ns_ in BENCHMARKED_OPS:
                        print('adding', alg_ns_, 'with', NB)

                        alg_df = operations_df.loc[operations_df['diag_blocksize'] == NB].copy(
                        )
                        alg_df['time'] = alg_df[alg_ns_]/alg_df['repetitions']
                        t = alg_df['time']*counts_
                        routine_time += t.item()
                    else:
                        print('skipping', alg_)
            times.append(routine_time)

        data_rou_b['time'] = times
        # Plot runtime of routine
        plt.plot(
            data_rou_b['diagonal_blocksize'],
            data_rou_b['time'],
            label=f'{int(b)}',
            marker='o',
        )

    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    plt.xscale('log', base=2)
    plt.xticks(operations_df['diag_blocksize'], [str(int(tick))
               for tick in operations_df['diag_blocksize']])
    plt.xlim(*lims[f'{routine}_{cluster}'])

    plt.ylabel('Runtime (secs)')
    plt.xlabel('$n_s$')

    plt.yscale('log', base=10)
    # plt.ylim(bottom=10e-2)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":

    for cluster in ['fritz', 'alex']:  # , 'fritz', alex
        # ['POBAF', 'POBASI', 'POBBAF', 'POBBASI',]:
        for routine in ['POBBASI']:
            filename = f"../jobs/{cluster}/results/operations.txt"
            filename_rou = f"../jobs/{cluster}/results/sc{routine.lower()}_64_df.txt"
            imgname = f"../jobs/{cluster}/images/operations_{routine}_theory.pdf"
            main(filename, filename_rou, imgname, routine, cluster)
