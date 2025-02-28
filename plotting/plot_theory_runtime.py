from serinv_utils.scaling.flops.const import OPERATIONS_FLOPS, ALG_OPERATION_COUNT, BENCHMARKED_OPS
from serinv_utils.scaling.storage.parameters import calculate_parameters_tri_diagonal
from const import PLT_PARAMS
import pandas as pd
import matplotlib.pyplot as plt


# plt.style.use("seaborn-v0_8-colorblind")
plt.rcParams.update(PLT_PARAMS)


lims = {
    'POBBAF_alex': (16, 4096),
    'POBBAF_fritz': (16, 2048),
    'POBBASI_alex': (16, 4096),
    'POBBASI_fritz': (16, 2048),
    'POBAF_alex': (256, 2048),
    'POBAF_fritz': (256, 2048),
    'POBASI_alex': (256, 2048),
    'POBASI_fritz': (256, 2048),
}

NB: int = 64
M: int = 65600


def main(filename, imgname, routine, cluster):
    data = pd.read_csv(filename)

    operations_df = data.groupby('diag_blocksize').sum().reset_index()

    plt.figure(figsize=(8, 6))

    routine_time = pd.DataFrame(0, index=operations_df.index, columns=['time'])
    for alg_ in ALG_OPERATION_COUNT[routine]:
        alg_f = alg_
        if '_1' in alg_f:
            alg_f = alg_f.replace('_1', '')
            print('replace', alg_, 'with', alg_f)
        if alg_f in BENCHMARKED_OPS:
            print('adding', alg_)
            alg_df = operations_df[[
                alg_f, 'repetitions', 'diag_blocksize']].copy()
            alg_df['time'] = alg_df[alg_f]/alg_df['repetitions']
            alg_df['counts'] = ALG_OPERATION_COUNT[routine][alg_](M-NB, 1)
            routine_time['time'] += alg_df['time']*alg_df['counts']

            plt.plot(
                alg_df['diag_blocksize'],
                alg_df['time'],
                label=OPERATIONS_FLOPS[alg_]['name'],
                linestyle='dashed',
                alpha=0.7,
            )

        else:
            alg_ns_ = alg_
            if 'b' in alg_:
                alg_ns_ = alg_ns_.replace('b', 's')
                print('replace', alg_, alg_ns_)

            if alg_ == 'ger_nb2':
                alg_ns_ = 'matrix_vector_nsnb'
                print('replace', alg_, alg_ns_)

            if alg_ns_ in BENCHMARKED_OPS:
                print('adding', alg_ns_, 'with', NB)
                alg_df = operations_df[[
                    alg_ns_, 'repetitions', 'diag_blocksize']].copy()

                alg_df = alg_df.loc[alg_df['diag_blocksize'] == NB]
                alg_df['time'] = alg_df[alg_ns_]/alg_df['repetitions']
                alg_df['counts'] = ALG_OPERATION_COUNT[routine][alg_](M-NB, 1)
                t = alg_df['time']*alg_df['counts']
                routine_time['time'] += t.item()
            else:
                print('skipping', alg_)

    # Plot runtime of routine
    print(routine_time.head())
    plt.plot(
        operations_df['diag_blocksize'],
        routine_time,
        label=f'{routine} theory',
        marker='o',
    )

    routine_ref = {
        'POBASI': 'POBTASI',
        'POBAF': 'POBTAF',

    }
    routine_time = pd.DataFrame(0, index=operations_df.index, columns=['time'])
    for alg_ in ALG_OPERATION_COUNT[routine_ref[routine]]:
        nts = []
        for ns in operations_df['diag_blocksize']:
            p = calculate_parameters_tri_diagonal(M, ns*2+1, NB)
            nts.append(p['parameters']['n_t'])

        if alg_ in BENCHMARKED_OPS:
            print('adding', alg_)
            alg_df = operations_df[[
                alg_, 'repetitions', 'diag_blocksize']].copy()
            alg_df['time'] = alg_df[alg_]/alg_df['repetitions']
            alg_df['counts'] = [
                ALG_OPERATION_COUNT[routine_ref[routine]][alg_](nt, 1) for nt in nts]
            routine_time['time'] += alg_df['time']*alg_df['counts']

            plt.plot(
                alg_df['diag_blocksize'],
                alg_df['time'],
                label=OPERATIONS_FLOPS[alg_]['name'],
                linestyle='dashdot',
                alpha=0.7,
            )

        else:
            alg_ns_ = alg_
            if 'b' in alg_:
                alg_ns_ = alg_ns_.replace('b', 's')
                print('replace', alg_, alg_ns_)

            if alg_ns_ in BENCHMARKED_OPS:
                print('adding', alg_ns_, 'with ', NB)
                alg_df = operations_df[[
                    alg_ns_, 'repetitions', 'diag_blocksize']].copy()
                alg_df['counts'] = [
                    ALG_OPERATION_COUNT[routine_ref[routine]][alg_](nt, 1) for nt in nts]

                alg_df = alg_df.loc[alg_df['diag_blocksize'] == NB]
                alg_df['time'] = alg_df[alg_ns_]/alg_df['repetitions']
                t = alg_df['time']*alg_df['counts']
                routine_time['time'] += t.item()
            else:
                print('skipping', alg_)

    # Plot runtime of routine
    print(routine_time.head())
    plt.plot(
        operations_df['diag_blocksize'],
        routine_time,
        label=f'{routine_ref[routine]} theory',
        marker='o',
        linestyle='-',
    )

    plt.ylabel('Runtime (secs)')

    # Customize the plot
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(loc='best', frameon=True, shadow=True)

    if 'block' in routine:
        plt.xscale('log', base=2)
    plt.xticks(operations_df['diag_blocksize'], [str(int(tick))
               for tick in operations_df['diag_blocksize']])
    plt.xlim(*lims[f'{routine}_{cluster}'])
    plt.xlabel('$n_s$')

    plt.yscale('log', base=10)
    # plt.ylim(bottom=10e-2)

    plt.tight_layout()
    plt.savefig(imgname, dpi=300)


if __name__ == "__main__":

    for cluster in ['fritz', 'alex']:  # , 'fritz', alex
        # ['POBAF', 'POBASI', 'POBBAF', 'POBBASI',]:
        for routine in ['POBAF']:
            filename = f"../jobs/{cluster}/results/operations.txt"
            imgname = f"../jobs/{cluster}/images/operations_{routine}_theory.pdf"
            main(filename, imgname, routine, cluster)
