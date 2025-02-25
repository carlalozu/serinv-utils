import pandas as pd

def main(filename, filename_ref):

        data = pd.read_csv(filename)
        data_ref = pd.read_csv(filename_ref)

        grouped_data = data.groupby('id').agg(
                time_mean=('time', 'mean'),
                time_std=('time', 'std')
        )

        grouped_data_ref = data_ref.groupby('id').agg(
                time_ref_mean=('time', 'mean'),
                time_ref_std=('time', 'std')
        )

        df = grouped_data_ref.join(grouped_data)
        df['speedup'] = df['time_ref_mean']/df['time_mean']
        df['bandwidth'] = [int(2**i) for i in df.index]
        df.set_index('bandwidth', inplace=True) 

        # Format
        print(filename)
        print(df.map(lambda x: f"{x:.4f}").to_csv(sep='&'))

if __name__ == "__main__":

        cluster = 'alex'
        for tile in ['bb', 'b']:  # bb=block n-diag, b=banded
                for alg in ['f', 'si']:   # f=cholesky, si=selected inversion
                        arrow = 64

                        filename = f'../jobs/{cluster}/results/scpo{tile}a{alg}_{arrow}.txt'
                        filename_ref = f'../jobs/{cluster}/results/pobta{alg}_{arrow}.txt'
                
                        main(filename, filename_ref)
