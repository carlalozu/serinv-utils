import pandas as pd

def main(filename, filename_ref, output_filename, tile, alg):

        data = pd.read_csv(filename)
        data_ref = pd.read_csv(filename_ref)

        all_grouped_data = pd.DataFrame()
        for i in sorted(data['id'].unique()):
                data_ = data[data['id']==i]
                data_ref_ = data_ref[data_ref['id']==i]

                grouped_data = data_.groupby('diagonal_blocksize').mean()
                grouped_data_ref = data_ref_.mean()
                grouped_data['run'] = data_.groupby('diagonal_blocksize').max()['run']

                grouped_data['speedup'] = grouped_data_ref[f'pobta{alg}_time']/grouped_data[f'scpo{tile}a{alg}_time']
                # Format
                
                # Append to the main DataFrame
                all_grouped_data = pd.concat([all_grouped_data, grouped_data])

        # Reset the index for cleanliness
        all_grouped_data.reset_index(inplace=True)
        print(f'scpo{tile}a{alg}')
        print(all_grouped_data)
        all_grouped_data.to_csv(output_filename, index=False)

if __name__ == "__main__":

        cluster = 'alex' # alex, fritz
        for tile in ['b', 'bb']:  # bb=block n-diag, b=banded
                for alg in ['f', 'si']:

                        arrow = 64

                        filename = f'../jobs/{cluster}/results/scpo{tile}asi_{arrow}.txt'
                        if tile=='bb':
                                filename = f'../jobs/{cluster}/results/scpo{tile}asi_blocks_{arrow}.txt'
                        filename_ref = f'../jobs/{cluster}/results/pobtasi_{arrow}.txt'
                        output_filename = f'../jobs/{cluster}/results/scpo{tile}asi_{arrow}_df.txt'
                
                        main(filename, filename_ref, output_filename, tile, alg)
