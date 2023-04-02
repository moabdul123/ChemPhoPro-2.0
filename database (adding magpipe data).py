import pandas as pd
import sqlite3

conn = sqlite3.connect('chepro.db')

# Convert cell line tables to pandas dataframes
mcf7_df = pd.read_sql_query("SELECT * FROM 'MCF-7'", conn)
ntera2_df = pd.read_sql_query("SELECT * FROM 'NTERA-2 clone D1'", conn)
hl60_df = pd.read_sql_query("SELECT * FROM 'HL-60'", conn)

# Import the TSV files into pandas dataframes
mcf7_tsv_df = pd.read_csv('MCF7.tsv', sep='\t')
ntera2_tsv_df = pd.read_csv('NTERA2.tsv', sep='\t')
hl60_tsv_df = pd.read_csv('HL60.tsv', sep='\t')

# merge the dataframes via perturbagen and substrates (how = left means that the left df will be merged with the right df)
mcf7_df = pd.merge(mcf7_df, mcf7_tsv_df, how='left',
                   on=['perturbagen', 'substrate'])
ntera2_df = pd.merge(ntera2_df, ntera2_tsv_df, how='left',
                     on=['perturbagen', 'substrate'])
hl60_df = pd.merge(hl60_df, hl60_tsv_df, how='left',
                   on=['perturbagen', 'substrate'])

# Import the merged dataframes back to the database
mcf7_df.to_sql('MCF-7', conn, if_exists='replace', index=False)
ntera2_df.to_sql('NTERA-2 clone D1', conn, if_exists='replace', index=False)
hl60_df.to_sql('HL-60', conn, if_exists='replace', index=False)


conn.close()
