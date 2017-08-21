
# coding: utf-8

# # New extract and matrixify in one notebook
# 
# Goal to reproduce the files needed.
# 
# Skip 3-extract
# skip parts of 4-matrixify
# 
# but produce the same output files

# In[1]:

import pandas as pd
import numpy as np

import os
import bz2
import itertools

import sys

sys.path.append("..")


# In[2]:

from src.extractor import MatrixFormattedGraph


# ---

# ## Get the matricies for all permutations

# In[3]:

def get_matrix(folder, index):
    perm_val = "" if index == 0 else "_perm-{}".format(index)
    nodes_fname = "{}/hetnet{}_nodes.csv".format(folder, perm_val)
    edges_fname = "{}/hetnet{}_edges.csv".format(folder, perm_val)

    return MatrixFormattedGraph(
        nodes_fname, edges_fname,
        start_kind="Compound", end_kind="Disease",
        max_length=4
    )


# In[4]:

folder = "../../integrate/data/import_csvs"

matricies = [
    get_matrix(folder, index) for index in range(6)
]


# ---

# ## Extract out the DWPCs

# In[5]:

dwpcs = [
    (mg
        .extract_dwpc(
            start_nodes="Compound", end_nodes="Disease",
            n_jobs=32
        )
        .rename(columns={"compound_id": "chemical_id"})
    )
    
    for mg in matricies
]


# ---

# ## Subset data

# In[6]:

partitions = (pd
    .read_csv("data/partitions.tsv", sep='\t')
    .assign(idx = lambda df: df["hetnet"].str[-1].astype(np.int64))
)


# In[7]:

dwpc_spread_df = pd.concat([
    dwpc_df.merge(
        partitions.query("idx == @val")[["hetnet", "chemical_id", "disease_id"]],
        how="right", on=["chemical_id", "disease_id"]
    )
    
    for (val, dwpc_df) in enumerate(dwpcs)
])


# In[8]:

dwpc_spread_df = dwpc_spread_df.dropna(axis=1, how="any")


# In[9]:

dwpc_spread_df.shape


# In[10]:

dwpc_spread_df.head()


# ---

# ## Filter out features which are all zero

# In[ ]:

bad_features = []
for col_name in dwpc_spread_df.columns:
    if dwpc_spread_df[col_name].dtype == "float64":
        if np.isclose(dwpc_spread_df[col_name], 0, atol=1e-6):
            bad_features.append(col_name)


# In[ ]:

dwpc_spread_df = dwpc_spread_df.drop(bad_features, axis=1)


# ---

# In[11]:

dwpc_spread_df["hetnet"].value_counts()


# ## Write to file

# In[12]:

path = 'data/matrix/dwpc.tsv.bz2'
with bz2.open(path, 'wt') as wf:
    dwpc_spread_df.to_csv(wf, index=False, sep='\t', float_format='%.5g')


# ---

# ## Calculate degree features

# In[13]:

url = "../../integrate/data/summary/metaedge-styles.tsv"

metaedge_style_df = pd.read_table(url)
metaedge_to_abbreviation = dict(zip(metaedge_style_df.metaedge, metaedge_style_df.abbreviation))

url = "../../integrate/data/summary/degrees.xlsx"

disease_degree_df = (pd
    .read_excel(url, sheetname='Disease')
    .rename(columns={'node_id': 'disease_id'})
    .drop('node_name', axis='columns')
    .rename(columns=metaedge_to_abbreviation)
)

compound_degree_df = (pd
    .read_excel(url, sheetname='Compound')
    .rename(columns={'node_id': 'chemical_id'})
    .drop('node_name', axis='columns')
    .rename(columns=metaedge_to_abbreviation)
)


# In[14]:

compound_degree_df.head(2)


# In[15]:

disease_degree_df.head(2)


# In[16]:

compound_degree_df.to_csv('data/matrix/compound_degree.tsv', index=False, sep='\t')
disease_degree_df.to_csv('data/matrix/disease_degree.tsv', index=False, sep='\t')


# ---

# ## Compute prior dataset

# In[17]:

# Read compound and disease degrees
compound_df = pd.read_table('../summary/compounds.tsv')
disease_df = pd.read_table('../summary/diseases.tsv')

total_pairs = len(compound_df) * len(disease_df)

nonzero_prior_pairs = sum(compound_df.treats > 0) * sum(disease_df.treats > 0)
total_pairs, nonzero_prior_pairs


# In[18]:

rows = list(itertools.product(compound_df.chemical_id, disease_df.disease_id))

prior_df = (pd
    .DataFrame(rows, columns=['chemical_id', 'disease_id'])
    .merge(
        pd.read_table('../prior/data/observation-prior.tsv')[['chemical_id', 'disease_id', 'prior_perm']],
        how='left'
    )
    .fillna(0)
    .rename(columns={'prior_perm': 'prior_prob'})
)

prior_df.head(2)


# In[19]:

sum(prior_df.prior_prob)


# In[20]:

(prior_df.prior_prob > 0).value_counts(True)


# In[21]:

prior_df.to_csv('data/matrix/prior.tsv', index=False, sep='\t', float_format='%.5g')


# ---

# ## Create a single matrix-like dataframe

# In[22]:

matrix_df = (partitions
    .drop("idx", axis=1)
    .merge(disease_df.iloc[:, :2])
    .merge(compound_df.iloc[:, :2])
    .merge(prior_df)
    .merge(compound_degree_df)
    .merge(disease_degree_df)
    .merge(dwpc_spread_df)
)


# In[23]:

matrix_df.head(2)


# In[24]:

df_creators = [
    {'feature_type': 'prior', 'feature': ['prior_prob']},
    {'feature_type': 'degree', 'feature': compound_degree_df.columns[1:]},
    {'feature_type': 'degree', 'feature': disease_degree_df.columns[1:]},
    
    # this line is super fragile
    # drops the "hetnet" column and the two identifier columns
    {
        'feature_type': 'dwpc',
        'feature': dwpc_spread_df.drop([
            "chemical_id", "disease_id", "hetnet"
        ], axis=1).columns
    },
]
feature_df = pd.concat(map(pd.DataFrame, df_creators))


# In[25]:

unperm_name = 'rephetio-v2.0'


# In[26]:

unperm_matrix_df = (matrix_df
    .query("hetnet == @unperm_name")
    .drop('hetnet', axis=1)
)


# In[27]:

feature_df['unperm_mean'] = list(
    unperm_matrix_df[feature_df.feature].mean()
)


# In[28]:

feature_df['unperm_sd'] = list(
    unperm_matrix_df[feature_df.feature].std()
)
feature_df.head(2)


# In[29]:

feature_df.to_csv('data/matrix/feature-type.tsv', index=False, sep='\t', float_format='%.5g')

path = 'data/matrix/features.tsv.bz2'
with bz2.open(path, 'wt') as wf:
    matrix_df.to_csv(wf, index=False, sep='\t', float_format='%.5g')


# In[30]:

# Save hetnet specific feature files
directory = os.path.join('data', 'matrix', unperm_name)
if not os.path.exists(directory):
    os.mkdir(directory)
path = os.path.join(directory, 'features.tsv.bz2')
with bz2.open(path, 'wt') as wf:
    unperm_matrix_df.to_csv(wf, index=False, sep='\t', float_format='%.5g')

