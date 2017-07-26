# # Compute the prior probability of treatment using permutation 

import itertools
import statistics

import pandas

from hetio.permute import permute_pair_list

from tqdm import tqdm

# Read treatments


treatment_df = (pandas
    .read_table('../summary/indications.tsv')
    .query("rel_type == 'TREATS_CtD'")
)

# Create node to degree dictionaries
compound_to_degree = dict(treatment_df.chemical_id.value_counts())
disease_to_degree = dict(treatment_df.disease_id.value_counts())


# In[6]:

# A degree (compound_degree, disease_degree) to all potential edges with that degree
degree_to_edges = dict()

rows = list()
for (c, c_deg), (d, d_deg) in itertools.product(compound_to_degree.items(), disease_to_degree.items()):
    rows.append((c, d, c_deg, d_deg))
    degree = c_deg, d_deg
    edge = c, d
    degree_to_edges.setdefault(degree, set()).add(edge)

pair_df = pandas.DataFrame(rows, columns=['chemical_id', 'disease_id', 'compound_treats', 'disease_treats'])
pair_df = pair_df.sort_values(['chemical_id', 'disease_id'])



# Not sure whether to filter this pair_df down to just the relation pairs which are found in the training set. Will leave it as is at the moment and return to determine if this contaminates the data at all.

# In[9]:

treatments = list(zip(treatment_df.chemical_id, treatment_df.disease_id))


# In[10]:

# Burn In
pair_list, stats = permute_pair_list(treatments, multiplier=10)
pandas.DataFrame(stats)


# In[11]:

# Set the multiplier based on the burn in stats
multiplier = 3


# In[12]:

# Calculate the number of perms
n_perm = treatment_df.chemical_id.nunique() * treatment_df.disease_id.nunique()
n_perm = int(n_perm * 25)
n_perm


# In[13]:

# Initialize a dictionary of degree to empirical probability list
degree_to_probs = {x: list() for x in degree_to_edges}

# Perform n_perm permutations
for i in tqdm(range(n_perm)):
    # Permute
    pair_list, stats = permute_pair_list(pair_list, multiplier=multiplier, seed=i)

    # Update
    pair_set = set(pair_list)

    # modifies the original degree_to_probs dictionary

    for degree, probs in degree_to_probs.items():
        edges = degree_to_edges[degree]
        probs.append(len(edges & pair_set) / len(edges))


# In[14]:

rows = list()
for (c_deg, d_deg), probs in tqdm(degree_to_probs.items()):
    mean = statistics.mean(probs)
    std_error = statistics.stdev(probs) / len(probs) ** 0.5
    rows.append((c_deg, d_deg, mean, std_error))
perm_df = pandas.DataFrame(rows, columns=['compound_treats', 'disease_treats', 'prior_perm', 'prior_perm_stderr'])
perm_df = perm_df.sort_values(['compound_treats', 'disease_treats'])


# In[15]:

# Add unpermuted treatment prevalence columns
rows = list()
treatment_set = set(treatments)
for (c_deg, d_deg), edges in degree_to_edges.items():
    n_treatments = len(edges & treatment_set)
    rows.append((c_deg, d_deg, n_treatments, len(edges)))
degree_prior_df = pandas.DataFrame(rows, columns=['compound_treats', 'disease_treats', 'n_treatments', 'n_possible'])
degree_prior_df = perm_df.merge(degree_prior_df)
degree_prior_df = degree_prior_df.sort_values(['compound_treats', 'disease_treats'])


# In[17]:

degree_prior_df.to_csv('data/degree-prior.tsv', sep='\t', index=False, float_format='%.6g')

# In[18]:

obs_prior_df = pair_df.merge(perm_df)

# In[21]:

obs_prior_df.to_csv('data/observation-prior.tsv', sep='\t', index=False, float_format='%.6g')
