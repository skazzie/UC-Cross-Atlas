import anndata, numpy as np, pandas as pd
a = anndata.read_h5ad('data/atlases/garrido.h5ad')
c = a.layers['counts']
n_genes  = np.asarray((c > 0).sum(axis=1)).ravel()
n_counts = np.asarray(c.sum(axis=1)).ravel()
cov = pd.DataFrame(index=a.obs_names); cov.index.name = 'cell_id'
cov['log_n_genes']  = np.log10(n_genes)
cov['log_n_counts'] = np.log10(n_counts)
donor = pd.get_dummies(a.obs['donor_id'].astype(str), prefix='donor', drop_first=True).astype(int)
donor.index = a.obs_names
cov = pd.concat([cov, donor], axis=1)
spd = a.obs.groupby('donor_id')['sample'].nunique().max()
if spd > 1:
    s = pd.get_dummies(a.obs['sample'].astype(str), prefix='sample', drop_first=True).astype(int)
    s.index = a.obs_names; cov = pd.concat([cov, s], axis=1)
    print(f'included sample one-hot (max {spd}/donor)')
else:
    print('skipped sample one-hot (1 sample/donor, collinear)')
cov.to_csv('data/atlases/garrido_covariates.tsv', sep='\t')
print('wrote data/atlases/garrido_covariates.tsv', cov.shape)
print('cols:', cov.columns.tolist()[:6], '...')
