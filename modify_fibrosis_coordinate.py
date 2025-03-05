

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

adata  = sc.read_h5ad('nico_celltype_annotation.h5ad')
print(adata)
#a = np.unique(ad1.obs['sample'])
#b=np.unique(ad1.obs['disease_status'])
#print(a)
#print('\n',b)

ind= np.where(adata.obs['disease_status']=='healthy')
ad= adata[ind[0],:].copy()
xy=ad.obsm['spatial']
print("healthy", xy.shape)

plt.plot(xy[:,0],xy[:,1],'b.')
plt.savefig('healthy.png')
plt.close('all')

ind= np.where(adata.obs['disease_status']=='fibrosis')[0]
adFib= adata[ind,:].copy()
shift = 5000*np.ones((len(ind),2))
print(shift)
adFib.obsm['spatial']+= shift

#xy=ad.obsm['spatial']
xy = adFib.obsm['spatial']
print('fibrosis',xy.shape,shift.shape)
plt.plot(xy[:,0],xy[:,1],'b.')
plt.savefig('fibrosis.png')
plt.close('all')


#shift coordinate of fibrosis
adata.obsm['spatial'][ind,:]=adFib.obsm['spatial']
xy=adata.obsm['spatial']
print('total',xy.shape)
plt.plot(xy[:,0],xy[:,1],'b.')
plt.savefig('total.png')
plt.close('all')

adata.write_h5ad('nico_celltype_annotation_mod.h5ad')
