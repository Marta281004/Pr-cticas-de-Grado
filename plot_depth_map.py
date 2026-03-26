#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# plot_depth_map

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse 

#definimos los args de la termianl
parser = argparse.ArgumentParser(description="Plot depht map")
parser.add_argument("--filename", required=True, type=str, help="Este el fichero de datos") #Selección de archivo
args = parser.parse_args()

#para leer el archivo
data = np.loadtxt(args.filename, skiprows=1) # pq los datos empiezan a partir de la fila dos, ya que en la 1 están las cabeceras
cell_id = data[:,0]
N = data[:,1]
ra = data[:,2] #cento del cuadrado
dec = data[:,3] #centro del cuadrado
maglim = data[:,4]
err_maglim = data[:,5]

#valores únicos de la rejilla
ra_unique = np.unique(ra)
dec_unique = np.unique(dec)

depth_map = np.full((len(dec_unique), len(ra_unique)), np.nan)

for i in range(len(ra)):
    ra_idx = np.where(ra_unique == ra[i])[0][0]
    dec_idx = np.where(dec_unique == dec[i])[0][0]
    depth_map[dec_idx, ra_idx] = maglim[i]
    
#plot
plt.figure()
plt.imshow(depth_map, origin="lower", extent=[ra_unique.min(), ra_unique.max(), dec_unique.min(), dec_unique.max()], aspect="auto")
plt.colorbar(label="Magnitud límite")
plt.xlabel("RA(º)")
plt.ylabel("DEC(º)")
plt.title("Depth Map")
plt.show()

