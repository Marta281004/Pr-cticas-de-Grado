#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#magnitud frente a error

#quiero representar magnitud frente a sus error, espcifiando flags y todo desde terminal


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse # importamos esto para poder seleccionar cosas externas al código cuando lo ejecutamos en terminal

parser = argparse.ArgumentParser(description="DES vs Rubin vs Euclid")
parser.add_argument("--filename", default="tablafinal.fits", type=str, help="Este el fichero de datos") #Selección de archivo
parser.add_argument("--flags", nargs="+", default=None, help="nombre de flags para aplicar. EJ: flags_footprint=1") #selección de flags
parser.add_argument("--mag", required=True, help="Nombre columna magnitud")
parser.add_argument("--err", required=True, help="Nombre columna error")
parser.add_argument("--color", default="blue", type=str, help="Color del scatter")
parser.add_argument("--xlim", nargs=2, type=float, default=None, help="Límites eje x: min max")
parser.add_argument("--ylim", nargs=2, type=float, default=None, help="Límites eje y: min max")
args = parser.parse_args()    


#Ya se ha seleccionado el archivo desde fuera de la terminal, ajora se lee

hdul = fits.open(args.filename)
datos = hdul[1].data


# extraer magnitudes

mag = np.array(datos[args.mag])
err = np.array(datos[args.err])


#mascara de flags

def crear_mascara_flags(datos, flag_dict=None):
    if flag_dict is None or len(flag_dict) == 0:
        return np.ones(len(datos), dtype=bool)  # No filtra nada

    mask = np.ones(len(datos), dtype=bool)

    for flag_name, flag_value in flag_dict.items():
        if flag_name not in datos.columns.names:
            print(f"Warning: {flag_name} no existe en el FITS")
            continue

        flag_array = np.array(datos[flag_name])
        mask &= (flag_array == flag_value)

    return mask


# definición de la función main

def main():
    # Leer FITS
    hdul = fits.open(args.filename)
    datos = hdul[1].data

    # Extraer columnas
    mag = np.array(datos[args.mag])
    err = np.array(datos[args.err])

    # Construir diccionario de flags dinámicamente
    flag_dict = {}
    if args.flags is not None:
        for f in args.flags:
            if "=" in f:
                key, value = f.split("=")
                if value.lower() == "true":
                    flag_dict[key] = True
                elif value.lower() == "false":
                    flag_dict[key] = False
                else:
                    flag_dict[key] = int(value)

    # Crear máscara
    mask_flags = crear_mascara_flags(datos, flag_dict)


    # Aplicar máscara
    mag_clean = mag[mask_flags]
    err_clean = err[mask_flags]
    
    # Representación diagrama densidad
    plt.figure()
    hb = plt.hexbin(mag_clean, err_clean, gridsize=60, cmap="viridis", mincnt=1, bins='log')
    plt.colorbar(hb, label="Número de objetos")
    if args.xlim is not None:
        plt.xlim(args.xlim[0], args.xlim[1])
    if args.ylim is not None:
        plt.ylim(args.ylim[0], args.ylim[1])
    plt.xlabel("Magnitud")
    plt.ylabel("Error")
    plt.title("Mag vs Magerr (densidad)")

    # Representación scatter
    plt.figure()
    plt.scatter(mag_clean, err_clean, s=5, alpha=0.5, color=args.color)
    if args.xlim is not None:
        plt.xlim(args.xlim[0], args.xlim[1])
    if args.ylim is not None:
        plt.ylim(args.ylim[0], args.ylim[1])
    plt.xlabel("Magnitud")
    plt.ylabel("Error")
    plt.title("Mag vs magerr")
    
    # Magnitud límite
    
    mask_valid = np.isfinite(mag_clean) & np.isfinite(err_clean)
    mag_clean = mag_clean[mask_valid]
    err_clean = err_clean[mask_valid]
    
    mask_err_range = (err_clean >= 0.1) & (err_clean <= 0.12)
    mag_err_lim = mag_clean[mask_err_range] # nos quedamos solo con las magnitudes en el rango de errores propuesto
    
    plt.figure()
    bins = np.arange(np.nanmin(mag_err_lim), np.nanmax(mag_err_lim) + 0.1, 0.1)
    counts, edges, _ = plt.hist(mag_err_lim, bins=bins, histtype="stepfilled", alpha=0.7, color=args.color)
    plt.xlabel("Magnitud")
    plt.ylabel("Número de objetos")
    plt.title("Histograma para 0.1<=err<=0.12")
    
    plt.show()
    
    bin_centers = (edges[:-1] + edges[1:]) / 2
    max_index = np.argmax(counts)
    mag_max = bin_centers[max_index]
    
    n_boot = 200 #terminar de mirar esto
    mag_lim_samples = []

    for i in range(n_boot):
        sample = np.random.choice(mag_err_lim, size=len(mag_err_lim), replace=True)
        counts_boot, edges_boot = np.histogram(sample, bins=bins)
        bin_centers_boot = (edges_boot[:-1] + edges_boot[1:]) / 2
        max_index_boot = np.argmax(counts_boot)
        mag_lim_boot = bin_centers_boot[max_index_boot]
        mag_lim_samples.append(mag_lim_boot)
        
    mag_lim_samples = np.array(mag_lim_samples)
    mag_lim_mean = np.mean(mag_lim_samples)
    mag_lim_std = np.std(mag_lim_samples)

    print("Magnitud límite:", mag_lim_mean)
    print("Incertidumbre (bootstrap):", mag_lim_std)

    
if __name__ == "__main__":
    main()    

