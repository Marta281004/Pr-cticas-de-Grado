#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#depth_map

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse 

#definimos los args de la termianl
parser = argparse.ArgumentParser(description="DES vs Rubin vs Euclid")
parser.add_argument("--filename", default="tablafinal.fits", type=str, help="Este el fichero de datos") #Selección de archivo
parser.add_argument("--flags", nargs="+", default=None, help="nombre de flags para aplicar. EJ: flags_footprint=1") #selección de flags
parser.add_argument("--mag", required=True, help="Nombre columna magnitud")
parser.add_argument("--ra", required=True, type=str, help="Nombre columna ascensión recta")
parser.add_argument("--dec", required=True, type=str, help="Nombre columna declinación")
parser.add_argument("--outfile", default="resultados.txt", type=str, help="Nombre del archivo txt donde guardar las tablas")
args = parser.parse_args()   


#eliminado de flags desde terminal

def crear_mascara_flags(datos, flag_dict=None):

    if flag_dict is None or len(flag_dict) == 0:
        return np.ones(len(datos), dtype=bool)

    mask = np.ones(len(datos), dtype=bool)

    for flag_name, (op, value) in flag_dict.items():

        if flag_name not in datos.columns.names:
            print(f"Warning: {flag_name} no existe en el FITS")
            continue

        flag_array = np.array(datos[flag_name])

        if op == "=":
            mask &= (flag_array == value)

        elif op == "<":
            mask &= (flag_array < value)

        elif op == ">":
            mask &= (flag_array > value)

        elif op == "<=":
            mask &= (flag_array <= value)

        elif op == ">=":
            mask &= (flag_array >= value)

    return mask

# función maglim reutilizable

def calcular_maglim(mag, bin_size=0.2, n_boot=100):
    mask = np.isfinite(mag)
    mag = mag[mask]
    
    if len(mag) == 0:
        return np.nan, np.nan
    
    bins = np.arange(np.nanmin(mag), np.nanmax(mag) + bin_size, bin_size)
    counts, edges = np.histogram(mag, bins=bins)
    bin_centers = (edges[:-1] + edges[1:]) / 2
    
    if len(counts) == 0:
        return np.nan, np.nan
    
    samples = []
    for i in range(n_boot):
        sample = np.random.choice(mag, size=len(mag), replace=True)
        counts_boot, _ = np.histogram(sample, bins=bins)
        max_idx_boot = np.argmax(counts_boot)
        samples.append(bin_centers[max_idx_boot])

    samples = np.array(samples)
    
    return np.mean(samples), np.std(samples)


# definición de la función main

def main():
    # Leer FITS
    hdul = fits.open(args.filename)
    datos = hdul[1].data

    # Extraer columnas
    mag = np.array(datos[args.mag])
    ra = np.array(datos[args.ra])
    dec = np.array(datos[args.dec])

    # Construir diccionario de flags dinámicamente
    flag_dict = {}

    if args.flags is not None:
        for f in args.flags:

            try:

                if "<=" in f:
                    key, value = f.split("<=")
                    flag_dict[key] = ("<=", float(value))

                elif ">=" in f:
                    key, value = f.split(">=")
                    flag_dict[key] = (">=", float(value))

                elif "<" in f:
                    key, value = f.split("<")
                    flag_dict[key] = ("<", float(value))

                elif ">" in f:
                    key, value = f.split(">")
                    flag_dict[key] = (">", float(value))

                elif "=" in f:
                    key, value = f.split("=")

                    if value.lower() == "true":
                        flag_dict[key] = ("=", True)

                    elif value.lower() == "false":
                        flag_dict[key] = ("=", False)

                    else:
                        flag_dict[key] = ("=", float(value))

                else:
                    flag_dict[f] = ("=", 0)

            except ValueError:
                print(f"Formato incorrecto en {f}")

                    
    # Crear y aplicar máscara
    mask_flags = crear_mascara_flags(datos, flag_dict)
    mag = mag[mask_flags]
    ra = ra[mask_flags]
    dec = dec[mask_flags]
    
    # definir tamaño de celda (en grados)
    cell_size = 0.07

    # crear bins espaciales
    ra_bins = np.arange(np.min(ra), np.max(ra) + cell_size, cell_size)
    dec_bins = np.arange(np.min(dec), np.max(dec) + cell_size, cell_size)
    
    # diccionario
    depths = {}
    cell_id = 0

    # mapa vacío
    # depth_map = np.full((len(ra_bins)-1, len(dec_bins)-1), np.nan)

    # loop sobre celdas
    for i in range(len(ra_bins)-1):
        for j in range(len(dec_bins)-1):

            mask_cell = ((ra >= ra_bins[i]) & (ra < ra_bins[i+1]) & (dec >= dec_bins[j]) & (dec < dec_bins[j+1]))

            mag_cell = mag[mask_cell]
            N_obj = len(mag_cell)

            # centro de la celda
            ra_center = 0.5 * (ra_bins[i] + ra_bins[i+1])
            dec_center = 0.5 * (dec_bins[j] + dec_bins[j+1])
            if N_obj == 0:
                maglim, err = np.nan, np.nan
            else:
                maglim, err = calcular_maglim(mag_cell) #calcular magnitud limite

            # guardar en diccionario
            depths[cell_id] = {"N": N_obj, "ra_center": ra_center, "dec_center": dec_center, "maglim": maglim, "err_maglim": err}

            cell_id += 1
          
    with open(args.outfile, "w") as f:
        f.write("cell_id N ra_center dec_center maglim err_maglim\n")
        for cell_id, data in depths.items():
            f.write(f"{cell_id} {data['N']} {data['ra_center']} {data['dec_center']} {data['maglim']} {data['err_maglim']}\n")
            
    print(f"Resultados guardados en {args.outfile}")        
                    

if __name__ == "__main__":
    main()        

