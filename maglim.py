#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Histograma: número de objetos frente a la magnitud medida en cada filtro

#con esto quiero sacar la magnitud límite


from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse # importamos esto para poder seleccionar cosas externas al código cuando lo ejecutamos en terminal

parser = argparse.ArgumentParser(description="DES vs Rubin vs Euclid")
parser.add_argument("--filename", default="tablafinal.fits", type=str, help="Este el fichero de datos") #Selección de archivo
parser.add_argument("--flags", nargs="+", default=None, help="nombre de flags para aplicar. EJ: flags_footprint=1") #selección de flags
parser.add_argument("--mag", required=True, help="Nombre columna magnitud")
parser.add_argument("--color", default="blue", type=str, help="Color del scatter")
parser.add_argument("--mag_range", nargs=2, type=float, default=None, help="Rango de mag: min max")
parser.add_argument("--bins", default=0.2, type=float, help="tamaño del bin en mag")
parser.add_argument("--yscale", default="linear", choices=["linear", "log"], help="Escala del eje Y")

args = parser.parse_args()   


#Ya se ha seleccionado el archivo desde fuera de la terminal, ajora se lee

hdul = fits.open(args.filename)
datos = hdul[1].data


# extraer magnitudes

mag = np.array(datos[args.mag])


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
    mask_valid = np.isfinite(mag_clean)
    mag_clean = mag_clean[mask_valid]
    if args.mag_range is not None:
        mag_min, mag_max = args.mag_range
        mag_clean = mag_clean[(mag_clean > mag_min) & (mag_clean < mag_max)]
    
    # Histograma
    min_mag = np.nanmin(mag_clean)
    max_mag = np.nanmax(mag_clean)

    bins = np.arange(min_mag, max_mag + args.bins, args.bins)

    counts, edges = np.histogram(mag_clean, bins=bins)

    bin_centers = (edges[:-1] + edges[1:]) / 2

    
    # Representación 1
    plt.figure()
    plt.plot(bin_centers, counts, color=args.color)
    plt.yscale(args.yscale)
    plt.xlabel("Magnitud")
    plt.ylabel("Número de objetos")
    plt.title(f"Curva suavizada de {args.mag}")

    
    # Representación 2
    # Histograma en modo barras
    plt.figure()
    plt.hist(mag_clean, bins=bins, color=args.color, histtype="stepfilled", alpha=0.7)
    plt.yscale(args.yscale)
    plt.xlabel("Magnitud")
    plt.ylabel("Número de objetos")
    plt.title(f"Histograma (barras) de {args.mag}")

    plt.show()
    
    #Magnitud limite
    max_index = np.argmax(counts)
    mag_limite = bin_centers[max_index]

    print("Magnitud donde el histograma alcanza el máximo:", mag_limite)
    
if __name__ == "__main__":
    main()    

