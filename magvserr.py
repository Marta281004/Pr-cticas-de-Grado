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

    # Representación
    plt.figure()
    plt.scatter(mag_clean, err_clean, s=5, alpha=0.5, color=args.color)
    if args.xlim is not None:
        plt.xlim(args.xlim[0], args.xlim[1])
    if args.ylim is not None:
        plt.ylim(args.ylim[0], args.ylim[1])
    plt.show()
    
if __name__ == "__main__":
    main()    

