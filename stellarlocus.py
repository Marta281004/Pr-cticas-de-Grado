#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#stellarlocus

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse # importamos esto para poder seleccionar cosas externas al código cuando lo ejecutamos en terminal

parser = argparse.ArgumentParser(description="DES vs Rubin vs Euclid")
parser.add_argument("--filename", default="tablafinal.fits", type=str, help="Este el fichero de datos") #Selección de archivo
parser.add_argument("--flags", nargs="+", default=None, help="nombre de flags para aplicar. EJ: flags_footprint=1") #selección de flags
args = parser.parse_args()   



#Ya se ha seleccionado el archivo desde fuera de la terminal, ajora se lee

hdul = fits.open(args.filename)
# hdul = fits.open("tablafinal.fits")
datos = hdul[1].data


#para extraer las magnitudes

def extraer_mag(datos, DES_g="psf_mag_aper_8_g_corrected", DES_r="psf_mag_aper_8_r_corrected", DES_i="psf_mag_aper_8_i_corrected",
                Rubin_g="g_psfmag", Rubin_r="r_psfmag", Rubin_i="i_psfmag"):
    mag_DES_g = np.array(datos[DES_g])
    mag_DES_r = np.array(datos[DES_r])
    mag_DES_i = np.array(datos[DES_i])
    mag_Rubin_g = np.array(datos[Rubin_g])
    mag_Rubin_r = np.array(datos[Rubin_r])
    mag_Rubin_i = np.array(datos[Rubin_i])
    return mag_DES_g, mag_DES_r, mag_DES_i, mag_Rubin_g, mag_Rubin_r, mag_Rubin_i


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

#cálculo de índices de color

def calc_color(mag_g,mag_r,mag_i):
    g_r = mag_g - mag_r
    r_i = mag_r - mag_i
    return g_r, r_i


#representar g_r vs r_i

def plot_stellar_locus(g_r, r_i, label="Survey", gridsize=60, cmap="viridis", xlim=None, ylim=None):

    # Limpiar NaNs / inf
    mask = np.isfinite(g_r) & np.isfinite(r_i)
    g_r = g_r[mask]
    r_i = r_i[mask]

    # Plot
    plt.figure()
    hb = plt.hexbin(g_r, r_i, gridsize=gridsize, cmap=cmap, mincnt=1, bins='log')  
    plt.colorbar(hb, label="Número de objetos")
    plt.xlabel("g - r")
    plt.ylabel("r - i")
    plt.title(f"Stellar locus ({label})")
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
        

def main():
    
    mag_DES_g, mag_DES_r, mag_DES_i, mag_Rubin_g, mag_Rubin_r, mag_Rubin_i = extraer_mag(datos, DES_g="psf_mag_aper_8_g_corrected", DES_r="psf_mag_aper_8_r_corrected", DES_i="psf_mag_aper_8_i_corrected",
                Rubin_g="g_psfmag", Rubin_r="r_psfmag", Rubin_i="i_psfmag")
    
    #para convertir las entradas de la terminal en un "diccionario"
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
                
    #crear máscara
    mask_flags = crear_mascara_flags(datos, flag_dict)

    mag_DES_g_clean = mag_DES_g[mask_flags]
    mag_DES_r_clean = mag_DES_r[mask_flags]
    mag_DES_i_clean = mag_DES_i[mask_flags]
    mag_Rubin_g_clean = mag_Rubin_g[mask_flags]
    mag_Rubin_r_clean = mag_Rubin_r[mask_flags]
    mag_Rubin_i_clean = mag_Rubin_i[mask_flags]
    
    #cálculo de los índices de color
    g_r_DES, r_i_DES = calc_color(mag_DES_g_clean, mag_DES_r_clean, mag_DES_i_clean)
    g_r_Rubin, r_i_Rubin = calc_color(mag_Rubin_g_clean, mag_Rubin_r_clean, mag_Rubin_i_clean)
    #salen valores muy raros de color, así que me quedo con estos intervalos
    mask_color_DES = (g_r_DES > -1) & (g_r_DES < 3) & (r_i_DES > -1) & (r_i_DES < 2)
    g_r_DES = g_r_DES[mask_color_DES]
    r_i_DES = r_i_DES[mask_color_DES]
    mask_color_Rubin = (g_r_Rubin > -1) & (g_r_Rubin < 3) & (r_i_Rubin > -1) & (r_i_Rubin < 2)
    g_r_Rubin = g_r_Rubin[mask_color_Rubin]
    r_i_Rubin = r_i_Rubin[mask_color_Rubin]
    
    #representación
    plot_stellar_locus(g_r_DES, r_i_DES, label="DES", gridsize=60, cmap="viridis", xlim=(-0.5, 2), ylim=(-0.5, 1.5))
    plot_stellar_locus(g_r_Rubin, r_i_Rubin, label="Rubin", gridsize=60, cmap="viridis", xlim=(-0.5, 1.5), ylim=(-0.5, 1.5))
    plt.show()
if __name__ == "__main__":
    main()

