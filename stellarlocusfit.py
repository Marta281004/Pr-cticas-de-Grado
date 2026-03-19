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
parser.add_argument("--g_r_min", type=float, default=-1, help="valor mínimo de g-r para el ajuste")
parser.add_argument("--g_r_max", type=float, default=3, help="valor máximo de g-r para el ajuste")
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

def plot_stellar_locus(g_r, r_i, label="Survey", gridsize=60, cmap="viridis", xlim=None, ylim=None, m=None, b=None):
    mask = np.isfinite(g_r) & np.isfinite(r_i)
    g_r = g_r[mask]
    r_i = r_i[mask]

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

    if m is not None and b is not None:
        x = np.linspace(min(g_r), max(g_r), 100)
        y = m*x + b
        plt.plot(x, y, color="red", lw=2, label="Fit")
        plt.legend()
        

#ajuste
def fit_stellar_locus(g_r, r_i):
    # Ajuste lineal: r_i = m * g_r + b
    p, cov = np.polyfit(g_r, r_i, 1, cov=True)
    m = p[0]
    b = p[1]
    m_err = np.sqrt(cov[0, 0])
    b_err = np.sqrt(cov[1, 1])

    # Calculamos r_i ajustado
    r_i_fit = m * g_r + b

    # Error perpendicular a la línea: sqrt((r_i - r_i_fit)^2 / (1 + m^2))
    sigma_perp = np.sqrt(np.mean((r_i - r_i_fit)**2 / (1 + m**2)))
    
    #R²
    ss_res = np.sum((r_i - r_i_fit)**2)
    ss_tot = np.sum((r_i - np.mean(r_i))**2)
    r2 = 1 - ss_res / ss_tot

    return m, b, m_err, b_err, sigma_perp, r2


        
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
    
    #ajuste y dispersión
    #filtrar por g-r para el ajuste
    fit_mask_DES = (g_r_DES >= args.g_r_min) & (g_r_DES <= args.g_r_max)
    g_r_DES_fit = g_r_DES[fit_mask_DES]
    r_i_DES_fit = r_i_DES[fit_mask_DES]
    fit_mask_Rubin = (g_r_Rubin >= args.g_r_min) & (g_r_Rubin <= args.g_r_max)
    g_r_Rubin_fit = g_r_Rubin[fit_mask_Rubin]
    r_i_Rubin_fit = r_i_Rubin[fit_mask_Rubin]
    #ajuste
    m_DES, b_DES, m_err_DES, b_err_DES, sigma_DES, r2_DES = fit_stellar_locus(g_r_DES_fit, r_i_DES_fit)
    m_Rubin, b_Rubin, m_err_Rubin, b_err_Rubin, sigma_Rubin, r2_Rubin = fit_stellar_locus(g_r_Rubin_fit, r_i_Rubin_fit)
    
    print("\n===== RESULTADOS =====")

    print(f"DES:")
    print(f"  pendiente = {m_DES:.3f} ± {m_err_DES:.3f}")
    print(f"  intercepto = {b_DES:.3f} ± {b_err_DES:.3f}")
    print(f"  dispersión = {sigma_DES:.3f}")
    print(f"  R^2 = {r2_DES:.4f}")

    print(f"\nRubin:")
    print(f"  pendiente = {m_Rubin:.3f} ± {m_err_Rubin:.3f}")
    print(f"  intercepto = {b_Rubin:.3f} ± {b_err_Rubin:.3f}")
    print(f"  dispersión = {sigma_Rubin:.3f}")
    print(f"  R^2 = {r2_Rubin:.4f}")
    
    #representación
    plot_stellar_locus(g_r_DES, r_i_DES, label="DES", gridsize=60, cmap="viridis", xlim=(-0.5, 2), ylim=(-0.5, 1.5), m=m_DES, b=b_DES)
    plot_stellar_locus(g_r_Rubin, r_i_Rubin, label="Rubin", gridsize=60, cmap="viridis", xlim=(-0.5, 2), ylim=(-0.5, 1.5), m=m_Rubin, b=b_Rubin)
    plt.show()
if __name__ == "__main__":
    main()

