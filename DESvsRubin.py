#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# rms.py 

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse # importamos esto para poder seleccionar cosas externas al código cuando lo ejecutamos en terminal

parser = argparse.ArgumentParser(description="DES vs Rubin")
parser.add_argument("--filename", default="tablafinal.fits", type=str, help="Este el fichero de datos") #Selección de archivo
parser.add_argument("--flags", nargs="+", default=None, help="nombre de flags para aplicar. EJ: flags_footprint=1") #selección de flags
# parser.add_argument("--flag_value", default="0", type=int, help="Valor de los flags seleccionados")
parser.add_argument("--save_tables", default=None, type=str, help="Nombre del archivo txt donde guardar las tablas")
args = parser.parse_args()    

#Ya se ha seleccionado el archivo desde fuera de la terminal, ajora se lee

hdul = fits.open(args.filename)
# hdul = fits.open("tablafinal.fits")
datos = hdul[1].data


# Para "extraer las magnitudes": Devuelve las magnitudes como arrays numpy

def extraer_mag(datos, DES="bdf_mag_r_corrected", Euclid="m_vis_sersic", Rubin="r_cmodel_mag"):
    mag_DES = np.array(datos[DES])
    mag_Euclid = np.array(datos[Euclid])
    mag_Rubin = np.array(datos[Rubin])
    return mag_DES, mag_Euclid, mag_Rubin


# ahora queremos que calcule la diferencia de magnitudes entre un survey y otro de referencia (en principio, tomamos q1)

def calc_delta(mag_survey,mag_ref):
    delta = mag_survey - mag_ref
    return delta


# Representación gráfica

def grafica_deltam(mag_ref, delta, label_survey="Survey", label_ref="Survey ref", xlim=(15,27), ylim=(-10,10)): 
    plt.figure()
    hb = plt.hexbin(mag_ref, delta, gridsize=80, extent=(xlim[0], xlim[1], ylim[0], ylim[1]), cmap="viridis", mincnt=1)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel(f"Mag {label_ref}")
    plt.ylabel(f"Mag {label_survey} - Mag {label_ref}")
    plt.title(f"Diferencia de magnitudes: {label_survey} vs {label_ref}")
    plt.colorbar(hb, label="Número de objetos")
    
    
# crear bines eje x

def calc_bin_stats(mag_ref, delta, bin_size=1):
    min_mag = 15
    max_mag = 27
    bins = np.arange(min_mag, max_mag + bin_size, bin_size)
    bin_center = (bins[:-1] + bins[1:]) / 2

    mean_list = []
    std_list = []
    rms_list = []
    n_list = []

    for i in range(len(bins)-1):
        mask = ((mag_ref >= bins[i]) & (mag_ref < bins[i+1]) & (delta >= -10) & (delta < 10))
        delta_bin = delta[mask]

        if len(delta_bin) > 0:
            mean_list.append(np.mean(delta_bin))
            std_list.append(np.std(delta_bin))
            rms_list.append(np.sqrt(np.mean(delta_bin**2)))
            n_list.append(len(delta_bin))
        else:
            mean_list.append(np.nan)
            std_list.append(np.nan)
            rms_list.append(np.nan)
            n_list.append(0)

    return (np.array(bin_center), np.array(mean_list), np.array(std_list), np.array(rms_list), np.array(n_list)) 



# Grafica de contornos    
    
def grafica_contorno(mag_ref, delta, label_survey="Survey", label_ref="Survey ref", xlim=(15,27), ylim=(-0.5,0.5), bins=100):

    plt.figure()

    # Histograma 2D
    H, xedges, yedges = np.histogram2d(mag_ref, delta, bins=bins, range=[[xlim[0], xlim[1]], [ylim[0], ylim[1]]])

    # Centros de los bins
    X = (xedges[:-1] + xedges[1:]) / 2
    Y = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(X, Y)

    # Transponer H porque histogram2d organiza diferente
    H = H.T

    # Dibujar contornos
    plt.contourf(X, Y, H, levels=15, cmap="viridis")
    plt.colorbar(label="Número de objetos")

    plt.xlabel(f"Mag {label_ref}")
    plt.ylabel(f"Mag {label_survey} - Mag {label_ref}")
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title(f"Contornos de densidad: {label_survey} vs {label_ref}")

    
    
# gráfica de bins

def grafica_bin(mag_ref, delta, bin_size=0.5, label_survey="Survey", label_ref="Referencia"):
    bin_center, mean, std, rms, n = calc_bin_stats(mag_ref, delta, bin_size)
    
    plt.errorbar(bin_center, mean, yerr=std, fmt="o", capsize=3, label=f"{label_survey} vs {label_ref}")
    # plt.xlabel(f"Mag {label_ref}")
    # plt.ylabel(f"Mag {label_survey} - Mag {label_ref}")
    # plt.title(f"Delta_mag bineada: {label_survey} vs {label_ref}")
    
    
    
# Contornos+bins

def grafica_contorno_con_bins(mag_ref, delta, label_survey="Survey", label_ref="Survey ref", xlim=(15,27), ylim=(-0.5,0.5), bins=100, bin_size=0.5):

    plt.figure()

    # ---- CONTORNOS ----
    H, xedges, yedges = np.histogram2d(mag_ref, delta, bins=bins, range=[[xlim[0], xlim[1]], [ylim[0], ylim[1]]])

    X = (xedges[:-1] + xedges[1:]) / 2
    Y = (yedges[:-1] + yedges[1:]) / 2
    X, Y = np.meshgrid(X, Y)
    H = H.T

    plt.contourf(X, Y, H, levels=15, cmap="viridis")
    plt.colorbar(label="Número de objetos")

    # ---- BINS ENCIMA ----
    bin_center, delta_mean, delta_std, rms, n = calc_bin_stats(mag_ref, delta, bin_size)

    plt.errorbar(bin_center, delta_mean, yerr=delta_std, fmt="o", color="red", capsize=3, label="Media por bin")

    plt.axhline(0, color="white", linestyle="--", linewidth=1)

    plt.xlim(xlim)
    plt.ylim(ylim)

    plt.xlabel(f"Mag {label_ref}")
    plt.ylabel(f"Mag {label_survey} - Mag {label_ref}")
    plt.title(f"{label_survey} vs {label_ref}")

    plt.legend()
    
#Tabla de bins  
  
def imprimir_tabla_bins(mag_ref, delta, titulo, bin_size=0.5, file_handle=None):

    bin_c, mean, std, rms, n = calc_bin_stats(mag_ref, delta, bin_size)

    header = f"\n===== {titulo} =====\n"
    header += "Bin_center   N_obj   Mean     Std      RMS\n"

    print(header)

    if file_handle is not None:
        file_handle.write(header)

    for i in range(len(bin_c)):
        if n[i] > 0:
            line = (f"{bin_c[i]:8.2f}   {n[i]:6d}   "
                    f"{mean[i]:7.4f}   {std[i]:7.4f}   {rms[i]:7.4f}\n")

            print(line.strip())

            if file_handle is not None:
                file_handle.write(line)

   
    
#eliminado de flags desde terminal

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


# Función completa con todos los pasos anteriores

def main():
    
    file_handle = None

    if args.save_tables is not None:
        file_handle = open(args.save_tables, "w")
        print(f"Guardando tablas en {args.save_tables}")
    
    mag_DES, mag_Euclid, mag_Rubin = extraer_mag(datos, DES="bdf_mag_r_corrected", Euclid="m_vis_sersic", Rubin="r_cmodel_mag")
    
    
    #para convertir las entradas de la terminal en un "diccionario"
    flag_dict = {}

    if args.flags is not None:
        for f in args.flags:
            if "=" in f:
                try:
                    key, value = f.split("=")

                    # Detectar booleanos
                    if value.lower() == "true":
                        flag_dict[key] = True
                    elif value.lower() == "false":
                        flag_dict[key] = False
                    else:
                        flag_dict[key] = int(value)

                except ValueError:
                    print(f"Formato incorrecto en {f}, debería ser flag=value")

            else:
                # Si no se especifica valor:
                # por defecto 0 (o False si fuera booleano en el FITS)
                flag_dict[f] = 0
                
    #crear máscara
    mask_flags = crear_mascara_flags(datos, flag_dict)

    print("Objetos totales:", len(mag_DES))
    print("Objetos tras aplicar flags:", np.sum(mask_flags))

    mag_DES_clean = mag_DES[mask_flags]
    mag_Euclid_clean = mag_Euclid[mask_flags]
    mag_Rubin_clean = mag_Rubin[mask_flags]

    
    #cálculo datos sin filtrar 
    delta= calc_delta(mag_Rubin, mag_DES)
    
    imprimir_tabla_bins(mag_DES, delta, "Rubin vs DES (SIN FLAGS)", file_handle=file_handle)
    
    grafica_contorno(mag_DES, delta, label_survey="Rubin", label_ref="DES", xlim=(15,27), ylim=(-0.5,0.5))  
    
    grafica_contorno_con_bins(mag_DES, delta, label_survey="Rubin", label_ref="DES")
    
    plt.figure()
    
    grafica_bin(mag_DES, delta, bin_size=0.5, label_survey="Rubin", label_ref="DES")
    
    plt.xlabel("mag DES")
    plt.ylabel("mag survey - mag DES")
    plt.title("Delta_mag bineada: Rubin vs DES")
    plt.legend()
    
    #Cálculo datos clean
    delta_clean = calc_delta(mag_Rubin_clean, mag_DES_clean)
    
    imprimir_tabla_bins(mag_DES_clean, delta_clean, "Rubin vs DES (CON FLAGS)", file_handle=file_handle)
    
    grafica_contorno(mag_DES_clean, delta_clean, label_survey="Rubin", label_ref="DES", xlim=(15,27), ylim=(-0.5,0.5)) 
    
    grafica_contorno_con_bins(mag_DES_clean, delta_clean, label_survey="Rubin", label_ref="DES")
    
    plt.figure()
    
    grafica_bin(mag_DES_clean, delta_clean, bin_size=0.5, label_survey="Rubin", label_ref="DES")
    
    plt.xlabel("mag DES")
    plt.ylabel("mag survey - mag DES")
    plt.title("Delta_mag bineada: Rubin vs DES")
    plt.legend()
    
    #comparación con y sin flags DES VS EUCLID
    
    plt.figure() #grafica de scatter comparando con y sin
    plt.scatter(mag_DES, delta, s=5, alpha=0.3, label="Sin flags")
    plt.scatter(mag_DES_clean, delta_clean, s=5, alpha=0.6, label="Con flags")
    plt.xlim(15,27)
    plt.ylim(-10,10)
    plt.xlabel("Mag DES")
    plt.ylabel("mag Rubin - Mag DES")
    plt.legend()
    
    plt.figure() # gráfica binneado con y sin flags
    bin_center_all, mean_list_all, std_list_all, rms_list_all, n_list_all = calc_bin_stats(mag_DES, delta, bin_size=0.5)
    plt.errorbar(bin_center_all, mean_list_all, yerr=std_list_all, fmt="o", capsize=3, label="Sin flags")
    bin_center_clean, mean_list_clean, std_list_clean, rms_list_clean, n_list_clean = calc_bin_stats(mag_DES_clean, delta_clean, bin_size=0.5)
    plt.errorbar(bin_center_clean, mean_list_clean, yerr=std_list_clean, fmt="o", capsize=3, label="Con flags")
    plt.xlabel("Mag DES")
    plt.ylabel("Mag Rubin - Mag DES")
    plt.title("Rubin vs DES: Δmag bineada con y sin flags")
    plt.legend()
    
    
    if file_handle is not None:
        file_handle.close()
    
    plt.show()
    
if __name__ == "__main__":
    main()

