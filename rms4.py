#!/usr/bin/env python
# coding: utf-8

# In[7]:


# rms.py 

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse # importamos esto para poder seleccionar cosas externas al código cuando lo ejecutamos en terminal

parser = argparse.ArgumentParser(description="DES vs Rubin vs Euclid")
parser.add_argument("--filename", default="tablafinal.fits", type=str, help="Este el fichero de datos") #Selección de archivo
parser.add_argument("--flags", nargs="+", default=None, help="nombre de flags para aplicar. EJ: flags_footprint=1") #selección de flags
# parser.add_argument("--flag_value", default="0", type=int, help="Valor de los flags seleccionados")
args = parser.parse_args()    

#Ya se ha seleccionado el archivo desde fuera de la terminal, ajora se lee

hdul = fits.open(args.filename)
# hdul = fits.open("tablafinal.fits")
datos = hdul[1].data


# Como output vamos a sacar dos gráficas mag_DES/Rubin-mag_Euclid en función de mag_euclid

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
            
    
            
# gráfica de bins

def grafica_bin(mag_ref, delta, bin_size=0.5, label_survey="Survey", label_ref="Referencia"):
    bin_center, mean, std, rms, n = calc_bin_stats(mag_ref, delta, bin_size)
    
    plt.errorbar(bin_center, mean, yerr=std, fmt="o", capsize=3, label=f"{label_survey} vs {label_ref}")
    # plt.xlabel(f"Mag {label_ref}")
    # plt.ylabel(f"Mag {label_survey} - Mag {label_ref}")
    # plt.title(f"Delta_mag bineada: {label_survey} vs {label_ref}")
    
    
    
#Tabla de bins  
  
def imprimir_tabla_bins(mag_ref, delta, titulo, bin_size=0.5):

    bin_c, mean, std, rms, n = calc_bin_stats(mag_ref, delta, bin_size)

    print(f"\n===== {titulo} =====")
    print("Bin_center   N_obj   Mean     Std      RMS")

    for i in range(len(bin_c)):
        if n[i] > 0:
            print(f"{bin_c[i]:8.2f}   {n[i]:6d}   "
                  f"{mean[i]:7.4f}   {std[i]:7.4f}   {rms[i]:7.4f}")


    
# eliminado de flags desde notebook

#def quitar_flags(datos):
    
 #   #Flags DES
 #   flag_footprint = np.array(datos["flags_footprint"])
 #   flag_foreground = np.array(datos["flags_foreground"])
 #   flag_gold = np.array(datos["flags_gold"])
    
    #Flags Euclid
 #   flag_extended = np.array(datos["extended_flag"])
 #   flag_spurious = np.array(datos["spurious_flag"])
 #   flag_det_quality = np.array(datos["det_quality_flag"])
    
    #Flags Rubin
    
    
    #máscara de flags
 #   mask = (flag_footprint == 1) & (flag_foreground == 0) & (flag_gold == 0) & (flag_spurious == 0) & (flag_det_quality == 0)
    
 #   return mask
    

    
    
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
    mag_DES, mag_Euclid, mag_Rubin = extraer_mag(datos, DES="bdf_mag_r_corrected", Euclid="m_vis_sersic", Rubin="r_cmodel_mag")
    #eliminado de flags desde notebook: volver a escribirlo
    
    
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
    delta_DES = calc_delta(mag_DES, mag_Euclid)
    delta_Rubin = calc_delta(mag_Rubin, mag_Euclid)
    
    imprimir_tabla_bins(mag_Euclid, delta_DES, "DES vs Euclid (SIN FLAGS)")
    imprimir_tabla_bins(mag_Euclid, delta_Rubin, "DES vs Euclid (SIN FLAGS)")
    
    grafica_deltam(mag_Euclid, delta_DES, label_survey="DES", label_ref="Euclid", xlim=(15,27), ylim=(-10,10)) 
    grafica_deltam(mag_Euclid, delta_Rubin, label_survey="Rubin", label_ref="Euclid", xlim=(15,27), ylim=(-10,10)) 
    
    plt.figure()
    
    grafica_bin(mag_Euclid, delta_DES, bin_size=0.5, label_survey="DES", label_ref="Euclid")
    grafica_bin(mag_Euclid, delta_Rubin, bin_size=0.5, label_survey="Rubin", label_ref="Euclid")
    
    plt.xlabel("mag Euclid")
    plt.ylabel("mag survey - mag Euclid")
    plt.title("Delta_mag bineada: DES/Rubin vs Euclid")
    plt.legend()
    
    #Cálculo datos clean
    delta_DES_clean = calc_delta(mag_DES_clean, mag_Euclid_clean)
    delta_Rubin_clean = calc_delta(mag_Rubin_clean, mag_Euclid_clean)
    
    imprimir_tabla_bins(mag_Euclid_clean, delta_DES_clean, "DES vs Euclid (CON FLAGS)")
    imprimir_tabla_bins(mag_Euclid_clean, delta_Rubin_clean, "DES vs Euclid (CON FLAGS)")
    
    grafica_deltam(mag_Euclid_clean, delta_DES_clean, label_survey="DES", label_ref="Euclid", xlim=(15,27), ylim=(-10,10)) 
    grafica_deltam(mag_Euclid_clean, delta_Rubin_clean, label_survey="Rubin", label_ref="Euclid", xlim=(15,27), ylim=(-10,10)) 
    
    plt.figure()
    
    grafica_bin(mag_Euclid_clean, delta_DES_clean, bin_size=0.5, label_survey="DES", label_ref="Euclid")
    grafica_bin(mag_Euclid_clean, delta_Rubin_clean, bin_size=0.5, label_survey="Rubin", label_ref="Euclid")
    
    plt.xlabel("mag Euclid")
    plt.ylabel("mag survey - mag Euclid")
    plt.title("Delta_mag bineada: DES/Rubin vs Euclid")
    plt.legend()
    
    #comparación con y sin flags DES VS EUCLID
    
    plt.figure() #grafica de scatter comparando con y sin
    plt.scatter(mag_Euclid, delta_DES, s=5, alpha=0.3, label="Sin flags")
    plt.scatter(mag_Euclid_clean, delta_DES_clean, s=5, alpha=0.6, label="Con flags")
    plt.xlim(15,27)
    plt.ylim(-10,10)
    plt.xlabel("Mag Euclid")
    plt.ylabel("mag DES - Mag Euclid")
    plt.legend()
    
    plt.figure() # gráfica binneado con y sin flags
    bin_center_all, mean_list_all, std_list_all, rms_list_all, n_list_all = calc_bin_stats(mag_Euclid, delta_DES, bin_size=0.5)
    plt.errorbar(bin_center_all, mean_list_all, yerr=std_list_all, fmt="o", capsize=3, label="Sin flags")
    bin_center_clean, mean_list_clean, std_list_clean, rms_list_clean, n_list_clean = calc_bin_stats(mag_Euclid_clean, delta_DES_clean, bin_size=0.5)
    plt.errorbar(bin_center_clean, mean_list_clean, yerr=std_list_clean, fmt="o", capsize=3, label="Con flags")
    plt.xlabel("Mag Euclid")
    plt.ylabel("Mag DES - Mag Euclid")
    plt.title("DES vs Euclid: Δmag bineada con y sin flags")
    plt.legend()
    
    #comparación con y sin flags RUBIN VS EUCLID
    
    plt.figure() #grafica de scatter comparando con y sin
    plt.scatter(mag_Euclid, delta_Rubin, s=5, alpha=0.3, label="Sin flags")
    plt.scatter(mag_Euclid_clean, delta_Rubin_clean, s=5, alpha=0.6, label="Con flags")
    plt.xlim(15,27)
    plt.ylim(-10,10)
    plt.xlabel("Mag Euclid")
    plt.ylabel("Mag Rubin - Mag Euclid")
    plt.legend()
    
    plt.figure() # gráfica binneado con y sin flags
    bin_center_all, mean_list_all, std_list_all, rms_list_all, n_list_all = calc_bin_stats(mag_Euclid, delta_Rubin, bin_size=0.5)
    plt.errorbar(bin_center_all, mean_list_all, yerr=std_list_all, fmt="o", capsize=3, label="Sin flags")
    bin_center_clean, mean_list_clean, std_list_clean, rms_list_clean, n_list_clean = calc_bin_stats(mag_Euclid_clean, delta_Rubin_clean, bin_size=0.5)
    plt.errorbar(bin_center_clean, mean_list_clean, yerr=std_list_clean, fmt="o", capsize=3, label="Con flags")
    plt.xlabel("Mag Euclid")
    plt.ylabel("Mag Rubin - Mag Euclid")
    plt.title("Rubin vs Euclid: Δmag bineada con y sin flags")
    plt.legend()
    
    plt.show()
    
if __name__ == "__main__":
    main()


# In[ ]:




