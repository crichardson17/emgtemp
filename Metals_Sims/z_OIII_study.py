import matplotlib.pyplot as plt
import numpy as np
import urllib
SDSS_File = '/Users/Sam/Documents/emgtemp/data/4363_gr_5_0_err.csv'
SDSS_Data = np.genfromtxt(SDSS_File,skip_header=2, delimiter = ',',dtype=float,unpack=True,names = True)
OI_6300 = SDSS_Data['Flux_OI_6300']
OIII_5006 = SDSS_Data['Flux_OIII_5006']
Ha_6562 = SDSS_Data['Flux_Ha_6562']
Hb_4861 = SDSS_Data['Flux_Hb_4861']
OIII_4363 = SDSS_Data['Flux_OIII_4363']
OIII_Hb = np.log10(OIII_5006/Hb_4861)
Temp_Ratio = np.log10(OIII_5006/OIII_4363)
z = SDSS_Data['z']
OIII_Hb = np.log10(OIII_5006/Hb_4861)
plt.figure()
