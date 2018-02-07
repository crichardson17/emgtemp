import matplotlib.pyplot as plt
import numpy as np
import urllib
SDSS_File = '/Users/Sam/Documents/emgtemp/data/4363_gr_5_0_err.csv'
SDSS_Data = np.genfromtxt(SDSS_File,skip_header=1, delimiter = ',',dtype=float,unpack=True,names = True)
OIII_5006 = SDSS_Data['Flux_OIII_5006']
Ha_6562 = SDSS_Data['Flux_Ha_6562']
Hb_4861 = SDSS_Data['Flux_Hb_4861']
OIII_4363 = SDSS_Data['Flux_OIII_4363']
OIII_Hb = np.log10(OIII_5006/Hb_4861)
Temp_Ratio = np.log10(OIII_5006/OIII_4363)
Ha_Hb_Ratio = np.log10(Ha_6562/Hb_4861)
z = SDSS_Data['z']
OIII_Hb = np.log10(OIII_5006/Hb_4861)

fig = plt.subplot(211)
#fig.subplots_adjust(wspace=0.4,hspace=0.4)
sp1 = plt.subplot(311)
plt.scatter(z, Temp_Ratio)
plt.xlabel('z')
plt.ylabel('OIII Ratio')

sp2 = plt.subplot(212)
plt.scatter(z, Ha_Hb_Ratio)
plt.xlabel('z')
plt.ylabel(r"H$\alpha$/H$\beta$")

#sp1 = plt.subplot(413)
#plt.scatter(z, Temp_Ratio)
#plt.xlabel('z')
#plt.ylabel('OIII Ratio')
plt.show()
plt.savefig("z_OIII_Ha_Study_plot.pdf")
