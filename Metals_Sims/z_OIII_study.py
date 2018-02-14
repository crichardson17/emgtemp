import matplotlib.pyplot as plt
import numpy as np
import urllib
SDSS_File = '/Users/Sam/Documents/emgtemp/data/4363_gr_5_0_err.csv'
Dered_File = '/Users/Sam/Documents/emgtemp/data/4363_gr_5_0_err_dered.csv'
Dered_Data = np.genfromtxt(Dered_File,skip_header=1, delimiter = ',',dtype=float,unpack=True,names = True)
SDSS_Data = np.genfromtxt(SDSS_File,skip_header=1, delimiter = ',',dtype=float,unpack=True,names = True)
OIII_5006 = Dered_Data['Flux_OIII_5006']
Ha_6562 = SDSS_Data['Flux_Ha_6562']
Hb_4861 = SDSS_Data['Flux_Hb_4861']
OIII_4363 = Dered_Data['Flux_OIII_4363']
z = SDSS_Data['z']
Temp_Ratio = np.log10(OIII_5006/OIII_4363)
Ha_Hb_Ratio = np.log10(Ha_6562/Hb_4861)
zlog = np.log10(z)
#####################################################################################################
def getColor(OIII_5006, OIII_4363):
        Temp_Color = 'k'
        if OIII_5006/OIII_4363<75:
            #Temp_Color = '0.25'
            Temp_Color = 'red'
       
        elif OIII_5006/OIII_4363>75:
            #Temp_Color = '0.75'
            Temp_Color = 'blue'

        #black = black + 1
        else:
            print ("error")
        return Temp_Color
#####################################################################################################

fig = plt.subplot(221)
#fig.subplots_adjust(wspace=0.4,hspace=0.4)
sp1 = plt.subplot(221)
plt.scatter(z, Temp_Ratio)
plt.xlabel('z')
plt.ylabel('OIII Ratio')

sp2 = plt.subplot(222)
plt.scatter(z, Ha_Hb_Ratio)
plt.xlabel('z')
plt.ylabel(r"H$\alpha$/H$\beta$")

sp3 = plt.subplot(223)
for i in range(0,len(SDSS_Data['z'])):
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    #print(Temp_Color)
    plt.scatter(z[i],Ha_Hb_Ratio[i], color = Temp_Color, edgecolor = 'none')
    #print (Temp_Color)    
plt.xlabel('z')
plt.ylabel(r"H$\alpha$/H$\beta$")
    
sp3 = plt.subplot(224)
for i in range(0,len(SDSS_Data['z'])):
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    #print(Temp_Color)
    plt.scatter(zlog[i],Ha_Hb_Ratio[i], color = Temp_Color, edgecolor = 'none')
    #print (Temp_Color)  
plt.xlabel('log z')
plt.ylabel(r"H$\alpha$/H$\beta$")
plt.show()
plt.savefig("z_OIII_Ha_Study_plot.pdf")
