import matplotlib.pyplot as plt
import numpy as np
import urllib
SDSS_File = '/Users/compastro/jenkins/SDSS_4363+3_z+0.04.csv'
SDSS_Data = np.genfromtxt(SDSS_File,skip_header=2, delimiter = ',',dtype=float,unpack=True)
NII_6583 = SDSS_Data[28,:]
Ha_6562 = SDSS_Data[27,:]
OIII_5006 = SDSS_Data[20,:]
Hb_4861 = SDSS_Data[18,:]
OIII_4363 = SDSS_Data[14,:]
OIII_Hb = np.log10(OIII_5006/Hb_4861)
NII_Ha = np.log10(NII_6583/Ha_6562)
plt.figure()
plt.xlim(-1.5,0.5)
plt.ylim(-1,1.5)
#plt.scatter(NII_Ha,OIII_Hb,s=30,c='b')
x=np.linspace(-1.5,0.3,50)
y=((.61/(x-.47))+1.19)
plt.plot(x,y,color='k')
x3=np.linspace(-1,-0.2,50)
y3=((.61/(x3-.05)+1.3))
plt.plot(x3,y3,linestyle='--',color='red')
counter=0
for i in range(0,len(SDSS_Data[0,:])):
    if OIII_5006[i]/OIII_4363[i]<100.0:
        cool = plt.scatter(NII_Ha[i],OIII_Hb[i],color='r')
        counter=counter+1
        #print ("madeit")
    elif OIII_5006[i]/OIII_4363[i]>100.0 and OIII_5006[i]/OIII_4363[i]<1000.0:
        mid_temp = plt.scatter(NII_Ha[i],OIII_Hb[i],color='g')
        counter=counter+1
        #print("k")
    elif OIII_5006[i]/OIII_4363[i]>1000.0:
        hot = plt.scatter(NII_Ha[i],OIII_Hb[i],color='k')
        counter=counter+1
        #print ("r")
    else:
        print ("error")
print(counter)
plt.legend((cool, mid_temp, hot), ('Low Te','Mid Te','High Te'), scatterpoints = 1, loc = 'lower left', ncol = 3, fontsize =8)
plt.ylabel(r"log ([OIII] $\lambda$5007/H$\beta$)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("BPT Diagram (AoN OIII_4363 > 3.0)")
plt.show()
#subplots showing two plots of the same file next to each other, one BPT and one 5007/4363 (temp) vs. NII/Hb (ionization)
#in notebook, write about trends shown in plots (low temp AGN, high temp mid ionization SF, etc.)
