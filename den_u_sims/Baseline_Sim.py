import matplotlib.pyplot as plt
import numpy as np
import urllib
Low_Temp_Color = 'k'
Mid_Temp_Color = 'g'
High_Temp_Color = 'r'
Cloudy_Sim_Color = 'cyan'
markersize = 40
SDSS_File = '/Users/compastro/jenkins/SDSS_4363+5_z+0.04_dered_nospace.csv'
SDSS_Data = np.genfromtxt(SDSS_File,skip_header=1, delimiter = ',',dtype=float,unpack=True,names=True)
NII_6584 = SDSS_Data['Flux_NII_6583']
Ha_6562 = SDSS_Data['Flux_Ha_6562']
OI_6300 = SDSS_Data['Flux_OI_6300']
OIII_5006 = SDSS_Data['Flux_OIII_5006']
Hb_4861 = SDSS_Data['Flux_Hb_4861']
OIII_4363 = SDSS_Data['Flux_OIII_4363']
SII_6716 = SDSS_Data['Flux_SII_6716']
SII_6731 = SDSS_Data['Flux_SII_6730']
OII_3727 = SDSS_Data['Flux_OII_3726'] + SDSS_Data['Flux_OII_3728']
OIII_Hb = np.log10(OIII_5006/Hb_4861)
NII_Ha = np.log10(NII_6584/Ha_6562)
Temp_Ratio = np.log10(OIII_5006/OIII_4363)
S_Ratio = np.log10(SII_6716/SII_6731)
NO_Ratio = np.log10(NII_6584/OII_3727)
OI_Ratio = np.log10(OI_6300/Ha_6562)
O_Ratio = np.log10(OIII_5006/OII_3727)
S_Ha_Ratio = np.log10((SII_6716+SII_6731)/Ha_6562)
Cloudy_File = '/Users/compastro/cloudy/c13.03/runs/Baseline_Sim1.csv'
Cloudy_Data = np.genfromtxt(Cloudy_File, delimiter = ',',dtype=float,unpack=True,names=True)
Cloudy_NII_6584 = Cloudy_Data['N__2__6584A']
Cloudy_Ha_6562 = Cloudy_Data['H__1__6563A']
Cloudy_OIII_5006 = Cloudy_Data['O__3__5007A']
Cloudy_Hb_4861 = Cloudy_Data['TOTL__4861A']
Cloudy_OIII_4363 = Cloudy_Data['TOTL__4363A']
Cloudy_SII_6716 = Cloudy_Data['S_II__6716A']
Cloudy_SII_6731 = Cloudy_Data['S_II__6731A']
Cloudy_OII_3727 = Cloudy_Data['TOTL__3727A']
Cloudy_OI_6300 = Cloudy_Data['O__1__6300A']
Cloudy_OIII_Hb = np.log10(Cloudy_OIII_5006/Cloudy_Hb_4861)
Cloudy_NII_Ha = np.log10(Cloudy_NII_6584/Cloudy_Ha_6562)
Cloudy_Temp_Ratio = np.log10(Cloudy_OIII_5006/Cloudy_OIII_4363)
Cloudy_S_Ratio = np.log10(Cloudy_SII_6716/Cloudy_SII_6731)
Cloudy_NO_Ratio = np.log10(Cloudy_NII_6584/Cloudy_OII_3727)
fig = plt.figure(1)
fig.subplots_adjust(wspace=0.4,hspace=0.4)
sp1 = plt.subplot(221)
red = 0
green = 0
black = 0
for i in range(0,len(SDSS_Data['z'])):
    if OIII_5006[i]/OIII_4363[i]<50:
        plt.scatter(NII_Ha[i],OIII_Hb[i],color=High_Temp_Color, s = markersize)
        red = red + 1 
    elif OIII_5006[i]/OIII_4363[i]>50 and OIII_5006[i]/OIII_4363[i]<100:
        plt.scatter(NII_Ha[i],OIII_Hb[i],color=Mid_Temp_Color, s = markersize)
        green = green + 1
    elif OIII_5006[i]/OIII_4363[i]>100:
        plt.scatter(NII_Ha[i],OIII_Hb[i],color=Low_Temp_Color, s = markersize)
        black = black + 1
    else:
        print ("error")
print(red)
print(green)
print(black)
#print(counter)
plt.xlim(-1.5,0.5)
plt.ylim(-1,1.3)
plt.ylabel(r"log([OIII] $\lambda$5007/H$\beta$)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("BPT Diagram")
plt.scatter(Cloudy_NII_Ha,Cloudy_OIII_Hb,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
plt.legend([plt.scatter(NII_Ha[i],OIII_Hb[i],color=Low_Temp_Color, s = markersize), plt.scatter(NII_Ha[i],OIII_Hb[i],color=Mid_Temp_Color, s = markersize), plt.scatter(NII_Ha[i],OIII_Hb[i],color=High_Temp_Color, s = markersize),plt.scatter(Cloudy_NII_Ha,Cloudy_OIII_Hb,c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
x=np.linspace(-1.5,0.3,50)
y=((.61/(x-.47))+1.19)
plt.plot(x,y,color=Low_Temp_Color)
x3=np.linspace(-1,-0.2,50)
y3=((.61/(x3-.05)+1.3))
plt.plot(x3,y3,linestyle='--',color='red')
counter=0
sp2 = plt.subplot(222)
for i in range(0,len(SDSS_Data['z'])):
    if OIII_5006[i]/OIII_4363[i]<50:
        hot = plt.scatter(NII_Ha[i],Temp_Ratio[i],color=High_Temp_Color, s = markersize)
        counter=counter+1
        #print ("madeit")
    elif OIII_5006[i]/OIII_4363[i]>50 and OIII_5006[i]/OIII_4363[i]<100:
        mid_temp = plt.scatter(NII_Ha[i],Temp_Ratio[i],color=Mid_Temp_Color, s = markersize)
        counter=counter+1
        #print("k")
    elif OIII_5006[i]/OIII_4363[i]>100:
        cool = plt.scatter(NII_Ha[i],Temp_Ratio[i],color=Low_Temp_Color, s = markersize)
        counter=counter+1
        #print ("r")
        #print(Temp_Ratio[i])
        #print(OIII_4363[i])
    else:
        print ("error")
#print(counter)
plt.ylabel(r"log([OIII] $\lambda$5007/4363)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("Temperature vs. Ionization")
plt.ylim(0,3)
sim = plt.scatter(Cloudy_NII_Ha,Cloudy_Temp_Ratio,color=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')
plt.legend((cool, mid_temp, hot, sim), (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"), scatterpoints = 1, loc = 'lower left',fontsize =8)
sp3 = plt.subplot(223)
for i in range(0,len(SDSS_Data['z'])):
    if OIII_5006[i]/OIII_4363[i]<50:
        plt.scatter(NII_Ha[i],S_Ratio[i],color=High_Temp_Color, s = markersize)
        counter=counter+1
        #print ("madeit")
    elif OIII_5006[i]/OIII_4363[i]>50 and OIII_5006[i]/OIII_4363[i]<100:
        plt.scatter(NII_Ha[i],S_Ratio[i],color=Mid_Temp_Color, s = markersize)
        counter=counter+1
        #print("k")
    elif OIII_5006[i]/OIII_4363[i]>100:
        plt.scatter(NII_Ha[i],S_Ratio[i],color=Low_Temp_Color, s = markersize)
        counter=counter+1
        #print ("r")
        #print(Temp_Ratio[i])
        #print(OIII_4363[i])
    else:
        print ("error")
plt.ylabel(r"log([SII] $\lambda$6717/6731)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.ylim(-1.0,1.0)
plt.title("Density vs. Ionization")
plt.scatter(Cloudy_NII_Ha,Cloudy_S_Ratio,color=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')
plt.legend([plt.scatter(NII_Ha[i],Temp_Ratio[i],color=Low_Temp_Color, s = markersize), plt.scatter(NII_Ha[i],Temp_Ratio[i],color=Mid_Temp_Color, s = markersize), plt.scatter(NII_Ha[i],Temp_Ratio[i],color=High_Temp_Color, s = markersize),plt.scatter(Cloudy_NII_Ha,Cloudy_S_Ratio,color=Cloudy_Sim_Color, s = markersize)], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
sp4 = plt.subplot(224)
for i in range(0,len(SDSS_Data['z'])):
    if OIII_5006[i]/OIII_4363[i]<50:
        plt.scatter(NII_Ha[i],NO_Ratio[i],color=High_Temp_Color, s = markersize)
        counter=counter+1
        #print ("madeit")
    elif OIII_5006[i]/OIII_4363[i]>50 and OIII_5006[i]/OIII_4363[i]<100:
        plt.scatter(NII_Ha[i],NO_Ratio[i],color=Mid_Temp_Color, s = markersize)
        counter=counter+1
        #print("k")
    elif OIII_5006[i]/OIII_4363[i]>100:
        plt.scatter(NII_Ha[i],NO_Ratio[i],color=Low_Temp_Color, s = markersize)
        counter=counter+1
        #print ("r")
        #print(Temp_Ratio[i])
        #print(OIII_4363[i])
    else:
        print ("error")
plt.ylabel(r"log([NII] $\lambda$6584/[OII] $\lambda$3727)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("Metallicity vs. Ionization")
plt.scatter(Cloudy_NII_Ha,Cloudy_NO_Ratio,color=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')
plt.legend([plt.scatter(NII_Ha[i],NO_Ratio[i],color=Low_Temp_Color, s = markersize), plt.scatter(NII_Ha[i],NO_Ratio[i],color=Mid_Temp_Color, s = markersize), plt.scatter(NII_Ha[i],NO_Ratio[i],color=High_Temp_Color, s = markersize),plt.scatter(Cloudy_NII_Ha,Cloudy_NO_Ratio,color=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'upper left',fontsize =8)
plt.show()
print (SII_6716/SII_6731)
print (OII_3727)
#plt.savefig("Subplots 4363 > 5")
