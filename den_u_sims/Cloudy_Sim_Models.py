import matplotlib.pyplot as plt
import numpy as np
import urllib
import matplotlib.cm as cm
Low_Temp_Color = 'k'
Mid_Temp_Color = 'g'
High_Temp_Color = 'r'
Temp_Color = 0.5
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
Cloudy_File = '/Users/compastro/cloudy/c13.03/runs/Complete_Sim1.csv'
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
Cloudy_OI_Ratio = np.log10(Cloudy_OI_6300/Cloudy_Ha_6562)
Cloudy_O_Ratio = np.log10(Cloudy_OIII_5006/Cloudy_OII_3727)
Cloudy_S_Ha_Ratio = np.log10((Cloudy_SII_6716+Cloudy_SII_6731)/Cloudy_Ha_6562)
Grid_File = '/Users/compastro/cloudy/c13.03/runs/Complete_Sim1_Grid.csv'
Grid_Data = np.genfromtxt(Grid_File,skip_header=1,delimiter = ',',dtype=float,unpack=True)
Cloudy_U = Grid_Data[6,:]
Cloudy_Den = Grid_Data[7,:]
Cloudy_NII_Ha_array = np.reshape(Cloudy_NII_Ha,(3,-1))
Cloudy_OI_Ratio_array = np.reshape(Cloudy_OI_Ratio,(3,-1))
Cloudy_OIII_Hb_array = np.reshape(Cloudy_OIII_Hb,(3,-1))
Cloudy_Temp_Ratio_array = np.reshape(Cloudy_Temp_Ratio,(3,-1))
Cloudy_S_Ratio_array = np.reshape(Cloudy_S_Ratio,(3,-1))
Cloudy_NO_Ratio_array = np.reshape(Cloudy_NO_Ratio,(3,-1))
Cloudy_O_Ratio_array = np.reshape(Cloudy_O_Ratio,(3,-1))
Cloudy_S_Ha_Ratio_array = np.reshape(Cloudy_S_Ha_Ratio,(3,-1))
Cloudy_NII_Ha_transpose = np.transpose(Cloudy_NII_Ha_array)
Cloudy_OI_Ratio_transpose = np.transpose(Cloudy_OI_Ratio_array)
Cloudy_OIII_Hb_transpose = np.transpose(Cloudy_OIII_Hb_array)
Cloudy_Temp_Ratio_transpose = np.transpose(Cloudy_Temp_Ratio_array)
Cloudy_S_Ratio_transpose = np.transpose(Cloudy_S_Ratio_array)
Cloudy_NO_Ratio_transpose = np.transpose(Cloudy_NO_Ratio_array)
Cloudy_O_Ratio_transpose = np.transpose(Cloudy_O_Ratio_array)
Cloudy_S_Ha_Ratio_transpose = np.transpose(Cloudy_S_Ha_Ratio_array)
hden_colors = [plt.cm.Greys(i) for i in np.linspace(0,3)]
u_colors = [plt.cm.Reds(i) for i in np.linspace(0,7)]
#sf_count = 0.0
#comp_count = 0.0
#agn_count = 0.0
#liner_count = 0.0
#amb_count = 0.0
shape = ['v']

#####################################################################################################
def getShape(NII_Ha, OIII_Hb, S_Ha_Ratio, OI_Ratio):
    # Star forming
    if OIII_Hb < 0.61/(NII_Ha-0.05)+1.3 and \
    OIII_Hb < 0.72/(S_Ha_Ratio-0.32)+1.30 and \
    OIII_Hb < 0.73/(OI_Ratio+0.59)+1.33:
            shape = 'x'
            #sf_count = sf_count+1

    # Composite
    elif 0.61/(NII_Ha-0.05)+1.3 < OIII_Hb and \
    0.61/(NII_Ha-0.47)+1.19 > OIII_Hb:
            shape = '+'
            #comp_count = comp_count+1

    # AGN
    elif 0.61/(NII_Ha-0.47)+1.19 < OIII_Hb and \
    0.72/(S_Ha_Ratio-0.32)+1.30 < OIII_Hb and \
    0.73/(OI_Ratio+0.59)+1.33 < OIII_Hb and \
    (1.89*S_Ha_Ratio)+0.76 < OIII_Hb and \
    (1.18*OI_Ratio)+1.30 < OIII_Hb:
            shape = 'D'
            #agn_count = agn_count+1

    # LINERs
    elif 0.61/(NII_Ha-0.47)+1.19 < OIII_Hb and \
    0.72/(S_Ha_Ratio-0.32)+1.30 < OIII_Hb and \
    OIII_Hb < (1.89*S_Ha_Ratio)+0.76 and \
    0.73/(OI_Ratio+0.59)+1.33 < OIII_Hb and \
    OIII_Hb < (1.18*OI_Ratio)+1.30:
            shape = 's'
            #liner_count = liner_count+1

    else:
        # Ambiguous
            shape = '*'
            #amb_count = amb_count+1
            
    return shape
#####################################################################################################    
   
#####################################################################################################
def getColor(OIII_5006, OIII_4363):
        Temp_Color = 'k'
        if OIII_5006/OIII_4363<50:
            Temp_Color = 'r'
        #red = red + 1 
        elif OIII_5006/OIII_4363>50 and OIII_5006/OIII_4363<100:
            Temp_Color = 'g'
        #green = green + 1
        elif OIII_5006/OIII_4363>100:
            Temp_Color = 'k'
        #black = black + 1
        else:
            print ("error")
        return Temp_Color
#####################################################################################################

fig = plt.figure(2)
fig.subplots_adjust(wspace=0.4,hspace=0.4)
sp1 = plt.subplot(221)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],OIII_Hb[i],color=Temp_Color, s = markersize, marker = shape)
    
    
#print(sf_count)
#print(comp_count)
#print(agn_count)
#print(liner_count)
#print(amb_count)
#print(red)
#print(green)
#print(black)
#print(counter)
plt.xlim(-2.5,0.5)
plt.ylim(-1,1.3)
plt.ylabel(r"log([OIII] $\lambda$5007/H$\beta$)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("BPT Diagram")
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_NII_Ha[i],Cloudy_OIII_Hb[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_NII_Ha_array[i],Cloudy_OIII_Hb_array[i],c='u_colors')
        plt.plot(Cloudy_NII_Ha_transpose[i],Cloudy_OIII_Hb_transpose[i],c='hden_colors')
plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
x=np.linspace(-1.5,0.3,50)
y=((.61/(x-.47))+1.19)
plt.plot(x,y,color=Low_Temp_Color)
x3=np.linspace(-1,-0.2,50)
y3=((.61/(x3-.05)+1.3))
plt.plot(x3,y3,linestyle='--',color='red')
#counter=0




sp2 = plt.subplot(222)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],Temp_Ratio[i],color=Temp_Color, s = markersize, marker = shape)
#print(counter)
plt.ylabel(r"log([OIII] $\lambda$5007/4363)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("Temperature vs. Ionization")
plt.ylim(0,3)
plt.xlim(-2.5,0.5)
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        sim = plt.scatter(Cloudy_NII_Ha[i],Cloudy_Temp_Ratio[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_NII_Ha_array[i],Cloudy_Temp_Ratio_array[i],c='k')
        plt.plot(Cloudy_NII_Ha_transpose[i],Cloudy_Temp_Ratio_transpose[i],c='k')
plt.legend([plt.scatter([],[],color='.75', s = markersize, marker = 'x', edgecolor = 'none'),plt.scatter([],[],color='0.75', s = markersize, marker = '+', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 'D', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 's', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = '*', edgecolor = 'none')], ("Star-Forming","Composite","AGN","LINER","Ambiguous"),scatterpoints = 1, loc = 'lower left',fontsize =8)


sp3 = plt.subplot(223)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],S_Ratio[i],color=Temp_Color, s = markersize, marker = shape)
plt.ylabel(r"log([SII] $\lambda$6717/6731)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.ylim(-1.0,1.0)
plt.xlim(-2.5,0.5)
plt.title("Density vs. Ionization")
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_NII_Ha[i],Cloudy_S_Ratio[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_NII_Ha_array[i],Cloudy_S_Ratio_array[i],c='k')
        plt.plot(Cloudy_NII_Ha_transpose[i],Cloudy_S_Ratio_transpose[i],c='k')
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)



sp4 = plt.subplot(224)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],NO_Ratio[i],color=Temp_Color, s = markersize, marker = shape)
plt.ylabel(r"log([NII] $\lambda$6584/[OII] $\lambda$3727)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("Metallicity vs. Ionization")
plt.xlim(-2.5,0.5)
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_NII_Ha[i],Cloudy_NO_Ratio[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_NII_Ha_array[i],Cloudy_NO_Ratio_array[i],c='k')
        plt.plot(Cloudy_NII_Ha_transpose[i],Cloudy_NO_Ratio_transpose[i],c='k')
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
plt.show()
plt.suptitle('2 < hden < 4, -0.5 < U < -3.5, Z = 1.2, T = 5.3, a(ox) = -1.42, a(uv) = -0.57, a(x) = -1.63')
#plt.savefig("First Complete Sim Plots.pdf")




fig2 = plt.figure(3)
sp5 = plt.subplot(221)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],OI_Ratio[i],color=Temp_Color, s = markersize, marker = shape)
plt.ylabel(r"log([OI] $\lambda$6300/H$\alpha$)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("OI_6300 vs. Ionization")
plt.xlim(-2.5,0.5)
plt.ylim(-2.5,0)
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_NII_Ha[i],Cloudy_OI_Ratio[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_NII_Ha_array[i],Cloudy_OI_Ratio_array[i],c='k')
        plt.plot(Cloudy_NII_Ha_transpose[i],Cloudy_OI_Ratio_transpose[i], c='k')
#for i in range (0,len(Cloudy_Data['O__1__6300A'])-1)    
plt.legend([plt.scatter([],[],color='.75', s = markersize, marker = 'x', edgecolor = 'none'),plt.scatter([],[],color='0.75', s = markersize, marker = '+', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 'D', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 's', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = '*', edgecolor = 'none')], ("Star-Forming","Composite","AGN","LINER","Ambiguous"),scatterpoints = 1, loc = 'lower left',fontsize =8)




sp6 = plt.subplot(222)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(OI_Ratio[i],OIII_Hb[i],color=Temp_Color, s = markersize, marker = shape)
plt.ylabel(r"log([OIII] $\lambda$5007/H$\beta$)")
plt.xlabel(r"log ([OI] $\lambda$6300/H$\alpha$)")
plt.title("OI_6300 vs. OIII_5007")
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_OI_Ratio[i],Cloudy_OIII_Hb[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_OI_Ratio_array[i],Cloudy_OIII_Hb_array[i],c='k')
        plt.plot(Cloudy_OI_Ratio_transpose[i],Cloudy_OIII_Hb_transpose[i],c='k')
x6 = np.linspace(-2.5,-0.6,50)
y6 = ((.73/(x6+0.59))+1.33)
plt.plot(x6,y6,color = 'k')
x7 = np.linspace(-1.125,0.25,50)
y7 = (1.18*x7) + 1.30
plt.plot(x7,y7, color = 'b')
plt.ylim(-1,1.5)
plt.xlim(-2.5,0.5)
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)



sp7 = plt.subplot(223)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(OI_Ratio[i],O_Ratio[i],color=Temp_Color, s = markersize, marker = shape)
plt.ylabel(r"log([OIII] $\lambda$5007/[OII]$\lambda$3727)")
plt.xlabel(r"log ([OI] $\lambda$6300/H$\alpha$)")
plt.title("Groves Diagram")
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_OI_Ratio[i],Cloudy_O_Ratio[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_OI_Ratio_array[i],Cloudy_O_Ratio_array[i],c='k')
        plt.plot(Cloudy_OI_Ratio_transpose[i],Cloudy_O_Ratio_transpose[i],c='k')
x1 = np.linspace(-2.0,-.25,50)
y1 = ((-1.701*x1)-2.163)
x2 = np.linspace(-1.05998,0,50)
y2 = x2 + 0.7
plt.plot(x2,y2, color = 'k')
plt.plot(x1,y1, color = 'k')
plt.xlim(-2.5,0)
plt.ylim(-1.5,1)
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)



sp8 = plt.subplot(224)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(S_Ha_Ratio[i],OIII_Hb[i],color=Temp_Color, s = markersize, marker = shape)
plt.ylabel(r"log([OIII] $\lambda$5007/H$\beta$)")
plt.xlabel(r"log ([SII]/H$\alpha$)")
plt.title("OIII_5007 vs. SII")
plt.ylim(-1,1.5)
x4 = np.linspace(-0.32,0.25,50)
y4 = ((1.89*x4)+0.76)
x5 = np.linspace(-1.5,0.25,50)
y5 = ((0.72/(x - 0.32))+1.3)
plt.plot(x5,y5,color = 'k')
plt.plot(x4,y4,color = 'b')
for i in range(0,len(Cloudy_Data['O__1__6300A'])):
        plt.scatter(Cloudy_S_Ha_Ratio[i],Cloudy_OIII_Hb[i],c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
for i in range(0,len(Cloudy_NII_Ha_array)):
        plt.plot(Cloudy_S_Ha_Ratio_array[i],Cloudy_OIII_Hb_array[i],c='k')
        plt.plot(Cloudy_S_Ha_Ratio_transpose[i],Cloudy_OIII_Hb_transpose[i],c='k')
plt.suptitle('2 < hden < 4, -0.5 < U < -3.5, Z = 1.2, T = 5.3, a(ox) = -1.42, a(uv) = -0.57, a(x) = -1.63')
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
#plt.savefig("First Complete Sim Plots1.pdf")