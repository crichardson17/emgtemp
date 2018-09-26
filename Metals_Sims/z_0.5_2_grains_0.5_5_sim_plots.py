import matplotlib.pyplot as plt
import numpy as np
import urllib
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cm as cm
Low_Temp_Color = 'k'
Mid_Temp_Color = 'g'
High_Temp_Color = 'r'
#Temp_Color = 0.5
Cloudy_Sim_Color = 'cyan'
markersize = 40
SDSS_File = '/Users/Sam/Documents/emgtemp/data/4363_gr_5_0_err_dered.csv'
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
Cloudy_File = '/Users/Sam/Documents/emgtemp/Grains_Sims/z_0.5_2.0_grains_sims.pun'
Cloudy_Data = np.genfromtxt(Cloudy_File, delimiter = '\t',dtype=float,unpack=True,names=True)
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
Grid_File = '/Users/Sam/Documents/emgtemp/Metals_Sims/Dusty_sims/z_0.6_2.0/z_0.6_2.0_sims.grd'
Grid_Data = np.genfromtxt(Grid_File,skip_header=1,delimiter = '\t',dtype=float,unpack=True)
Cloudy_Metals = Grid_Data[8,:]
Cloudy_Den = Grid_Data[6,:]
Cloudy_NII_Ha_array = np.reshape(Cloudy_NII_Ha,(6,-1))
Cloudy_OI_Ratio_array = np.reshape(Cloudy_OI_Ratio,(6,-1))
Cloudy_OIII_Hb_array = np.reshape(Cloudy_OIII_Hb,(6,-1))
Cloudy_Temp_Ratio_array = np.reshape(Cloudy_Temp_Ratio,(6,-1))
Cloudy_S_Ratio_array = np.reshape(Cloudy_S_Ratio,(6,-1))
Cloudy_NO_Ratio_array = np.reshape(Cloudy_NO_Ratio,(6,-1))
Cloudy_O_Ratio_array = np.reshape(Cloudy_O_Ratio,(6,-1))
Cloudy_S_Ha_Ratio_array = np.reshape(Cloudy_S_Ha_Ratio,(6,-1))
Cloudy_NII_Ha_transpose = np.transpose(Cloudy_NII_Ha_array)
Cloudy_OI_Ratio_transpose = np.transpose(Cloudy_OI_Ratio_array)
Cloudy_OIII_Hb_transpose = np.transpose(Cloudy_OIII_Hb_array)
Cloudy_Temp_Ratio_transpose = np.transpose(Cloudy_Temp_Ratio_array)
Cloudy_S_Ratio_transpose = np.transpose(Cloudy_S_Ratio_array)
Cloudy_NO_Ratio_transpose = np.transpose(Cloudy_NO_Ratio_array)
Cloudy_O_Ratio_transpose = np.transpose(Cloudy_O_Ratio_array)
Cloudy_S_Ha_Ratio_transpose = np.transpose(Cloudy_S_Ha_Ratio_array)
#cold_data_colors = [plt.cm.Blues(i) for i in np.linspace(0,1,len(SDSS_Data['z']))]
#mid_data_colors = [plt.cm.Greens(i) for i in np.linspace(0,1,len(SDSS_Data['z']))]
#hot_data_colors = [plt.cm.Reds(i) for i in np.linspace(0,1,len(SDSS_Data['z']))]
grains_colors = [plt.cm.Reds(i) for i in np.linspace(0.25,1,10)]
z_colors = [plt.cm.Blues(i) for i in np.linspace(0.25,1,6)]
def truncate_colormap(cmap, minval=0.15, maxval=1.0, n=100):
  	new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
  	return new_cmap
grains_colors_map = truncate_colormap(cm.Reds)
z_colors_map = truncate_colormap(cm.Blues)

#This is bad^ 3 and 7 are the number of densities and ionization parameters used, but ideally this wouldn't be hardcoded.

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
            #Temp_Color = '0.25'
            Temp_Color = plt.cm.gray(0.2)
        #red = red + 1 
        elif OIII_5006/OIII_4363>50 and OIII_5006/OIII_4363<100:
            #Temp_Color = '0.5'
            Temp_Color = plt.cm.gray(0.5)
        #green = green + 1
        elif OIII_5006/OIII_4363>100:
            #Temp_Color = '0.75'
            Temp_Color = plt.cm.gray(0.75)

        #black = black + 1
        else:
            print ("error")
        return Temp_Color
#####################################################################################################
agn_count = 0
liner_count = 0
sf_count = 0
comp_count = 0
amb_count = 0
fig = plt.figure(131)
fig.subplots_adjust(wspace=0.4,hspace=0.4)
sp1 = plt.subplot(221)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    if shape == 'D':
        agn_count = agn_count +1
    elif shape == 's':
        liner_count = liner_count +1
    elif shape == '+':
        comp_count = comp_count +1
    elif shape == 'x':
        sf_count = sf_count+1
    elif shape == '*':
        amb_count = amb_count +1
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    #print(Temp_Color)
    plt.scatter(NII_Ha[i],OIII_Hb[i],s = markersize, marker = shape, color = Temp_Color, edgecolor = 'none')
    #print (Temp_Color)    

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
#plt.scatter(Cloudy_NII_Ha,Cloudy_OIII_Hb,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp1.set_color_cycle(grains_colors)
plt.plot(Cloudy_NII_Ha_array,Cloudy_OIII_Hb_array, lw = '2')
sp1.set_color_cycle(z_colors)
plt.plot(Cloudy_NII_Ha_transpose,Cloudy_OIII_Hb_transpose, lw = '2',linestyle = '--')
plt.legend([plt.scatter([],[],color='.75', s = markersize, marker = 'x', edgecolor = 'none'),plt.scatter([],[],color='0.75', s = markersize, marker = '+', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 'D', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 's', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = '*', edgecolor = 'none')], ("Star-Forming","Composite","AGN","LINER","Ambiguous"),scatterpoints = 1, loc = 'lower left',fontsize =8)
x=np.linspace(-1.5,0.3,50)
y=((.61/(x-.47))+1.19)
plt.plot(x,y,color=Low_Temp_Color)
x3=np.linspace(-1,-0.2,50)
y3=((.61/(x3-.05)+1.3))
plt.plot(x3,y3,linestyle='--',color='k')
#counter=0

sm = plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0.5, vmax=5.0),cmap=grains_colors_map)
sm._A = []
smaxes = inset_axes(sp1, width=0.06, height=0.4, loc=3, bbox_to_anchor=(0.31, 0.1), bbox_transform=sp1.figure.transFigure)
#smaxes = inset_axes(sp1, width="3%", height="20%", loc=3, bbox_to_anchor=(0.1, 0.1), bbox_transform=ax.figure.transFigure)
cbar = plt.colorbar(sm,cax=smaxes)
cbar.ax.set_title('Grains',fontsize=8)
cbar.set_ticks([0.5,5.0])
cbar.set_ticklabels([0.5,5.0])
cbar.ax.tick_params(labelsize=8) 


sp2 = plt.subplot(222)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],Temp_Ratio[i], s = markersize, marker = shape, color = Temp_Color, edgecolor = 'none')
#print(counter)
plt.ylabel(r"log([OIII] $\lambda$5007/4363)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("Temperature")
plt.ylim(0,3)
plt.xlim(-2.5,0.5)
#plt.scatter(Cloudy_NII_Ha,Cloudy_Temp_Ratio,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp2.set_color_cycle(grains_colors)
plt.plot(Cloudy_NII_Ha_array,Cloudy_Temp_Ratio_array, lw = '2')
sp2.set_color_cycle(z_colors)
plt.plot(Cloudy_NII_Ha_transpose,Cloudy_Temp_Ratio_transpose, lw = '2',linestyle = '--')
plt.legend([plt.scatter([],[],color='0.75', s = markersize), plt.scatter([],[],color='0.5', s = markersize), plt.scatter([],[],color='0.25', s = markersize)], (r"T$_e$<1.17*10$^4$",r"1.17*10$^4$<T$_e$<1.54*10$^4$",r"T$_e$>1.54*10$^4$"),scatterpoints = 1, loc = 'lower left',fontsize =8)

sm = plt.cm.ScalarMappable(norm=colors.Normalize(vmin=0.5, vmax=2.0),cmap=z_colors_map)
sm._A = []
smaxes = inset_axes(sp2, width=0.06, height=0.4, loc=3, bbox_to_anchor=(0.14, .1), bbox_transform=sp2.figure.transFigure)
cbar = plt.colorbar(sm,cax=smaxes)
cbar.ax.set_title('Z',fontsize=8)
cbar.set_ticks([0.5,2.0])
cbar.set_ticklabels([0.5,2.0])
cbar.ax.tick_params(labelsize=8) 

sp3 = plt.subplot(223)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],S_Ratio[i], s = markersize, marker = shape, c = Temp_Color, edgecolor = 'none')
plt.ylabel(r"log([SII] $\lambda$6717/6731)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.ylim(-1.0,1.0)
plt.xlim(-2.5,0.5)
plt.title("Density")
#plt.scatter(Cloudy_NII_Ha,Cloudy_S_Ratio,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp3.set_color_cycle(grains_colors)
plt.plot(Cloudy_NII_Ha_array,Cloudy_S_Ratio_array, lw = '2')
sp3.set_color_cycle(z_colors)
plt.plot(Cloudy_NII_Ha_transpose,Cloudy_S_Ratio_transpose, lw = '2',linestyle = '--')
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)



sp4 = plt.subplot(224)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],NO_Ratio[i], s = markersize, marker = shape, c = Temp_Color, edgecolor = 'none')
plt.ylabel(r"log([NII] $\lambda$6584/[OII] $\lambda$3727)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("Metallicity")
plt.xlim(-2.5,0.5)
#plt.scatter(Cloudy_NII_Ha,Cloudy_NO_Ratio,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp4.set_color_cycle(grains_colors)
plt.plot(Cloudy_NII_Ha_array,Cloudy_NO_Ratio_array, lw = '2')
sp4.set_color_cycle(z_colors)
plt.plot(Cloudy_NII_Ha_transpose,Cloudy_NO_Ratio_transpose, lw = '2',linestyle = '--')
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
plt.show()
#plt.suptitle('hden = 2.4, U = -1.5, 0.5 < Z < 2.0, 0.5 < grains < 5.0')
plt.savefig("z_0.5_2_grains_0.5_5_plots.pdf")




fig2 = plt.figure(132)
sp5 = plt.subplot(221)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(NII_Ha[i],OI_Ratio[i], s = markersize, marker = shape, c = Temp_Color, edgecolor = 'none')
plt.ylabel(r"log([OI] $\lambda$6300/H$\alpha$)")
plt.xlabel(r"log ([NII] $\lambda$6584/H$\alpha$)")
plt.title("OI_6300")
plt.xlim(-2.5,0.5)
plt.ylim(-2.5,0)
#plt.scatter(Cloudy_NII_Ha,Cloudy_OI_Ratio,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp5.set_color_cycle(grains_colors)
plt.plot(Cloudy_NII_Ha_array,Cloudy_OI_Ratio_array, lw = '2')
sp5.set_color_cycle(z_colors)
plt.plot(Cloudy_NII_Ha_transpose,Cloudy_OI_Ratio_transpose, lw = '2',linestyle = '--')  
plt.legend([plt.scatter([],[],color='.75', s = markersize, marker = 'x', edgecolor = 'none'),plt.scatter([],[],color='0.75', s = markersize, marker = '+', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 'D', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = 's', edgecolor = 'none'), plt.scatter([],[],color='.75', s = markersize, marker = '*', edgecolor = 'none')], ("Star-Forming","Composite","AGN","LINER","Ambiguous"),scatterpoints = 1, loc = 'lower left',fontsize =8)




sp6 = plt.subplot(222)
for i in range(0,len(SDSS_Data['z'])):
    shape = getShape(NII_Ha[i], OIII_Hb[i], S_Ha_Ratio[i], OI_Ratio[i])
    Temp_Color = getColor(OIII_5006[i], OIII_4363[i])
    plt.scatter(OI_Ratio[i],OIII_Hb[i], s = markersize, marker = shape, c = Temp_Color, edgecolor = 'none')
plt.ylabel(r"log([OIII] $\lambda$5007/H$\beta$)")
plt.xlabel(r"log ([OI] $\lambda$6300/H$\alpha$)")
plt.title("OI_6300 vs. OIII_5007")
#plt.scatter(Cloudy_OI_Ratio,Cloudy_OIII_Hb,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp6.set_color_cycle(grains_colors)
plt.plot(Cloudy_OI_Ratio_array,Cloudy_OIII_Hb_array, lw = '2')
sp6.set_color_cycle(z_colors)
plt.plot(Cloudy_OI_Ratio_transpose,Cloudy_OIII_Hb_transpose, lw = '2',linestyle = '--')
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
    plt.scatter(OI_Ratio[i],O_Ratio[i], s = markersize, marker = shape, c = Temp_Color, edgecolor = 'none')
plt.ylabel(r"log([OIII] $\lambda$5007/[OII]$\lambda$3727)")
plt.xlabel(r"log ([OI] $\lambda$6300/H$\alpha$)")
plt.title("Groves Diagram")
#plt.scatter(Cloudy_OI_Ratio,Cloudy_O_Ratio,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp7.set_color_cycle(grains_colors)
plt.plot(Cloudy_OI_Ratio_array,Cloudy_O_Ratio_array, lw = '2')
sp7.set_color_cycle(z_colors)
plt.plot(Cloudy_OI_Ratio_transpose,Cloudy_O_Ratio_transpose, lw = '2',linestyle = '--')
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
    plt.scatter(S_Ha_Ratio[i],OIII_Hb[i], s = markersize, marker = shape, c = Temp_Color, edgecolor = 'none')
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
#plt.scatter(Cloudy_S_Ha_Ratio,Cloudy_OIII_Hb,c=Cloudy_Sim_Color, s = markersize, edgecolor ='none')
sp8.set_color_cycle(grains_colors)
plt.plot(Cloudy_S_Ha_Ratio_array,Cloudy_OIII_Hb_array, lw = '2')
sp8.set_color_cycle(z_colors)
plt.plot(Cloudy_S_Ha_Ratio_transpose,Cloudy_OIII_Hb_transpose, lw = '2',linestyle = '--')
plt.suptitle('hden = 2.4, U = -1.5, 0.5 < Z < 2.0, 0.5 < grains < 5.0')
#plt.legend([plt.scatter([],[],color=Low_Temp_Color, s = markersize), plt.scatter([],[],color=Mid_Temp_Color, s = markersize), plt.scatter([],[],color=High_Temp_Color, s = markersize),plt.scatter([],[],c=Cloudy_Sim_Color, s = markersize, edgecolor = 'none')], (r"$\frac{OIII[5007]}{OIII[4363]}$<50.0",r"$50.0<\frac{OIII[5007]}{OIII[4363]}<100.0$",r"$\frac{OIII[5007]}{OIII[4363]}$>100.0","Cloudy Simulation"),scatterpoints = 1, loc = 'lower left',fontsize =8)
#plt.savefig("Metallicity Sim Plots1.pdf")
plt.show()
#print(agn_count)
#print(liner_count)
#print(comp_count)
#print(sf_count)
#print (amb_count)
#Not sure why I'm only getting five metallicity lines. I should be getting six.