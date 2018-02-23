import matplotlib.pyplot as plt
import numpy as np
import urllib
import matplotlib.cm as cm
import csv
Overview_File = '/Users/Sam/Documents/emgtemp/Studies/Basic_Grain5_Sim.ovr'
Overview_Data = np.genfromtxt(Overview_File, skip_header = 0, delimiter = '\t',dtype=float,names=True)
Heating_File = '/Users/Sam/Documents/emgtemp/Studies/Basic_Grain5_Sim.het'
Heating_Data_Floats = np.genfromtxt(Heating_File, skip_header = 0, delimiter = '\t',dtype=float,names=True)
Heating_Data_Strings = np.genfromtxt(Heating_File, skip_header = 0, delimiter = '\t',dtype=None,names=True)
heating_file = open('/Users/Sam/Documents/emgtemp/Studies/Basic_Grain5_Sim.het')
reader = csv.reader(heating_file)
##Heating_File = '/Users/compastro/jenkins/Metals_Sims/Baseline_Sim.het'
##Heating_Data = np.genfromtxt(Heating_File, skip_header = 0, delimiter = '\t',dtype=float,names=True)
temp = Overview_Data['Te']
depth = Overview_Data['depth']
O3 = Overview_Data['O3']
O2 = Overview_Data['O2']
O1 = Overview_Data['O1']
O4 = Overview_Data['O4']
H1 = Overview_Data['HI']
H2 = Overview_Data['HII']
hden = Overview_Data['hden']
eden = Overview_Data['eden']
heat_fracs1 = Heating_Data_Strings['heat_fracs1']
heat_val1 = Heating_Data_Floats['values1']
heat_fracs2 = Heating_Data_Strings['heat_fracs2']
heat_val2 = Heating_Data_Floats['values2']
heat_fracs3 = Heating_Data_Strings['heat_fracs3']
heat_val3 = Heating_Data_Floats['values3']
heat_fracs4 = Heating_Data_Strings['heat_fracs4']
heat_val4 = Heating_Data_Floats['values4']
##primary_heating = Heating_Data[]
O3_Ratio = 10**(O3)
O2_Ratio = 10**(O2)
O1_Ratio = 10**(O1)
O4_Ratio = 10**(O4)
H1_Ratio = 10**(H1)
H2_Ratio = 10**(H2)
##print(OIII)
#print(O3_Ratio)


fig = plt.figure(233)
fig.subplots_adjust(wspace=0,hspace=0)
sp1 = plt.subplot(311)
teline, = plt.plot(depth, temp, label = "T$_e$", c = 'blue')
hdenline, = plt.plot(depth, hden, label = "n$_h$", c = 'k')
edenline, = plt.plot(depth, eden, label = "n$_e$", c = 'r')
plt.legend(handles = [teline, hdenline, edenline], loc = 1)
plt.xlabel("Depth")
plt.ylabel("T$_e$")
#plt.title("T$_e$ vs. Depth")

sp2 = plt.subplot(312)
line1, = plt.plot(depth,O1_Ratio, linestyle = '--', c = 'k', label = "OI")
line2, = plt.plot(depth,O2_Ratio, c = 'k', label = "OII")
line3, = plt.plot(depth,O3_Ratio, linestyle = ":", linewidth = 2, c = 'k', label = "OIII")
line4, = plt.plot(depth,O4_Ratio, linestyle = '-.', linewidth = 1.5, c = 'k', label = "OIV")
line5, = plt.plot(depth,H1_Ratio, linestyle = '--', c = 'r', label = "HI")
line6, = plt.plot(depth,H2_Ratio, c = 'r', label = "HII")
plt.xlabel("Depth")
plt.ylabel("O Frac")
#plt.title("O Frac vs. Depth")
plt.legend(handles = [line1, line2, line3, line4, line5, line6], loc = 1)
print(H1_Ratio)
#print(O1_Ratio)

sp3 = plt.subplot(313)
He2_Array = np.zeros(len(hden))
H1_Array = np.zeros(len(hden))
He1_Array = np.zeros(len(hden))
Grain_Array = np.zeros(len(hden))
for i in range (0,len(hden)):
    #print (heat_fracs1[i])
    if heat_fracs1[i] == "He 2":
        He2_Array[i] = heat_val1[i]
    elif heat_fracs1[i] == "H  1":
        H1_Array[i] = heat_val1[i]
    elif heat_fracs1[i] == "He 1":
        He1_Array[i] = heat_val1[i]
    elif heat_fracs1[i] == "GrnP":
        Grain_Array[i] = heat_val1[i]
    if heat_fracs2[i] == "He 2":
        He2_Array[i] = heat_val2[i]
    elif heat_fracs2[i] == "H  1":
        H1_Array[i] = heat_val2[i]
    elif heat_fracs2[i] == "He 1":
        He1_Array[i] = heat_val2[i]
    elif heat_fracs2[i] == "GrnP":
        Grain_Array[i] = heat_val2[i]
    if heat_fracs3[i] == "He 2":
        He2_Array[i] = heat_val3[i]
    elif heat_fracs3[i] == "H  1":
        H1_Array[i] = heat_val3[i]
    elif heat_fracs3[i] == "He 1":
        He1_Array[i] = heat_val3[i]
    elif heat_fracs3[i] == "GrnP":
        Grain_Array[i] = heat_val3[i]
    if heat_fracs4[i] == "He 2":
        He2_Array[i] = heat_val4[i]
    elif heat_fracs4[i] == "H  1":
        H1_Array[i] = heat_val4[i]
    elif heat_fracs4[i] == "He 1":
        He1_Array[i] = heat_val4[i]
    elif heat_fracs4[i] == "GrnP":
        Grain_Array[i] = heat_val4[i]
#print(He2_Array)
line8, = plt.plot(depth,H1_Array, c = 'k', label = "HI")
line9, = plt.plot(depth,He1_Array, c = 'r', label = "HeI")
line10, = plt.plot(depth,He2_Array, c = 'b', label = "HeII")
line11, = plt.plot(depth,Grain_Array, c = 'g', label = "Grains")
plt.xlabel("Depth")
plt.ylabel("Heating Fraction")
#plt.title("Heating Fraction vs. Depth")
plt.legend(handles = [line8, line9, line10, line11], loc = 1)
plt.suptitle("Grains 5x")
#print(Grain_Array)
#print(H1_Array)
##for row in reader:
  ##  if row == "He 1"
 #   he2 = primary_heating
#elif primary_heating = "H  1":
 #   h1 = primary_heating
#elif primary_heating = "He 1":
 #   he1 = primary_heating
#if second_heating = "H  1":
#plt.plot()
plt.show()
plt.savefig("Basic Grains Heating Study.pdf")