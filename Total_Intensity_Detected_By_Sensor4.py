# -*- coding: utf-8 -*-
"""
Created on Thu May 28 17:12:03 2020

@author: elham.rahmati
"""


# -*- coding: utf-8 -*-
"""
Created on Mon May 25 10:29:32 2020

@author: elham.rahmati
"""


# -*- coding: utf-8 -*-
"""
Created on Mon May 25 10:02:38 2020

@author: elham.rahmati

Total sun irradiance that is detected with GEN3 & GEN8 after passing the stack that contains
OLED display, IR filter,low-temp GRN_400-1100nm,SIC_1000, Angular filter
vs OLED display, PassBandKeihin,low-temp GRN_400-1100nm,SIC_1000, Angular filter and GEN3 & GEN8 sensors.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.integrate import simps
import xlsxwriter

#from sympy import pretty_print as pp, latex

#function for Reading different excel files
def ReadFromExcel(file_url, rowmin, rowmax, SheetName, Columns, flip):
    df = pd.read_excel(file_url,sheet_name=SheetName, skiprows=range(1,rowmin-1), usecols=Columns)
    data_list = df.values.tolist()
    data_list_part=data_list[0:rowmax+1-rowmin]
    data0=np.asarray(data_list_part)
    if (flip):
        data =np.flipud(data0)
    else:
        data=data0
    return data

#main
   
############      read data from excel     ###########
    
folder=r'C:\Users\elham.rahmati\Documents\AFMM\Aindex=50_1.01228Defgree_Inpstruct_600\HOmeWork\Excel'

#Harmonics 301 and coating thickness 0 Incident angle 0
file1=folder+'\Dataset_0_H301_h0.csv'
file1_sh1=ReadFromExcel(file1, 1,125, 'Dataset_0_H301_h0', 'A,B',False)

#Harmonics 301 and coating thickness 1 Incident angle 0
file2=folder+'\Dataset_0_H301_h1.csv'
file2_sh1=ReadFromExcel(file2, 1,130, 'Dataset_0_H301_h1', 'A,B',False)

#Harmonics 301 and coating thickness 20 Incident angle 0
file3=folder+'\Dataset_0_H301_h20.csv'
file3_sh1=ReadFromExcel(file3, 1,130, 'Dataset_0_H301_h20', 'A,B',False)

#Harmonics 301 and coating thickness 250 Incident angle 0
file4=folder+'\Dataset_0_H301_h250.csv'
file4_sh1=ReadFromExcel(file4, 1,127, 'Dataset_0_H301_h250', 'A,B',False)

#Harmonics 301 and coating thickness 350 Incident angle 0
file5=folder+'\Dataset_0_H301_h350.csv'
file5_sh1=ReadFromExcel(file5, 1,127, 'Dataset_0_H301_h350', 'A,B',False)

#Harmonics 401 and coating thickness 0 Incident angle 0
file6=folder+'\Dataset_0_H401_h0.csv'
file6_sh1=ReadFromExcel(file6, 1,128, 'Dataset_0_H401_h0', 'A,B',False)

#Harmonics 401 and coating thickness 1 Incident angle 0
file7=folder+'\Dataset_0_H401_h1.csv'
file7_sh1=ReadFromExcel(file7, 1,126, 'Dataset_0_H401_h1', 'A,B',False)

#Harmonics 401 and coating thickness 20 Incident angle 0
file8=folder+'\Dataset_0_H401_h20.csv'
file8_sh1=ReadFromExcel(file8, 1,127, 'Dataset_0_H401_h20', 'A,B',False)

#Harmonics 401 and coating thickness 250 Incident angle 0
file9=folder+'\Dataset_0_H401_h250.csv'
file9_sh1=ReadFromExcel(file9, 1,128, 'Dataset_0_H401_h250', 'A,B',False)

#Harmonics 401 and coating thickness 350 Incident angle 0
file10=folder+'\Dataset_0_H401_h350.csv'
file10_sh1=ReadFromExcel(file10, 1,128, 'Dataset_0_H401_h350', 'A,B',False)

#Harmonics 501 and coating thickness 0 Incident angle 0
file11=folder+'\Dataset_0_H501_h0.csv'
file11_sh1=ReadFromExcel(file11, 1,127, 'Dataset_0_H501_h0', 'A,B',False)

#Harmonics 501 and coating thickness 1 Incident angle 0
file12=folder+'\Dataset_0_H501_h1.csv'
file12_sh1=ReadFromExcel(file12, 1,128, 'Dataset_0_H501_h1', 'A,B',False)

#Harmonics 501 and coating thickness 20 Incident angle 0
file13=folder+'\Dataset_0_H501_h20.csv'
file13_sh1=ReadFromExcel(file13, 1,132, 'Dataset_0_H501_h20', 'A,B',False)

#Harmonics 501 and coating thickness 250 Incident angle 0
file14=folder+'\Dataset_0_H501_h250.csv'
file14_sh1=ReadFromExcel(file14, 1,129, 'Dataset_0_H501_h250', 'A,B',False)

#Harmonics 501 and coating thickness 350 Incident angle 0
file15=folder+'\Dataset_0_H501_h350.csv'
file15_sh1=ReadFromExcel(file15, 1,129, 'Dataset_0_H501_h350', 'A,B',False)

#Harmonics 301 and coating thickness 0 Incident angle 1.0123
file16=folder+'\Dataset_H301_h0.csv'
file16_sh1=ReadFromExcel(file16, 1,123, 'Dataset_H301_h0', 'A,B',False)

#Harmonics 301 and coating thickness 1 Incident angle 1.0123
file17=folder+'\Dataset_H301_h1.csv'
file17_sh1=ReadFromExcel(file17, 1,124, 'Dataset_H301_h1', 'A,B',False)

#Harmonics 301 and coating thickness 20 Incident angle 1.0123
file18=folder+'\Dataset_H301_h20.csv'
file18_sh1=ReadFromExcel(file18, 1,123, 'Dataset_H301_h20', 'A,B',False)

#Harmonics 301 and coating thickness 250 Incident angle 1.0123
file19=folder+'\Dataset_H301_h250.csv'
file19_sh1=ReadFromExcel(file19, 1,122, 'Dataset_H301_h250', 'A,B',False)

#Harmonics 301 and coating thickness 350 Incident angle 1.0123
file20=folder+'\Dataset_H301_h350.csv'
file20_sh1=ReadFromExcel(file20, 1,124, 'Dataset_H301_h350', 'A,B',False)

#Harmonics 401 and coating thickness 0 Incident angle 1.0123
file21=folder+'\Dataset_H401_h0.csv'
file21_sh1=ReadFromExcel(file21, 1,123, 'Dataset_H401_h0', 'A,B',False)

#Harmonics 401 and coating thickness 1 Incident angle 1.0123
file22=folder+'\Dataset_H401_h1.csv'
file22_sh1=ReadFromExcel(file22, 1,121, 'Dataset_H401_h1', 'A,B',False)

#Harmonics 401 and coating thickness 20 Incident angle 1.0123
file23=folder+'\Dataset_H401_h20.csv'
file23_sh1=ReadFromExcel(file23, 1,120, 'Dataset_H401_h20', 'A,B',False)

#Harmonics 401 and coating thickness 250 Incident angle 1.0123
file24=folder+'\Dataset_H401_h250.csv'
file24_sh1=ReadFromExcel(file24, 1,121, 'Dataset_H401_h250', 'A,B',False)

#Harmonics 401 and coating thickness 350 Incident angle 1.0123
file25=folder+'\Dataset_H401_h350.csv'
file25_sh1=ReadFromExcel(file25, 1,121, 'Dataset_H401_h350', 'A,B',False)

#Harmonics 501 and coating thickness 0 Incident angle 1.0123
file26=folder+'\Dataset_H501_h0.csv'
file26_sh1=ReadFromExcel(file26, 1,123, 'Dataset_H501_h0', 'A,B',False)

#Harmonics 501 and coating thickness 1 Incident angle 1.0123
file27=folder+'\Dataset_H501_h1.csv'
file27_sh1=ReadFromExcel(file27, 1,121, 'Dataset_H501_h1', 'A,B',False)

#Harmonics 501 and coating thickness 20 Incident angle 1.0123
file28=folder+'\Dataset_H501_h20.csv'
file28_sh1=ReadFromExcel(file28, 1,123, 'Dataset_H501_h20', 'A,B',False)

#Harmonics 501 and coating thickness 250 Incident angle 1.0123
file29=folder+'\Dataset_H501_h250.csv'
file29_sh1=ReadFromExcel(file29, 1,123, 'Dataset_H501_h250', 'A,B',False)

#Harmonics 501 and coating thickness 350 Incident angle 1.0123
file30=folder+'\Dataset_H501_h350.csv'
file30_sh1=ReadFromExcel(file30, 1,123, 'Dataset_H501_h350', 'A,B',False)

################################ Interpolation #####################################
XMin= -0.00083  #is valid for all data
XMax=0.00083  #is valid for all data
dX=0.00001
XArray = np.arange(XMin,XMax, dX)

#Harmonics 301 and coating thickness 0 Incident angle 0
file1_sh1_interp = interpolate.splrep(file1_sh1[:,0],file1_sh1[:,1],s=0)
H301_0_h0 = interpolate.splev(XArray, file1_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 1 Incident angle 0
file2_sh1_interp = interpolate.splrep(file2_sh1[:,0],file2_sh1[:,1],s=0)
H301_0_h1 = interpolate.splev(XArray, file2_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 20 Incident angle 0
file3_sh1_interp = interpolate.splrep(file3_sh1[:,0],file3_sh1[:,1],s=0)
H301_0_h20 = interpolate.splev(XArray, file3_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 250 Incident angle 0
file4_sh1_interp = interpolate.splrep(file4_sh1[:,0],file4_sh1[:,1],s=0)
H301_0_h250 = interpolate.splev(XArray, file4_sh1_interp, der=0,ext=3) #IR filter sheet 34

#Harmonics 301 and coating thickness 350 Incident angle 0
file5_sh1_interp = interpolate.splrep(file5_sh1[:,0],file5_sh1[:,1],s=0)
H301_0_h350 = interpolate.splev(XArray, file5_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 0 Incident angle 0
file6_sh1_interp = interpolate.splrep(file6_sh1[:,0],file6_sh1[:,1],s=0)
H401_0_h0 = interpolate.splev(XArray, file6_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 1 Incident angle 0
file7_sh1_interp = interpolate.splrep(file7_sh1[:,0],file7_sh1[:,1],s=0)
H401_0_h1 = interpolate.splev(XArray, file7_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 20 Incident angle 0
file8_sh1_interp = interpolate.splrep(file8_sh1[:,0],file8_sh1[:,1],s=0)
H401_0_h20 = interpolate.splev(XArray, file8_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 250 Incident angle 0
file9_sh1_interp = interpolate.splrep(file9_sh1[:,0],file9_sh1[:,1],s=0)
H401_0_h250 = interpolate.splev(XArray, file9_sh1_interp, der=0,ext=3) #IR filter sheet 34

#Harmonics 401 and coating thickness 350 Incident angle 0
file10_sh1_interp = interpolate.splrep(file10_sh1[:,0],file10_sh1[:,1],s=0)
H401_0_h350 = interpolate.splev(XArray, file10_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 0 Incident angle 0
file11_sh1_interp = interpolate.splrep(file11_sh1[:,0],file11_sh1[:,1],s=0)
H501_0_h0 = interpolate.splev(XArray, file11_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 1 Incident angle 0
file12_sh1_interp = interpolate.splrep(file12_sh1[:,0],file12_sh1[:,1],s=0)
H501_0_h1 = interpolate.splev(XArray, file12_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 20 Incident angle 0
file13_sh1_interp = interpolate.splrep(file13_sh1[:,0],file13_sh1[:,1],s=0)
H501_0_h20 = interpolate.splev(XArray, file13_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 250 Incident angle 0
file14_sh1_interp = interpolate.splrep(file14_sh1[:,0],file14_sh1[:,1],s=0)
H501_0_h250 = interpolate.splev(XArray, file14_sh1_interp, der=0,ext=3) #IR filter sheet 34

#Harmonics 501 and coating thickness 350 Incident angle 0
file15_sh1_interp = interpolate.splrep(file15_sh1[:,0],file15_sh1[:,1],s=0)
H501_0_h350 = interpolate.splev(XArray, file15_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 0 Incident angle 1.0123
file16_sh1_interp = interpolate.splrep(file16_sh1[:,0],file16_sh1[:,1],s=0)
H301_h0 = interpolate.splev(XArray, file16_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 1 Incident angle 1.0123
file17_sh1_interp = interpolate.splrep(file17_sh1[:,0],file17_sh1[:,1],s=0)
H301_h1 = interpolate.splev(XArray, file17_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 20 Incident angle 1.0123
file18_sh1_interp = interpolate.splrep(file18_sh1[:,0],file18_sh1[:,1],s=0)
H301_h20 = interpolate.splev(XArray, file18_sh1_interp, der=0,ext=3)

#Harmonics 301 and coating thickness 250 Incident angle 1.0123
file19_sh1_interp = interpolate.splrep(file19_sh1[:,0],file19_sh1[:,1],s=0)
H301_h250 = interpolate.splev(XArray, file19_sh1_interp, der=0,ext=3) #IR filter sheet 34

#Harmonics 301 and coating thickness 350 Incident angle 1.0123
file20_sh1_interp = interpolate.splrep(file20_sh1[:,0],file20_sh1[:,1],s=0)
H301_h350 = interpolate.splev(XArray, file20_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 0 Incident angle 1.0123
file21_sh1_interp = interpolate.splrep(file21_sh1[:,0],file21_sh1[:,1],s=0)
H401_h0 = interpolate.splev(XArray, file21_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 1 Incident angle 1.0123
file22_sh1_interp = interpolate.splrep(file22_sh1[:,0],file22_sh1[:,1],s=0)
H401_h1 = interpolate.splev(XArray, file22_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 20 Incident angle 1.0123
file23_sh1_interp = interpolate.splrep(file23_sh1[:,0],file23_sh1[:,1],s=0)
H401_h20 = interpolate.splev(XArray, file23_sh1_interp, der=0,ext=3)

#Harmonics 401 and coating thickness 250 Incident angle 1.0123
file24_sh1_interp = interpolate.splrep(file24_sh1[:,0],file24_sh1[:,1],s=0)
H401_h250 = interpolate.splev(XArray, file24_sh1_interp, der=0,ext=3) #IR filter sheet 34

#Harmonics 401 and coating thickness 350 Incident angle 1.0123
file25_sh1_interp = interpolate.splrep(file25_sh1[:,0],file25_sh1[:,1],s=0)
H401_h350 = interpolate.splev(XArray, file25_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 0 Incident angle 1.0123
file26_sh1_interp = interpolate.splrep(file26_sh1[:,0],file26_sh1[:,1],s=0)
H501_h0 = interpolate.splev(XArray, file26_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 1 Incident angle 1.0123
file27_sh1_interp = interpolate.splrep(file27_sh1[:,0],file27_sh1[:,1],s=0)
H501_h1 = interpolate.splev(XArray, file27_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 20 Incident angle 1.0123
file28_sh1_interp = interpolate.splrep(file28_sh1[:,0],file28_sh1[:,1],s=0)
H501_h20 = interpolate.splev(XArray, file28_sh1_interp, der=0,ext=3)

#Harmonics 501 and coating thickness 250 Incident angle 1.0123
file29_sh1_interp = interpolate.splrep(file29_sh1[:,0],file29_sh1[:,1],s=0)
H501_h250 = interpolate.splev(XArray, file29_sh1_interp, der=0,ext=3) #IR filter sheet 34

#Harmonics 501 and coating thickness 350 Incident angle 1.0123
file30_sh1_interp = interpolate.splrep(file30_sh1[:,0],file30_sh1[:,1],s=0)
H501_h350 = interpolate.splev(XArray, file30_sh1_interp, der=0,ext=3)

workbook = xlsxwriter.Workbook('H501_h350.xlsx')
worksheet = workbook.add_worksheet()
H501_h350_array = [H501_h350]
row = 0
for col, data in enumerate(H501_h350_array):
    worksheet.write_column(row, col, data)
workbook.close()


# Compute the area under' the EM curves using the composite Simpson's rule in different number of Harmonics and different thickness of coating.

##############  Area under the H301_0_h0  ##################
area_H301_0_h0 = simps(H301_0_h0, dx=dX)
print("area_H301_0_h0 =", '{:.2f}'.format(area_H301_0_h0))

##############  Area under the H301_0_h1  ##################
area_H301_0_h1 = simps(H301_0_h1, dx=dX)
print("area_H301_0_h1 =", '{:.2f}'.format(area_H301_0_h1))

##############  Area under the H301_0_h20  ##################
area_H301_0_h20 = simps(H301_0_h20, dx=dX)
print("area_H301_0_h20 =", '{:.2f}'.format(area_H301_0_h20))

##############  Area under the H301_0_h250  ##################
area_H301_0_h250 = simps(H301_0_h250, dx=dX)
print("area_H301_0_h250 =", '{:.2f}'.format(area_H301_0_h250))

##############  Area under the H301_0_h350  ##################
area_H301_0_h350 = simps(H301_0_h350, dx=dX)
print("area_H301_0_h350 =", '{:.2f}'.format(area_H301_0_h350))

##############  Area under the H401_0_h0  ##################
area_H401_0_h0 = simps(H301_0_h0, dx=dX)
print("area_H401_0_h0 =", '{:.2f}'.format(area_H401_0_h0))

##############  Area under the H401_0_h1  ##################
area_H401_0_h1 = simps(H401_0_h1, dx=dX)
print("area_H401_0_h1 =", '{:.2f}'.format(area_H401_0_h1))

##############  Area under the H401_0_h20  ##################
area_H401_0_h20 = simps(H401_0_h20, dx=dX)
print("area_H401_0_h20 =", '{:.2f}'.format(area_H401_0_h20))

##############  Area under the H401_0_h250  ##################
area_H401_0_h250 = simps(H401_0_h250, dx=dX)
print("area_H401_0_h250 =", '{:.2f}'.format(area_H401_0_h250))

##############  Area under the H401_0_h350  ##################
area_H401_0_h350 = simps(H401_0_h350, dx=dX)
print("area_H401_0_h350 =", '{:.2f}'.format(area_H401_0_h350))

##############  Area under the H501_0_h0  ##################
area_H501_0_h0 = simps(H301_0_h0, dx=dX)
print("area_H501_0_h0 =", '{:.2f}'.format(area_H501_0_h0))

##############  Area under the H501_0_h1  ##################
area_H501_0_h1 = simps(H501_0_h1, dx=dX)
print("area_H501_0_h1 =", '{:.2f}'.format(area_H501_0_h1))

##############  Area under the H501_0_h20  ##################
area_H501_0_h20 = simps(H501_0_h20, dx=dX)
print("area_H501_0_h20 =", '{:.2f}'.format(area_H501_0_h20))

##############  Area under the H501_0_h250  ##################
area_H501_0_h250 = simps(H501_0_h250, dx=dX)
print("area_H501_0_h250 =", '{:.2f}'.format(area_H501_0_h250))

##############  Area under the H501_0_h350  ##################
area_H501_0_h350 = simps(H501_0_h350, dx=dX)
print("area_H501_0_h350 =", '{:.2f}'.format(area_H501_0_h350))


##############  Area under the H301_h0  ##################
area_H301_h0 = simps(H301_0_h0, dx=dX)
print("area_H301_h0 =", '{:.2f}'.format(area_H301_h0))

##############  Area under the H301_h1  ##################
area_H301_h1 = simps(H301_h1, dx=dX)
print("area_H301_h1 =", '{:.2f}'.format(area_H301_h1))

##############  Area under the H301_h20  ##################
area_H301_h20 = simps(H301_h20, dx=dX)
print("area_H301_h20 =", '{:.2f}'.format(area_H301_h20))

##############  Area under the H301_h250  ##################
area_H301_h250 = simps(H301_h250, dx=dX)
print("area_H301_h250 =", '{:.2f}'.format(area_H301_h250))

##############  Area under the H301_h350  ##################
area_H301_h350 = simps(H301_h350, dx=dX)
print("area_H301_h350 =", '{:.2f}'.format(area_H301_h350))

##############  Area under the H401_h0  ##################
area_H401_h0 = simps(H301_h0, dx=dX)
print("area_H401_h0 =", '{:.2f}'.format(area_H401_h0))

##############  Area under the H401_h1  ##################
area_H401_h1 = simps(H401_h1, dx=dX)
print("area_H401_h1 =", '{:.2f}'.format(area_H401_h1))

##############  Area under the H401_h20  ##################
area_H401_h20 = simps(H401_h20, dx=dX)
print("area_H401_h20 =", '{:.2f}'.format(area_H401_h20))

##############  Area under the H401_h250  ##################
area_H401_h250 = simps(H401_h250, dx=dX)
print("area_H401_h250 =", '{:.2f}'.format(area_H401_h250))

##############  Area under the H401_h350  ##################
area_H401_h350 = simps(H401_h350, dx=dX)
print("area_H401_h350 =", '{:.2f}'.format(area_H401_h350))

##############  Area under the H501_h0  ##################
area_H501_h0 = simps(H301_h0, dx=dX)
print("area_H501_h0 =", '{:.2f}'.format(area_H501_h0))

##############  Area under the H501_h1  ##################
area_H501_h1 = simps(H501_h1, dx=dX)
print("area_H501_h1 =", '{:.2f}'.format(area_H501_h1))

##############  Area under the H501_h20  ##################
area_H501_h20 = simps(H501_h20, dx=dX)
print("area_H501_h20 =", '{:.2f}'.format(area_H501_h20))

##############  Area under the H501_h250  ##################
area_H501_h250 = simps(H501_h250, dx=dX)
print("area_H501_h250 =", '{:.2f}'.format(area_H501_h250))

##############  Area under the H501_h350  ##################
area_H501_h350 = simps(H501_0_h350, dx=dX)
print("area_H501_h350 =", '{:.2f}'.format(area_H501_h350))

################### plot interpolated data ###################

plt.figure(1)
plt.plot(XArray,H301_0_h0, color='m', label='sun_Irradiance')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Irradiance (\u03bc W/cm²/nm)')
plt.title('Sun_Irradiance')
plt.legend(loc='best')
plt.savefig('sun_I.png', dpi=300)
plt.show()