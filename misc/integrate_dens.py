#!/usr/bin/env python3
import sys
import os
import numpy as np
import random as ran
from scipy.signal import argrelmin
import matplotlib.pyplot as plt

    



filename=str(sys.argv[1])

original = open(filename, 'r')
original_list = original.readlines()
original.close()
data=original_list

Ncol = len(data[3].split())

#It is assumed that the first column contains the z coordinates and the last column the total density
#the other densities are in the dens_elems directory
dens_elems= {}
for x in range(1,Ncol-1,1):
    dens_elems["El"+str(x)] = []

#Get data from the input file into the different lists
zcor=[]
totdens=[]
for x in range(1,len(data),1):
    zcor.append(float(data[x].split()[0]))
    totdens.append(float(data[x].split()[-1]))
    for y in range(1,Ncol-1,1):
        dens_elems["El"+str(y)].append(data[x].split()[y])


#find all minima in the total denisty
totdens=np.array(totdens,dtype=float)
minima_all=argrelmin(totdens,order=10)[0]


#filter the minima 
minima = []
for x in range(len(minima_all)):
    if totdens[minima_all[x]] > 0.01:
        minima.append(minima_all[x])


out=open("concs.dat","w")

nl="\n"
nt="\t"

out.write("Elemental concentrations calculated by integrate.py!" + nl)
out.write("DISCLAMER: Some concentrations may be unphysical or wrong depending on the size of the cell!"+nl)
out.write(nl)

out.write("Area/Layer"+nt+nt)
for x in range(1,Ncol-1,1):
    out.write("Element_"+str(x)+nt+nt) 
out.write(nl)


#now the different denisties shall be calculted

#first layer for each element!
def integrate(index,index2=None,mode="s"):
    work_dic = {}
    for x in range(1,Ncol-1,1):
        work_dic["El"+str(x)] = []


    totdens_work = []
    x_list = []
    for x in range(len(totdens)):
        if mode == "s":
            if x <= minima[index]:    
                totdens_work.append(totdens[x])
                x_list.append(x)
                for  y in range(1,Ncol-1,1):
                      work_dic["El"+str(y)].append(dens_elems["El"+str(y)][x])
        elif mode == "g":
            if x >= minima[index]:
                totdens_work.append(totdens[x])
                x_list.append(x)
                for  y in range(1,Ncol-1,1):
                      work_dic["El"+str(y)].append(dens_elems["El"+str(y)][x])
        elif mode == "c":
            if x >= minima[index] and x <= minima[index2]:
                totdens_work.append(totdens[x])
                x_list.append(x)
                for  y in range(1,Ncol-1,1):
                      work_dic["El"+str(y)].append(dens_elems["El"+str(y)][x])

#convert work lists to np.arrays
    totdens_work=np.array(totdens_work,dtype=float)
    x_list=np.array(x_list,dtype=int)
    for x in range(1,Ncol-1,1):
        work_dic["El"+str(x)] = np.array(work_dic["El"+str(x)],dtype=float)

    results  = []
    Itot = np.trapz(totdens_work)
    for x in range(1,Ncol-1,1):
        dummy = np.trapz(work_dic["El"+str(x)])
        result = (dummy/Itot)*100
        results.append(result)

    return results

def integrate2(index,length,mode="s"):
    work_dic = {}
    for x in range(1,Ncol-1,1):
        work_dic["El"+str(x)] = []


    totdens_work = []
    x_list = []
    for x in range(len(totdens)):
        if mode == "s":
            if zcor[x] <= zcor[minima[index]] + length:
                totdens_work.append(totdens[x])
                x_list.append(x)
                for  y in range(1,Ncol-1,1):
                    work_dic["El"+str(y)].append(dens_elems["El"+str(y)][x])
        elif mode == "g":
            if zcor[x] >= zcor[minima[index]] - length:
                totdens_work.append(totdens[x])
                x_list.append(x)
                for  y in range(1,Ncol-1,1):
                      work_dic["El"+str(y)].append(dens_elems["El"+str(y)][x])

#convert work lists to np.arrays
    totdens_work=np.array(totdens_work,dtype=float)
    x_list=np.array(x_list,dtype=int)
    for x in range(1,Ncol-1,1):
        work_dic["El"+str(x)] = np.array(work_dic["El"+str(x)],dtype=float)

    results  = []
    Itot = np.trapz(totdens_work)
    for x in range(1,Ncol-1,1):
        dummy = np.trapz(work_dic["El"+str(x)])
        result = (dummy/Itot)*100
        results.append(result)

    return results

def calc_final(array1,array2,area):
    results = []
    for x in range(len(array1)):
        results.append((array1[x]+array2[x])/2)

    out.write(area+nt+nt+nt)
    for x in range(len(results)):
        out.write(str("{:9.7f}".format(results[x]))+nt+nt)
    out.write(nl) 


calc_final(integrate(0),integrate(-1,mode="g"),area="1")
calc_final(integrate(0,1,"c"),integrate(-2,-1,"c"),area="2")
calc_final(integrate(0,2,"c"),integrate(-3,-1,"c"),area="3")
calc_final(integrate(1),integrate(-2,mode="g"),area="1+2")
calc_final(integrate(2),integrate(-3,mode="g"),area="1+2+3")

out.write(nl)


calc_final(integrate2(0,5),integrate2(-1,5,mode ="g"),area="1+5A")
calc_final(integrate2(0,10),integrate2(-1,10,mode ="g"),area="1+10A")
calc_final(integrate2(0,20),integrate2(-1,20,mode ="g"),area="1+20A")

print("first Layers")
print(integrate(0),zcor[minima[0]])
print(integrate(-1,mode="g"),zcor[minima[-1]])
print()
print("Second layers")
print(integrate(1),zcor[minima[1]])
print(integrate(-2,mode="g"),zcor[minima[-2]])
print()
print("only second layer")
print(integrate(0,1,"c"))
print(integrate(-2,-1,"c"))
print()
print("first layer plus 5")
print(integrate2(0,5),zcor[minima[0]] + 5)
print(integrate2(-1,5,mode ="g"),zcor[minima[-1]] - 5)
print()
print("first layer plus 10")
print(integrate2(0,10),zcor[minima[0]] + 10)
print(integrate2(-1,10,mode ="g"),zcor[minima[-1]] - 10)
print()
print("first layer plus 20")
print(integrate2(0,20),zcor[minima[0]] + 20)
print(integrate2(-1,20,mode ="g"),zcor[minima[-1]] - 20)







