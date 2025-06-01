#!/usr/bin/env python3
#   
#    eval_vasp_ml: Evaluate the content of a VASP on the fly
#      learning training set file ML_AB and generate a number 
#      of plots to enable a detailed analysis of the learning 
#    Part of VASP4CLINT
#     Andreas Mölkner, 2025 (andreas.moelkner@fau.de)
#

#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import math

###################################################################################################
#
#Read in the Informations of the Header of the ML_AB file.
#Only ML_AB with 6 or less Elements can be analyzed at the moment.
#
print('''
This script reads in a VASP on the fly learning training set 
file (ML_AB) and evaluates the progress of the learning from 
its contents. A number of plot files is written based on the 
Python matplotlib library.
The analysis contains: Basis functions per structure/atom
CTIFOR (Bayesian error threshold), Etot, Force magnitude ...
Furthermore, a new ML_AB file is generated in which the basis 
functions of high energy structures are removed.
''')

# Size of the font for plotting
plt.rcParams.update({'font.size': 12})

print("Read in the ML_AB file...")

ML_AB = open('ML_AB','r')
ML_AB_list = ML_AB.readlines()
ML_AB.close()

num_conf=int(ML_AB_list[4].split()[0])
num_elems=int(ML_AB_list[8].split()[0])

if num_elems <= 3:
    elems_list = ML_AB_list[12].split()
    max_atoms = int(ML_AB_list[16].split()[0])
    max_atoms_per = int(ML_AB_list[20].split()[0])
    num_bsets_list =  ML_AB_list[32].split()
    end_header=35
elif num_elems > 3  and num_elems <= 6 :
    elems_list = ML_AB_list[12].split()
    elems_list2 = ML_AB_list[13].split()
    num_bsets_list = ML_AB_list[35].split()
    num_bsets_list2 = ML_AB_list[36].split()

    for x in range(len(elems_list2)):
        elems_list.append(elems_list2[x])
        num_bsets_list.append(num_bsets_list2[x])

    max_atoms = int(ML_AB_list[17].split()[0])
    max_atoms_per = int(ML_AB_list[21].split()[0])

    end_header=39

else:
    print("Too many Elements")




print("Number of configurations")
print(num_conf)
print("Number of Elements")
print(num_elems)
print("Elements")
print(elems_list)
print("Max. number of atoms per structure")
print(max_atoms)
print("Max. number of atoms per element")
print(max_atoms_per)
print("Number of basis sets per element")
print(num_bsets_list)
print(ML_AB_list[end_header].strip())
print("")
###################################################################################################
#
#Now the Basis sets for the respective elements are read and analyzed!
#At first lists for the respective Elements are generated and stored in a dictionary.

basis_set_dic = {}
for x in range(len(elems_list)): 
    basis_set_dic[elems_list[x]+"_basis_list"] = []
    basis_set_dic[elems_list[x]+"_atoms_list"] = []
#print(basis_set_dic)


#
#now the data is read for each element
#

start=end_header+1
for x in range(len(elems_list)):
    for y in range(start,start+int(num_bsets_list[x]),1):
        basis_set_dic[elems_list[x]+"_basis_list"].append(int(ML_AB_list[y].split()[0]))
        basis_set_dic[elems_list[x]+"_atoms_list"].append(int(ML_AB_list[y].split()[1]))
    start=start+int(num_bsets_list[x])+3

#
#Now different histograms for the data of each element are gerated
#to see how often which Basis set is used and how many basis sets each atom has.
#

for x in range(len(elems_list)):
    fig,ax= plt.subplots()
    ax.hist(basis_set_dic[elems_list[x]+"_basis_list"],bins=num_conf)
    ax.set_title(elems_list[x]+' basis set histogram')
    ax.set_xlabel("Configuration Number")
    ax.set_ylabel("count [a.u.]")
    ax.set_xbound(lower=0,upper=num_conf)
    plt.savefig(elems_list[x]+'_basis_set.jpg',dpi=300)

    fig,ax= plt.subplots()
    ax.hist(basis_set_dic[elems_list[x]+"_atoms_list"],bins=max_atoms_per)
    ax.set_xlabel("Atom Number")
    ax.set_ylabel("count [a.u.]")
    #ax.set_xbound(lower=0)
    ax.set_title(elems_list[x]+' Atoms histogram')
    plt.savefig(elems_list[x]+'_Atoms.jpg',dpi=300)

###################################################################################################
#
# Extract and plot CTIFOR and energy values

ctifor_list=[]
for x in range(len(ML_AB_list)):
    if "CTIFOR" in ML_AB_list[x]:
        ctifor_list.append(float(ML_AB_list[x+2].strip()))

toten_list=[]
for x in range(len(ML_AB_list)):
    if "Total energy (eV)" in ML_AB_list[x]:
        toten_list.append(float(ML_AB_list[x+2].strip()))


#plot the results!
fig,ax=plt.subplots()
#plt.rcParams['font.size'] = '18'
ax.plot(ctifor_list)
ax.set_xlabel("Number of Configurations")
ax.set_ylabel("CTIFOR")
ax.set_xbound(lower=0,upper=num_conf)
ax.set_title("Progression of CTIFOR")
plt.savefig("CTIFOR.jpg",dpi=300)
plt.cla()


ax.plot(toten_list)
ax.set_xlabel("Number of Configurations")
ax.set_ylabel("Total energy [eV]")
ax.set_xbound(lower=0,upper=num_conf)
ax.set_title("Progression of Etot")
plt.savefig("Etot.jpg",dpi=300)
plt.cla()

nl="\n"
nt="\t"


print("test")
out=open("E_ctifor.dat","w")
out.write("# Etot      CTIFOR")
for x in range(len(ctifor_list)):
    out.write(str("{:9.7f}".format(toten_list[x])) + nt + str("{:9.7f}".format(ctifor_list[x])) + nl)
out.close()

###################################################################################################
#
# Extract z cooridnates of the atom of the configurations
#

zcoor_dic = {}
for x in range(1,max_atoms+1,1):
    zcoor_dic["atom"+str(x)]=[]


for x in range(len(ML_AB_list)):
    if "Atomic positions (ang.)" in ML_AB_list[x]:
        i=0
        for y in range(x+2,x+2+max_atoms,1):
            i=i+1
            zcoor_dic["atom"+str(i)].append(float(ML_AB_list[y].split()[2].strip()))

mean_zcoor = []
for x in range(1,max_atoms+1,1):
    mean = sum(zcoor_dic["atom"+str(x)])/len(zcoor_dic["atom"+str(x)])
    mean_zcoor.append(mean)


#
# Extract Froces and calculate the magnitude
#
force_dic = {}
for x in range(1,max_atoms+1,1):
    force_dic["force"+str(x)]=[]

mean_force_structure = []
for x in range(len(ML_AB_list)):
    if "Forces (eV ang.^-1)" in ML_AB_list[x]:
        i=0
        mag_sum = 0
        for y in range(x+2,x+2+max_atoms,1):
            i=i+1
            dummy=ML_AB_list[y].split()
            dummy[0]=float(dummy[0])
            dummy[1]=float(dummy[1])
            dummy[2]=float(dummy[2])
            magnitude=math.sqrt(dummy[0]**2 + dummy[1]**2 + dummy[2]**2)
            force_dic["force"+str(i)].append(magnitude)
            mag_sum  = mag_sum + magnitude
        mean_force_structure.append(mag_sum/float(max_atoms))
#
#Analyze the mean force and Etot of each structure for outliers
#


toten_array = np.array(toten_list)
q1 = np.percentile(toten_array, 25)
q3 = np.percentile(toten_array, 75)
iqr = q3 - q1
lower_fence = q1 - 1.5 * iqr
upper_fence = q3 + 1.5 * iqr
E_outliers = np.where((toten_array < lower_fence) | (toten_array > upper_fence))
E_outliers = np.array(E_outliers)
dummy = E_outliers.tolist()
E_outliers = dummy[0]
E_outliers = [ x+1 for x in E_outliers]
#print(E_outliers)

force_array = np.array(mean_force_structure)
q1 = np.percentile(force_array, 25)
q3 = np.percentile(force_array, 75)
iqr = q3 - q1
lower_fence = q1 - 1.5 * iqr
upper_fence = q3 + 1.5 * iqr
F_outliers = np.where((force_array < lower_fence) | (force_array > upper_fence))
F_outliers = np.array(F_outliers)
dummy = F_outliers.tolist()
F_outliers = dummy[0]
F_outliers = [ x+1 for x in F_outliers]
#print(F_outliers)

E_F_outliers = E_outliers + F_outliers
E_F_outliers = list(set(E_F_outliers))
E_F_outliers.sort()

print("The following structures have a very High energy or force magnitude:")
print(E_F_outliers)
print("Number of structures")
print(len(E_F_outliers))
#
#Calculate the mean magnitude of forces for each structure and each atom and plot it:
#
mean_force_atom = []
for x in range(1,max_atoms+1,1):
    mean = sum(force_dic["force"+str(x)])/len(force_dic["force"+str(x)])
    mean_force_atom.append(mean)

ax.plot(mean_force_structure)
ax.set_xlabel("Number of Configurations")
ax.set_ylabel("Mean Force Magnitude [eV ang.^-1]")
ax.set_xbound(lower=0,upper=num_conf)
ax.set_title("Progression of Mean Force Magnitude")
plt.savefig("Force_structure.jpg",dpi=300)
plt.cla()

ax.plot(mean_force_atom)
ax.set_xlabel("Atom Number")
ax.set_ylabel("Mean Force Magnitude [eV ang.^-1]")
ax.set_xbound(lower=0,upper=max_atoms)
ax.set_title("Mean Force Magnitude of each Atom")
plt.savefig("Force_atoms.jpg",dpi=300)
plt.cla()

ax.plot(mean_zcoor,mean_force_atom,marker='.',linestyle='None')
ax.set_xlabel("Mean Z-Coordinate [Ang.]")
ax.set_ylabel("Mean Force Magnitude [eV ang.^-1]")
ax.set_title("Mean Force Magnitude vs Mean Zcoord of each Atom")
plt.savefig("Force_Z.jpg",dpi=300)
plt.cla()



#
#Find the highest force magnitude of each atom and store it in list.
#Get the N atoms with highest forces and print the forces and zcoordiantes for them.
#
max_force_list = []
for x in range(1,max_atoms+1,1):
    max_force_list.append(max(force_dic["force"+str(x)]))

N=10
res =sorted(range(len(max_force_list)), key = lambda sub: max_force_list[sub])[-N:]
res = [x + 1 for x in res]
res = [str(x) for x in res]

print("The Atoms with the 10 highest froce magnitudes are:")
print(res)


#
#Print zcoordinates and magnitudes of given Atoms:
#Or Print the stuff for the atoms with the Highest force magnitudes
#

a_list=res
#a_list = input ("Give Number of the to analyze Atoms (e.g., 1,4,70): ").split(",")


max_force_print_list= []
for x in range(len(a_list)):
    ax.plot(zcoor_dic["atom"+a_list[x]])
    ax.set_xlabel("Number of Configurations")
    ax.set_ylabel("Z-Coordinate [Ang.]")
    ax.set_xbound(lower=0,upper=num_conf)
    ax.set_title("Progression of the Z-Coordinate of Atom"+a_list[x])
    plt.savefig("atom"+a_list[x]+"_zcoord.jpg",dpi=300)
    plt.cla()

    ax.plot(force_dic["force"+a_list[x]])
    ax.set_xlabel("Number of Configurations")
    ax.set_ylabel("Magnitude of Force [eV ang.^-1]")
    ax.set_xbound(lower=0,upper=num_conf)
    ax.set_title("Progression of the Force Magnitude of Atom"+a_list[x])
    plt.savefig("atom"+a_list[x]+"_force.jpg",dpi=300)
    plt.cla()

    eval_list = []
    eval_list.append(a_list[x])
    index=force_dic["force"+a_list[x]].index(max(force_dic["force"+a_list[x]]))
    eval_list.append(index)
    eval_list.append(zcoor_dic["atom"+a_list[x]][index])
    eval_list.append(force_dic["force"+a_list[x]][index])

    max_force_print_list.append(eval_list)
print("")
print("Data of this atoms: Atomnumber, Structure, Z-coordiante, Force magnitude")
print(max_force_print_list)


###################################################################################################
#
#Print the structure of outliers into XDATCAR_outliers file
#
#Start with the Header:
header_atoms = []
header_vectors = []
for x in range(len(ML_AB_list)):
    if "Atom types and atom numbers" in ML_AB_list[x]:
        for y in range(x+2,x+max_atoms,1):
            if "==================================================" in ML_AB_list[y]:
                break
            else:
                header_atoms.append(ML_AB_list[y].split())
    if "Primitive lattice vectors (ang.)" in ML_AB_list[x]:
        for y in range(x+2,x+max_atoms,1):
            if "==================================================" in ML_AB_list[y]:
                break
            else:
                header_vectors.append(ML_AB_list[y].split())

        break


output = open("OUTLIERCAR",'w')
output.write("Structures with High Energy and Force magnitude!\n")
output.write("1\n")
for x in range(len(header_vectors)):
    output.write(header_vectors[x][0]+"\t"+header_vectors[x][1]+"\t"+header_vectors[x][2]+"\n")
for x in range(len(header_atoms)):
    output.write(header_atoms[x][0]+"\t")
output.write("\n")
for x in range(len(header_atoms)):
    output.write(header_atoms[x][1]+"\t")
output.write("\n")


x=0
for i in E_F_outliers:
    for x in range(x,len(ML_AB_list),1):
        if 'Configuration num.' in ML_AB_list[x]:
            number=ML_AB_list[x].strip()[-6:].strip()
            if number == str(i):
                #system_coord_list = []
                for y in range(x,len(ML_AB_list),1):
                    if "Atomic positions (ang.)" in ML_AB_list[y]:
                        output.write("Cartesian configuration= "+str(i)+"\n")
                        for z in  range(y+2,y+2+max_atoms,1):
                            #system_coord_list.append(ML_AB_list[z])
                            output.write(ML_AB_list[z])
                        #print(ML_AB_list[x],system_coord_list[0],number,i)
                        break

                break

output.close()

###################################################################################################
#
#Generate ML_ABN file with the OUtliers removed.
#
#Adept the list of basis functions from the beginning:
#

print("")
print("Two ML_AB files will be generated:")
print("ML_AB_vasp: High energy and force strucutres are removed from the basis set but the structures")
print("            are kept in the file. This can directly be used by VASP")
print("ML_AB_cutted: Additionally to the basis sets the structures are also removed")

removed_num = []
for x in range(len(elems_list)):
    removed_list = []
    i=0
    for y in range(len(basis_set_dic[elems_list[x]+"_basis_list"])):
        if basis_set_dic[elems_list[x]+"_basis_list"][y] in E_F_outliers:
            i=i+1
            removed_list.append(y)
    removed_list.reverse()
    for y in removed_list:        
        basis_set_dic[elems_list[x]+"_basis_list"].pop(y)
        basis_set_dic[elems_list[x]+"_atoms_list"].pop(y)
    removed_num.append(len(removed_list))
    #print(len(removed_list),i)


#
#Start writing the file ML_AB_vasp
#

ML_ABN= open("ML_AB_vasp","w")

seperator1="**************************************************\n"
seperator2="--------------------------------------------------\n"


string=""
if num_elems <= 3:
    for x in range(0,end_header-3,1):
        ML_ABN.write(ML_AB_list[x])
    for x in range(len(elems_list)):
        string=string+"      "+str(int(num_bsets_list[x])-removed_num[x])
    ML_ABN.write(string+"\n")

elif num_elems > 3 and num_elems <= 6:
    for x in range(0,end_header-4,1):
        ML_ABN.write(ML_AB_list[x])
    ML_ABN.write("      "+str(int(num_bsets_list[0])-removed_num[0])+"      "+str(int(num_bsets_list[1])-removed_num[1])+"      "+str(int(num_bsets_list[2])-removed_num[2])+"\n")
    for x in range(3,len(elems_list),1):
        string=string+"      "+str(int(num_bsets_list[x])-removed_num[x])
    ML_ABN.write(string+"\n")

for x in range(len(elems_list)):
    ML_ABN.write(seperator1)
    ML_ABN.write("     Basis set for "+elems_list[x]+"\n")
    ML_ABN.write(seperator2)
    for y in range(len(basis_set_dic[elems_list[x]+"_basis_list"])):
        ML_ABN.write("   "+str(basis_set_dic[elems_list[x]+"_basis_list"][y])+"   "+str(basis_set_dic[elems_list[x]+"_atoms_list"][y])+"\n")

ML_ABN.write(seperator1)

#
#Determine length of structure block:
#

length_list = []
for x in range(len(ML_AB_list)):
    if 'Configuration num.' in ML_AB_list[x]:
        length_list.append(x)
    if len(length_list) == 2:
        break

length = length_list[1]-length_list[0]

for x in range(length_list[0],len(ML_AB_list),1):
    ML_ABN.write(ML_AB_list[x])

ML_ABN.close()

###################################################################################################
#Printing ML_AB_cutted

ML_ABN= open("ML_AB_cutted","w")

for x in range(0,4,1):
    ML_ABN.write(ML_AB_list[x])

ML_ABN.write("       "+str(num_conf-len(E_F_outliers))+"\n")


string=""
if num_elems <= 3:
    for x in range(5,end_header-3,1):
        ML_ABN.write(ML_AB_list[x])
    for x in range(len(elems_list)):
        string=string+"      "+str(int(num_bsets_list[x])-removed_num[x])
    ML_ABN.write(string+"\n")
    
elif num_elems > 3 and num_elems <= 6:
    for x in range(5,end_header-4,1):
        ML_ABN.write(ML_AB_list[x])
    ML_ABN.write("      "+str(int(num_bsets_list[0])-removed_num[0])+"      "+str(int(num_bsets_list[1])-removed_num[1])+"      "+str(int(num_bsets_list[2])-removed_num[2])+"\n")
    for x in range(3,len(elems_list),1):
        string=string+"      "+str(int(num_bsets_list[x])-removed_num[x])
    ML_ABN.write(string+"\n")

for x in range(len(elems_list)):
    ML_ABN.write(seperator1)
    ML_ABN.write("     Basis set for "+elems_list[x]+"\n")
    ML_ABN.write(seperator2)
    for y in range(len(basis_set_dic[elems_list[x]+"_basis_list"])):
        ML_ABN.write("   "+str(basis_set_dic[elems_list[x]+"_basis_list"][y])+"   "+str(basis_set_dic[elems_list[x]+"_atoms_list"][y])+"\n")

ML_ABN.write(seperator1)

x = length_list[0]
while x < len(ML_AB_list):
    if 'Configuration num.' in ML_AB_list[x]:
        number=ML_AB_list[x].strip()[-6:].strip()
        if int(number) in E_F_outliers:
            x = x +length
        else:
            ML_ABN.write(ML_AB_list[x])
            x=x+1
    else:
        ML_ABN.write(ML_AB_list[x])
        x=x+1

ML_ABN.close()
