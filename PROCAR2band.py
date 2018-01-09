#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import os
import time

ylims=(-10,10)

# the common case
#import subprocess
#B=subprocess.getoutput('grep band PROCAR')

def get_case():
    '''
    Find What the calculation is.
    '''
    ispin=os.popen('grep ISPIN OUTCAR').read().split()[2]
    soc=os.popen('grep LSORBIT OUTCAR').read().split()[2]
    if soc == 'T':
        case = 1
        nspin = 1
    else:
        if ispin=='2':
            case = 2
            nspin = 2
        else:
            case = 0
            nspin = 1
    return case,nspin

def get_fermi(filename='DOSCAR'):
    with open(filename,'r') as file:
        fermi=file.readlines()[5].split()[3]
    return float(fermi)

def bandinfo(filename='PROCAR'):
    with open(filename,'r') as file:
        bandinfo=file.readlines()[1].split()
    kptsnum=int(bandinfo[3])
    bandsnum=int(bandinfo[7])
    ionsnum=int(bandinfo[-1])
    return kptsnum,bandsnum,ionsnum

def get_klabels(filename='KPOINTS'):
    '''
    Find the klabels and the number of k-points per line.
    '''
    file=open(filename,'r')
    i=0
    kls=[]
    for line in file.readlines():
        if i == 1:
            knump=int(line.split()[0])
        elif i >= 4:
            key=line.split()
            if (key != []):
                kl=key[4]
                kls.append(kl)
        i=i+1
    num=int(len(kls[1:-1])/2)
    klabels=[kls[0]]
    for i in range(num):
        if kls[2*i+1]==kls[2*i+2]:
            klabels.append(kls[2*i+1])
        else:
            klabels.append(kls[2*i+1]+'|'+kls[2*i+2])
            #print(kls[2*i+1],kls[2*i+2],klabels)
    klabels.append(kls[-1])
    return knump,klabels

(case,nspin)=get_case()
(kptsnum,bandsnum,ionsnum)=bandinfo()

def get_kvec(filename='PROCAR'):
    procar=open(filename,'r')
    kvec=[]
    for line in procar.readlines():
        key=line.split()
        #print(key)
        if (key != []) and (key[0]=='k-point'):
            kvec.append([float(key[3]),float(key[4]),float(key[5])])
    procar.close()
    kvec=np.array(kvec[:kptsnum])
    return kvec

def get_kdist():
    kvec=get_kvec()
    dkvec=np.zeros((len(kvec)-1,3))
    for kv in range(len(kvec)-1):
        dkvec[kv]=kvec[kv+1]-kvec[kv]
    dk=np.zeros(len(dkvec))
    for kv in range(len(dkvec)):
        dk[kv]=np.sqrt(dkvec[kv,0]**2+dkvec[kv,1]**2+dkvec[kv,2]**2)
    kdist=[0,]
    value=0
    for i in range(len(dk)):
        value=value+dk[i]
        kdist.append(value)
        
    return kdist

def get_eigs(filename='PROCAR'):
    fermi=get_fermi()
    procar=open(filename,'r')
    eigs=[]
    for line in procar.readlines():
        key=line.split()
        #print(key)
        if (key != []) and (key[0]=='band'):
            eigs.append(float(key[4]))
    procar.close()
    eigs=np.array(eigs)-fermi
    
    def get_gap(eigs):
        eigs_pos=[]
        eigs_neg=[]
        for i in range(len(eigs)):
            if eigs[i] >=0:
                eigs_pos.append(eigs[i])
            else:
                eigs_neg.append(eigs[i])
            
        bandgap=np.array(eigs_pos).min()-np.array(eigs_neg).max()
        return bandgap
        

    if nspin == 1:
        bandgap=get_gap(eigs)
        eigs=eigs.reshape((kptsnum,bandsnum))
        results=(eigs,bandgap)
        
    else:
        eigs_up=eigs[0:bandsnum*kptsnum]
        eigs_dn=eigs[bandsnum*kptsnum:]
        
        gap_up=get_gap(eigs_up)
        gap_dn=get_gap(eigs_dn)

        eigs_up=eigs_up.reshape((kptsnum,bandsnum))
        eigs_dn=eigs_dn.reshape((kptsnum,bandsnum))
        results=(eigs_up,eigs_dn,gap_up,gap_dn)
    return results

(knump,klabels)=get_klabels()
kdist=np.array(get_kdist())
kdists=[kdist[0]]
xlims=(kdist.min(),kdist.max())

for i in range(len(klabels)-1):
    value=kdist[(i+1)*knump-1]
    kdists.append(value)

eigs_res=get_eigs()

def band_dat():
    
    file=open('bandspy.dat','w+')

    for j in range(bandsnum):
        for i in range(kptsnum+1):
            if i < len(kdist):
                #file.write('{:2.8f}  {:4.8f}\n'.format(kdist[i],eigs1[i,j]))
                if nspin == 1:
                    file.write('%2.8f    %1.8f\n' %(kdist[i],eigs_res[0][i,j]))
                else:
                    file.write('%2.8f    %1.8f    %1.8f\n' %(kdist[i],eigs_res[0][i,j],eigs_res[1][i,j]))
                    
            else:
                file.write('\n')
    file.close()
            

def plotband():
    if nspin == 1:
        for i in range(bandsnum):
            plt.plot(kdist,eigs_res[0][:,i],c='blue')
            if eigs_res[1] < 1e-3:
                plt.title('BE CARE! Small Band Gap: %.6f eV' %eigs_res[1])
            else:
                plt.title('Band Gap: %.6f eV' %eigs_res[1])
    else:
        for i in range(bandsnum):
            plt.plot(kdist,eigs_res[0][:,i],c='blue')
            plt.plot(kdist,eigs_res[1][:,i],c='red')
            #plt.legend(loc=0)
            plt.title('gap_up(blue): %.6f, gap_dn(red): %.6f (eV)' %(eigs_res[2],eigs_res[3]))
        
    plt.xlim(xlims)
    plt.ylim(ylims)
    plt.hlines(0,xlims[0],xlims[1],colors='green',linestyles='--',linewidths=0.8)
    for i in range(len(kdists[1:-1])):
        plt.vlines(kdists[i+1],ylims[0],ylims[1],colors='black',linestyles='--',linewidths=0.8)

    plt.xticks(kdists,klabels)
    plt.savefig('bands.pdf')
    plt.show()

def report():
    localtime = time.asctime(time.localtime(time.time()))
    with open('bands.log','w+') as file:
        file.write('*********************************************************\n')
        file.write(localtime)
        file.write('\n')
        file.write('Xiangru Kong; Python3, numpy, matplotlibv2+;20180109\n')
        file.write('This version is only tested with PBE calculations by vasp5.3+\n')
        file.write('*********************************************************\n\n')
        
        if case == 0:
            file.write('The non-spin polorized calculation without SOC:\n\n')
        elif case == 1:
            file.write('The non-colliear calculation with SOC:\n\n')
        else:
            file.write('The spin polorized calculation without SOC:\n\n')
            
        file.write('number of ions: %d\n' %ionsnum)
        file.write('number of bands: %d\n' %bandsnum)
        file.write('number of k-points: %d\n' %kptsnum)
        file.write('number of k-points per high symmetry line: %d\n' %knump)
        file.write('number of high symmetry lines: %d\n' %(int(kptsnum/knump)))
        file.write('the high symmetry points: ' + str(klabels) +'\n')
        file.write('the position of high symmetry points: ' + str(kdists) +'\n\n')
        if nspin == 1:
            if eigs_res[1] < 1e-3:
                file.write('BE CARE! The Band Gap: %.6f eV is too small to tell its properties. Please check bands.dat carefully.\n' %eigs_res[1])
            else:
                file.write('This calculation is insulating with band gap %.6f eV.\n' %eigs_res[1])
        else:
            file.write('gap_up(blue): %.6f, gap_dn(red): %.6f (eV)\n' %(eigs_res[2],eigs_res[3]))

        file.write('\n')
        file.write('Finished!')

def __main__():
    band_dat()
    plotband()
    report()

__main__()
