# -*- coding: utf-8 -*-
"""
@author: eddybrandon
"""

import pylab
import numpy as np
import matplotlib.pyplot as plt

AU=149597870700

NPlaneta='Jupiter'		#Cuerpo a analizar (Júpiter en este ejemplo)
hh=1                
NP= 5000*hh
data = pylab.loadtxt('RK{}Mond.dat'.format(NPlaneta), dtype=float)  	#Archivo de datos con potencial modificado.
jdata = pylab.loadtxt('JPL{}.dat'.format(NPlaneta), dtype=float)		#Archivo de datos con información de HORIZONS JPL.
ndata = pylab.loadtxt('RK{}Newton.dat'.format(NPlaneta), dtype=float)	#Archivo de datos con potencial newtoniano.

xm = data[:,0]
ym = data[:,1]
zm = data[:,2]
xn = ndata[:,0]
yn = ndata[:,1]
zn = ndata[:,2]
xj = jdata[:,0]*1000/AU
yj = jdata[:,1]*1000/AU
zj = jdata[:,2]*1000/AU

#MOND:
rminM=rminN=100    #Para el perihelio
rmaxM=rmaxN=0      #Para el afelio
perM=perN=10      #Para el periodo
diasM=diasN=0      #Num. de días
#JPL:
rminj=100    #Para el perihelio
rmaxj=0      #Para el afelio
perj=10     #Para el periodo
diasj=0      #Num. de días

for i in range(NP):
    #Lectura de los datos de salida del método numérico
    dx=xm[i]
    dy=ym[i]
    dz=zm[i]
    px=xm[i]-xm[0]
    py=ym[i]-ym[0]
    pz=zm[i]-zm[0]
    ndx=xn[i]
    ndy=yn[i]
    ndz=zn[i]
    npx=xn[i]-xn[0]
    npy=yn[i]-yn[0]
    npz=zn[i]-zn[0]
    jdx=xj[i]
    jdy=yj[i]
    jdz=zj[i]
    jpx=(xj[i]-xj[0])
    jpy=(yj[i]-yj[0])
    jpz=(zj[i]-zj[0])
    rM=np.sqrt(dx**2+dy**2+dz**2)
    rpM=np.sqrt(px**2+py**2+pz**2)
    rN=np.sqrt(ndx**2+ndy**2+ndz**2)
    rpN=np.sqrt(npx**2+npy**2+npz**2)
    rj=np.sqrt(jdx**2+jdy**2+jdz**2)
    rpj=np.sqrt(jpx**2+jpy**2+jpz**2)
    if rM<rminM:
        rminM=rM
    elif rM>rmaxM:
        rmaxM=rM
    if i>100:
        if rpM<perM:
            perM=rpM
            diasM=i
    if rN<rminN:
        rminN=rN
    elif rN>rmaxN:
        rmaxN=rN
    if i>100:
        if rpN<perN:
            perN=rpN
            diasN=i
    if rj<rminj:
        rminj=rj
    elif rj>rmaxj:
        rmaxj=rj
    if i>100:
        if rpj<perj:
            perj=rpj
            diasj=i
            
Datos=open("JPL{}RK.dat".format(NPlaneta), "w")

### Mond
saM=(rminM+rmaxM)/2
eccM=(rmaxM-rminM)/(rminM+rmaxM)
Datos.write('{} \n \n'.format(NPlaneta)+'MOND: \n \n'+'a={} [UA]; a={} [km]'.format(saM, saM*AU)+'ecc= {}'.format(eccM))      
Datos.write('\n'+ 'Perihelio: en [UA]={}; en [km]={} \n'.format(rminM, rminM*AU))
Datos.write('Afelio: en [UA]={}; en [km]={} \n'.format(rmaxM, rmaxM*AU))
Datos.write('Período: {} días; {} años; ({} [UA]) \n \n'.format(diasM/hh, diasM/365.2425, perM))

### Newton
saN=(rminN+rmaxN)/2
eccN=(rmaxN-rminN)/(rminN+rmaxN)
Datos.write('NEWTON: \n \n'+'a={} [UA]; a={} [km]'.format(saN, saN*AU)+'ecc= {}'.format(eccN))      
Datos.write('\n'+ 'Perihelio: en [UA]={}; en [km]={} \n'.format(rminN, rminN*AU))
Datos.write('Afelio: en [UA]={}; en [km]={} \n'.format(rmaxN, rmaxN*AU))
Datos.write('Período: {} días; {} años; ({} [UA]) \n \n'.format(diasN/hh, diasN/365.2425, perN))

### JPL
saj=(rminj+rmaxj)/2
eccj=(rmaxj-rminj)/(rminj+rmaxj)
Datos.write('JPL: \n \n'+'a={} [UA]; a={} [km]'.format(saj, saj*AU)+'ecc= {}'.format(eccj))      
Datos.write('\n'+ 'Perihelio: en [UA]={}; en [km]={} \n'.format(rminj, rminj*AU))
Datos.write('Afelio: en [UA]={}; en [km]={} \n'.format(rmaxj, rmaxj*AU))
Datos.write('Período: {} días; {} años; ({} [UA]) \n \n'.format(diasj/hh, diasj/(hh*365.2425), perj))

##############
#Errores porcentuales Newton vs JPL
errAN=100*(saj-saN)/saj
erreN=100*(eccj-eccN)/eccj
errPerN=100*(diasj-diasN)/diasj
Datos.write('Diferencias Newton vs. JPL: \n \n'+'D a = {} % \n'.format(errAN))    
Datos.write('D ecc = {} % \n'.format(erreN))  
Datos.write('D Per = {} % \n \n'.format(errPerN))  

#Errores porcentuales Mond vs JPL
errA=100*(saj-saM)/saj
erre=100*(eccj-eccM)/eccj
errPer=100*(diasj-diasM)/diasj
Datos.write('Diferencias MOND vs. JPL: \n \n'+'D a = {} % \n'.format(errA))    
Datos.write('D ecc = {} % \n'.format(erre))  
Datos.write('D Per = {} % \n \n'.format(errPer))

#Errores porcentuales Mond vs Newton
errAMN=100*(saN-saM)/saN
erreMN=100*(eccN-eccM)/eccN
errPerMN=100*(diasN-diasM)/diasN
Datos.write('Diferencias MOND vs. NEWTON: \n \n'+'D a = {} % \n'.format(errAMN))    
Datos.write('D ecc = {} % \n'.format(erreMN))  
Datos.write('D Per = {} % \n'.format(errPerMN))   

Datos.closed
