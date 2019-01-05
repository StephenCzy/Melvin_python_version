#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:29:22 2018

@author: chenzy
"""


##This is the toobox part of the program, including experimental device functions for OAM and polarization entanglement, 
##post-selection functions and SRV functions.

##I use the "replace" function to simulate the role of the "/." function of the melvin program in mathematica.
##You can find the introduce of "replace" function in 
##https://docs.sympy.org/latest/modules/core.html?highlight=replace#sympy.core.basic.Basic.replace
##But its execution efficiency is much lower than the "/." function and I will try to speed up this code.
##(I found it when I used python and mathematica to run the same setup to get the final quantum state time).





import numpy as np
import sympy as sp
#from scipy.sparse import dok_matrix
#from numpy.linalg import matrix_rank
from sympy.matrices import SparseMatrix  ##Once I used some SparseMatrix module of numpy, it is strange that their speed is lower than those module in sympy.matrices
from sympy import IndexedBase, Idx, symbols, Symbol,exp,simplify,sqrt

V,H=symbols('V H',cls=Idx) ## H, V represents different polarizations. In the current post-selection function, I assume that all the final states have the same polarization. This is obviously wrong, but at the moment I only consider this situation.
zero=Symbol('zero')  ##  A symbol variable is used to make it easier for me to implement coincidence
psi=Symbol('psi')  ## A symbol represents quantum state
a,b,c,d,FF1,FF2,FF3,FF4,HH,GG1,GG2,GG3,GG4=map(IndexedBase,['a','b','c','d','FF1','FF2','FF3','FF4','HH','GG1','GG2','GG3','GG4']) ##represent Path
P,l,l1,l2,l3,l4,l5,l6,P1,P2,P3,P4,k=map(sp.Wild,['P','l','l1','l2','l3','l4','l5','l6','P1','P2','P3','P4','k']) ## Generally,l represents the OAM dimension and P represents the polarization dimension. k represents the number of coefficients
st2=simplify(1/sqrt(2))  ## st2 equals to sqrt(2)
l_list=np.arange(-5,5,1).tolist() ## l_list define the dimension of OAM I used.
lnum=len(l_list)


## The name of the function represents the corresponding experimental device.
##BS_fun and PBS_fun is created to facilitate the definition of the BS and PBS function.
def Reflection(psi,a): 
    psi0=psi.replace(a[l,P],lambda l,P: -1j*a[-l,H] if P==H else 1j*a[-l,V] )
    return psi0

def BS_fun(expr,a,b):
    if expr.base==a:return expr.replace(a[l,P],lambda l,P:st2*(b[l,P]+Reflection(a[l,P],a)))
    else: return expr.replace(b[l,P],lambda l,P:st2*(a[l,P]+Reflection(b[l,P],b)))
    
def BS(psi,a,b):
    psi0=sp.expand(psi.replace(lambda expr: expr.base in [a,b], lambda expr: BS_fun(expr,a,b)))
    return psi0

def PBS_fun(expr,a,b):
    if expr.base==a:return expr.replace(a[l,P],lambda l,P:b[l,H] if P==H else 1j*a[-l,V])
    else: return expr.replace(b[l,P],lambda l,P:a[l,H] if P==H else 1j*b[-l,V])

def PBS(psi,a,b):
    psi0=psi.replace(lambda expr: expr.base in [a,b], lambda expr: PBS_fun(expr,a,b))
    return psi0

def HWP(psi,a):
     psi0=psi.replace(a[l,P],lambda l,P: a[l,V] if P==H else -a[l,H])
     return psi0
 
def OAMHolo(psi,a,n):
    psi0=psi.replace(a[l,P],a[l+n,P])
    return psi0

def OAMHoloSP(psi,a,n):
    psi0=psi.replace(a[l,P],st2*(a[l,P]+a[l+n,P]))
    return psi0

def DP(psi,a,n):
    psi0=psi.replace(a[l,P],lambda l,P: simplify(exp(1j*sp.pi/n*l))*Reflection(a[l,P],a))
    return psi0

def LI(psi,a,b):
    psi0=BS(Reflection(Reflection(DP(Reflection(BS(psi,a,b),a),a,1),b),b),a,b)
    return psi0

def Coincidence(psi,FFl=[a,b,c,d]):  ##Reserved state that can meet the 4-fold in FFl, generally FFl = [a, b, c, d] 
    psi0=sp.expand(psi*zero)
    psi0=psi0.replace(k*FFl[0][l1,P1]*FFl[1][l2,P2]*FFl[2][l3,P3]*FFl[3][l4,P4]*zero,k*FFl[0][l1]*FFl[1][l2]*FFl[2][l3]*FFl[3][l4])
    psi0=psi0.subs(zero,0)
    return psi0

def Trigger(psi,a,nlist): ##must be used after Coincidence. Set a Trigger in a[l],l in nlist.
    lt=sp.Wild('lt',exclude=nlist)
    psi0=psi.replace(a[lt],0)
    psi0=psi0.replace(a[l],1)
    return psi0

def Postselect(psi,FFl,nlist): ## combine the Coincidence and Trigger Function
    #psi0=sp.simplify(psi)
    psi0=psi
    psi0=Trigger(Coincidence(psi0,FFl),FFl[0],nlist)
    psi0=psi0.subs([(FFl[0],FF1),(FFl[1],FF2),(FFl[2],FF3),(FFl[3],FF4)])
    return psi0

def toHH(psi):  ## change the state form to calculation the SRV
    psi1=psi.replace(k*FF2[l1]*FF3[l2]*FF4[l3],sp.conjugate(k)*GG2[l1]*GG3[l2]*GG4[l3])
    rho0=sp.expand(psi*psi1)
    rho0=rho0.replace(k*FF2[l1]*FF3[l2]*FF4[l3]*GG2[l4]*GG3[l5]*GG4[l6],k*HH[l1,l2,l3,l4,l5,l6])
    return rho0
    
    

def partialtrace3(rho,n): ##calculate the partial trace of the n_th photon 
    
    ## First, we get the coefficents of our term ---our matrix elemen valuet of the Density Matrix.
    dictadd=rho.as_coefficients_dict()
    dict0={}
    HH0=0
    for term in dictadd:
        listmul=term.as_ordered_factors()
        coeff=dictadd[term]
        for factor in listmul:
            if factor.is_Indexed: HH0=factor 
            else: coeff=coeff*factor
        dict0[HH0]=coeff    
    ## Second, we get Density Matrix after of our partial operation and calculate the rank finally.
    S = SparseMatrix(lnum**2, lnum**2, {(0, 0): 0})
    for term in dict0:
        if term.is_Indexed:
          if term.indices[n-1]==term.indices[n+2]:
            ll=[term.indices[0],term.indices[1],term.indices[2],term.indices[3],term.indices[4],term.indices[5]]
            del(ll[n-1],ll[n+1])
            for i in range(4):
              if ll[i]<0 : ll[i]=ll[i]+lnum
            if np.abs(dict0[term])>10e-5:
                S=SparseMatrix(lnum**2, lnum**2, {(ll[0]+ll[1]*lnum,ll[2]+ll[3]*lnum):dict0[term]})+S
        else: print('erro',term)
    return S.rank()
  
def SRV(psi):   ## The function to calculation the SRV of our quantum state
    rho=toHH(psi)
    #print(rho)
    if rho==0 :return [0,0,0]
    else: return sorted([partialtrace3(rho,1),partialtrace3(rho,2),partialtrace3(rho,3)],reverse=True)
