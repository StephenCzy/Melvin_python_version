#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 17:29:22 2018

@author: chenzy
"""

import numpy as np
import sympy as sp
#from scipy.sparse import dok_matrix
#from numpy.linalg import matrix_rank
from sympy.matrices import SparseMatrix
from sympy import IndexedBase, Idx, symbols, Symbol,exp,simplify,sqrt

V,H=symbols('V H',cls=Idx)
zero=Symbol('zero')
psi=Symbol('psi')
a,b,c,d,FF1,FF2,FF3,FF4,HH,GG1,GG2,GG3,GG4=map(IndexedBase,['a','b','c','d','FF1','FF2','FF3','FF4','HH','GG1','GG2','GG3','GG4'])
P,l,l1,l2,l3,l4,l5,l6,P1,P2,P3,P4,k=map(sp.Wild,['P','l','l1','l2','l3','l4','l5','l6','P1','P2','P3','P4','k'])
st2=simplify(1/sqrt(2))
l_list=np.arange(-5,5,1).tolist()
lnum=len(l_list)


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

def Coincidence(psi,FFl):
    psi0=sp.expand(psi*zero)
    psi0=psi0.replace(k*FFl[0][l1,P1]*FFl[1][l2,P2]*FFl[2][l3,P3]*FFl[3][l4,P4]*zero,k*FFl[0][l1]*FFl[1][l2]*FFl[2][l3]*FFl[3][l4])
    psi0=psi0.subs(zero,0)
    return psi0

def Trigger(psi,a,nlist): ##must be used after Coincidence
    lt=sp.Wild('lt',exclude=nlist)
    psi0=psi.replace(a[lt],0)
    psi0=psi0.replace(a[l],1)
    return psi0

def Postselect(psi,FFl,nlist):
    #psi0=sp.simplify(psi)
    psi0=psi
    psi0=Trigger(Coincidence(psi0,FFl),FFl[0],nlist)
    psi0=psi0.subs([(FFl[0],FF1),(FFl[1],FF2),(FFl[2],FF3),(FFl[3],FF4)])
    return psi0

def toHH(psi):
    psi1=psi.replace(k*FF2[l1]*FF3[l2]*FF4[l3],sp.conjugate(k)*GG2[l1]*GG3[l2]*GG4[l3])
    rho0=sp.expand(psi*psi1)
    rho0=rho0.replace(k*FF2[l1]*FF3[l2]*FF4[l3]*GG2[l4]*GG3[l5]*GG4[l6],k*HH[l1,l2,l3,l4,l5,l6])
    return rho0
    
    

def partialtrace3(rho,n):
    dictadd=rho.as_coefficients_dict()
    dict0={}
    HH0=0
    for term in dictadd:
        listmul=term.as_ordered_factors()
        #print(listmul)
        coeff=dictadd[term]
        for factor in listmul:
            if factor.is_Indexed: HH0=factor 
            else: coeff=coeff*factor
        dict0[HH0]=coeff     
    #print(dict0)
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
  
def SRV(psi):
    rho=toHH(psi)
    #print(rho)
    if rho==0 :return [0,0,0]
    else: return sorted([partialtrace3(rho,1),partialtrace3(rho,2),partialtrace3(rho,3)],reverse=True)
    
'''
psi=sp.expand((a[0,V]*b[0,V]+a[1,V]*b[-1,V]+a[-1,V]*b[1,V]+c[0,V]*d[0,V]+c[1,V]*d[-1,V]+c[-1,V]*d[1,V])**2)
#psi=a[0, V]*b[0, V]*c[0, V]*d[0, V] + a[0, V]*b[0, V]*c[1, V]*d[-1, V] + a[0, V]*b[0, V]*c[-1, V]*d[1, V] + a[1, V]*b[-1, V]*c[0, V]*d[0, V] +a[-1, V]*b[1, V]*c[0, V]*d[0, V] + a[1, V]*b[-1, V]*c[1, V]*d[-1, V] + a[1, V]*b[-1, V]*c[-1, V]*d[1, V] + a[-1, V]*b[1, V]*c[1, V]*d[-1, V] + a[-1, V]*b[1, V]*c[-1, V]*d[1, V]
r1=LI(psi,b,c)

r2=Reflection(r1,a)
r3=OAMHolo(r2,a,-2)
r4=BS(r3,a,c)
print(r4)

psi0=a[-1, V]**2*b[1, V]**2 + 2*a[-1, V]*a[0, 0]*b[0, 0]*b[1, V] + 2*a[-1, V]*a[1, V]*b[-1, V]*b[1, V] + 2*a[-1, V]*b[1, V]*c[0, 0]*d[0, 0] + a[0, 0]**2*b[0, 0]**2 + 2*a[0, 0]*a[1, V]*b[-1, V]*b[0, 0] + 2*a[0, 0]*b[0, 0]*c[0, 0]*d[0, 0] + a[1, V]**2*b[-1, V]**2 + 2*a[1, V]*b[-1, V]*c[0, 0]*d[0, 0] + c[0, 0]**2*d[0, 0]**2
r1=DP(psi0,d,1)
print(r1)
r2=OAMHoloSP(r1,a,1)
print(r2)
r4=OAMHoloSP(LI(OAMHoloSP(BS(r2,c,d),a,1),a,c),b,1)
print(r4)
r5=Postselect(r4,[a,b,c,d],[1,0])
print(r5)

print(SRV(r5))

'''