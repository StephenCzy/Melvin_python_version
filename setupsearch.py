#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 18:48:48 2018

@author: chenzy
"""
import numpy as np
from numpy.random import randint
import toolbox_2018_12_26 as toolbox
import sympy as sp


a,b,c,d=map(sp.IndexedBase,['a','b','c','d'])
V,H=sp.symbols('V H',cls=sp.Idx)
psi=sp.Symbol('psi')
Possiblepath=[]
Possiblesinglepath=[]

fundict={"BS":toolbox.BS,"LI":toolbox.LI,"Reflection":toolbox.Reflection,"OAMHolo":toolbox.OAMHolo,"DP":toolbox.DP,"OAMHoloSP":toolbox.OAMHoloSP}
pdict={'a':a,'b':b,'c':c,'d':d,'H':H,'V':V}
mydict=dict(fundict,**pdict)

for ii in range(ord('a'), ord('f') + 1):
    Possiblesinglepath=Possiblesinglepath+[chr(ii)]
    for jj in range(ii, ord('f') + 1):
        Possiblepath = [(chr(ii),chr(jj))]+Possiblepath

Possiblepathsmall=[('a','b'),('a','c'),('a','d'),('b','c'),('b','d'),('c','d')]
Possiblesinglepathsmall = ['a', 'b', 'c', 'd'];
l_list=np.arange(-5,5,1).tolist()
nHOM=5
HOM_list=np.arange(-nHOM,-1,1).tolist()+np.arange(1,nHOM,1).tolist()


def DefineActions():
    action=[]
    for pp in Possiblepathsmall:
        action=action+["BS(XXX,"+pp[0]+","+pp[1]+")"]
        action=action+["LI(XXX,"+pp[0]+","+pp[1]+")"]
    for psp in Possiblesinglepathsmall:
        action=action+["Reflection(XXX,"+psp+")"]
        for ii in HOM_list:
             action=action+["OAMHolo(XXX,"+psp+","+str(ii)+")"]
             action=action+["DP(XXX,"+psp+","+str(ii)+")"]
             action=action+["OAMHoloSP(XXX,"+psp+","+str(ii)+")"]
    return action

def Createsetup(setup,psi):
    expr=str(psi)
    for element in setup:
        expr=element.replace("XXX",expr)
    return expr

def Creatinitial_state(DC):
    initial_state=a[0,0]*b[0,0]+c[0,0]*d[0,0]
    for ii in range(1,DC+1):
        initial_state=a[ii,V]*b[-ii,V]+a[-ii,V]*b[ii,V]+initial_state
    initial_state=sp.expand(initial_state**2)
    return initial_state
    
action=DefineActions()
n_action=len(action)
n_try=10
n_device=3
exp_list=randint(n_action, size=(n_try,n_device))
DC=1
target_srv=(3,3,3)

initial_state=Creatinitial_state(DC)



for ii in range(n_try):
    setup0=Createsetup([action[jj] for jj in exp_list[ii]],initial_state)
    state0=sp.expand(sp.sympify(setup0,locals=mydict))
    final_state=toolbox.Postselect(state0,[a,b,c,d],[0,1])
    srv=toolbox.SRV(final_state)
    for jj in l_list:
        for kk in l_list:
            final_state=toolbox.Postselect(state0,[a,b,c,d],[jj,kk])
            srv=toolbox.SRV(final_state)
            print(final_state,srv)
            if srv==target_srv: print('get it:',setup0,final_state,[jj,kk])
    

        
             
             
        
        
    
