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

##This is the search part of this program, which implements a simple process of finding the target SRV state, 
##without optimization, such as parallel computing, deleting useless states, and so on. I may add these optimizations recently, but this is not my goal. 
##How to find good search criteria is the most interesting question for me.

##define some variable
a,b,c,d=map(sp.IndexedBase,['a','b','c','d'])
V,H=sp.symbols('V H',cls=sp.Idx)
psi=sp.Symbol('psi')
Possiblepath=[]     ##Search path list of (a,b), it will be used by BS,PBS and so on.
Possiblesinglepath=[]  ##Search path list of a, it will be used by HWP,Reflection and so on

fundict={"BS":toolbox.BS,"LI":toolbox.LI,"Reflection":toolbox.Reflection,"OAMHolo":toolbox.OAMHolo,"DP":toolbox.DP,"OAMHoloSP":toolbox.OAMHoloSP} ## A dictionary used for connect the string to functions we defined in toolbox
pdict={'a':a,'b':b,'c':c,'d':d,'H':H,'V':V}  ##A dictionary used for connect the string to existed variable
mydict=dict(fundict,**pdict) ## combine two dictionary


##create Possiblepath  and Possiblesinglepath
for ii in range(ord('a'), ord('f') + 1):
    Possiblesinglepath=Possiblesinglepath+[chr(ii)]
    for jj in range(ii, ord('f') + 1):
        Possiblepath = [(chr(ii),chr(jj))]+Possiblepath

## Simplification of the Possiblepath  and Possiblesinglepath
Possiblepathsmall=[('a','b'),('a','c'),('a','d'),('b','c'),('b','d'),('c','d')]
Possiblesinglepathsmall = ['a', 'b', 'c', 'd'];

## the number in DP, OAMHoloSPP and OAMHolo we can use
nHOM=5
HOM_list=np.arange(-nHOM,-1,1).tolist()+np.arange(1,nHOM,1).tolist()

def DefineActions(): ##create the strings of our function
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


def Createsetup(setup,psi): ## create whole strings from device number
    expr=str(psi)
    for element in setup:
        expr=element.replace("XXX",expr)
    return expr

def Creatinitial_state(DC):  ## create the initial state,DC is a parameter
    initial_state=a[0,0]*b[0,0]+c[0,0]*d[0,0]
    for ii in range(1,DC+1):
        initial_state=a[ii,V]*b[-ii,V]+a[-ii,V]*b[ii,V]+initial_state
    initial_state=sp.expand(initial_state**2)
    return initial_state
    
action=DefineActions()
n_action=len(action)
n_try=10    ## the number that we search
n_device=3  ## the number of device in every setup
exp_list=randint(n_action, size=(n_try,n_device)) ## create random setup 
DC=1
target_srv=(3,3,3)  ## the SRV we want to get

initial_state=Creatinitial_state(DC)



for ii in range(1):
    setup0=Createsetup([action[jj] for jj in exp_list[ii]],initial_state)
    state0=sp.expand(sp.sympify(setup0,locals=mydict))  ## get the express from strings, 'local' is necessary.
    final_state=toolbox.Postselect(state0,[a,b,c,d],[0,1])
    srv=toolbox.SRV(final_state)
    for jj in l_list:  ## try every trigger in a path
        for kk in l_list:
            final_state=toolbox.Postselect(state0,[a,b,c,d],[jj,kk])
            srv=toolbox.SRV(final_state)
            print(final_state,srv) ##print for check
            if srv==target_srv: print('get it:',setup0,final_state,[jj,kk])  ## if we find the setup, print it
    

        
             
             
        
        
    
