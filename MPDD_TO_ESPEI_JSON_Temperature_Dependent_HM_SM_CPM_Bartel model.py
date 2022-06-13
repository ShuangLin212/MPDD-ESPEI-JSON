#!/usr/bin/env python
# coding: utf-8

# In[2]:


import math
from itertools import combinations
from pymongo import MongoClient
import os
from pymatgen.analysis.phase_diagram import PDEntry, PDPlotter, PhaseDiagram
import pymatgen.core as mg
import json
import math
import ssl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# In[3]:


# Connect with your MPDD Mongodb which contains temperature dependent data
client_string = 'mongodb+srv://shuangMPDD:FRWUXrZ8DhKGgCAB@mpdd.hklel.mongodb.net/MPDD?retryWrites=true&w=majority'
client = MongoClient(client_string)
MPDDshuangTest = client['MPDD']['shuangTest']


# In[5]:


def subsystems(comp):
    eList = comp.chemical_system.split('-')
    out = []
    for nComponents in range(eList.__len__()):
        for c in combinations(eList, 1+nComponents):
            t = list(c)
            t.sort()
            out.append('-'.join(t))
    return out


# In[6]:


# Define the subsystem 
sysList = subsystems(mg.Composition('NbMoW'))
print(sysList)


# In[34]:


# Enthalpy with anonymizedFormula : 'ABC'
for e in client['MPDD']['shuangTest'].find(
        {'material.system': {'$in': sysList},
         'material.anonymizedFormula': 'ABC',
         }
        ):
    reducedFormula = e['material']['reducedFormula']
    gSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['G_SISSO']
    sSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['S_SISSO']
    temperatures = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['T_SISSO']
    mpddid = str(e['_id'])
    TS = []
    for i in range(0, len(sSISSO)):
        TS.append(sSISSO[i] * temperatures[i])
    H_ev = gSISSO+TS
    H = [i*96491.5666 for i in H_ev] #Unit change from eV/atom to J/mol, 1 eV/atom=96491.5666 J/mol
    #print(H)
    dictionary ={
    "components":["NB", "Mo", "W"],
    "phases":reducedFormula,
    "solver":{
	  "mode": "manual",
	  "sublattice_site_ratios": [1],
	  "sublattice_configurations": [["NB", "MO", "W"]] 
  },
  "conditions": {
	  "P": 101325,
	  "T": temperatures,
  },
  "output": 'HMFORM',
  "values":   [[[H]]],
  "reference": "MPDDid:"+mpddid,
  "comment":"Bartel Model"
            }
    with open(mpddid+"_HM_MPDD_Bartel_Model.json", "w") as outfile:
        json.dump(dictionary, outfile,indent=1)


# In[31]:


#Entropy
for e in client['MPDD']['shuangTest'].find(
        {'material.system': {'$in': sysList},
         'material.anonymizedFormula': 'ABC',
         }
        ):
    reducedFormula = e['material']['reducedFormula']
    gSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['G_SISSO']
    sSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['S_SISSO']
    cpSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['Cp_SISSO']
    temperatures = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['T_SISSO']
    mpddid = str(e['_id'])
    s = [i*96491.5666 for i in sSISSO]
    dictionary ={
    "components":["NB", "Mo", "W"],
    "phases":reducedFormula,
    "solver":{
	  "mode": "manual",
	  "sublattice_site_ratios": [1],
	  "sublattice_configurations": [["NB", "MO", "W"]] 
  },
  "conditions": {
	  "P": 101325,
	  "T": temperatures,
  },
  "output": 'SM',
  "values":   [[[s]]],
  "reference": "MPDDid:"+mpddid,
  "comment":"Bartel Model"
            }
    with open(mpddid+"_SM_MPDD_Bartel_Model.json", "w") as outfile:
        json.dump(dictionary, outfile,indent=1)


# In[33]:


#heat capacity
for e in client['MPDD']['shuangTest'].find(
        {'material.system': {'$in': sysList},
         'material.anonymizedFormula': 'ABC',
         }
        ):
    reducedFormula = e['material']['reducedFormula']
    gSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['G_SISSO']
    sSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['S_SISSO']
    cpSISSO = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['Cp_SISSO']
    temperatures = e['properties']['GFT_SIPFENN_NMM+Bartel2018']['T_SISSO']
    mpddid = str(e['_id'])
    CP = [i*96491.5666 for i in cpSISSO]
    dictionary ={
    "components":["NB", "Mo", "W"],
    "phases":reducedFormula,
    "solver":{
	  "mode": "manual",
	  "sublattice_site_ratios": [1],
	  "sublattice_configurations": [["NB", "MO", "W"]] 
  },
  "conditions": {
	  "P": 101325,
	  "T": temperatures,
  },
  "output": 'CPM',
  "values":   [[[CP]]],
  "reference": "MPDDid:"+mpddid,
  "comment":"Bartel Model"
            }
    with open(mpddid+"_CPM_MPDD_Bartel_Model.json", "w") as outfile:
        json.dump(dictionary, outfile,indent=1)


# In[ ]:





# In[ ]:




