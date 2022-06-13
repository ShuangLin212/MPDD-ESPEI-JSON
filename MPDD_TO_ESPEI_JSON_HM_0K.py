#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


# Connect with your MPDD Mongodb
MPDD_client_string = 'mongodb+srv://MPDD_Client:mVjfZG8LIySUpxny@mpdd.hklel.mongodb.net/MPDD?retryWrites=true&w=majority'
MPDDclient = MongoClient(MPDD_client_string)
MPDDcurated = MPDDclient['MPDD']['curated']


# In[3]:


def subsystems(comp):
    eList = comp.chemical_system.split('-')
    out = []
    for nComponents in range(eList.__len__()):
        for c in combinations(eList, 1+nComponents):
            t = list(c)
            t.sort()
            out.append('-'.join(t))
    return out


# In[4]:


#Define a subsystem you are interest in
sysList = subsystems(mg.Composition('NbW'))

for e in MPDDcurated.find(
        {'material.system': {'$in': sysList},
         'properties.Stability_SIPFENN_AMK2020_NMM': 0},
        {'material.reducedFormula': 1, 'metadata.parentDatabaseID': 1,
         'properties.FormationEnergy_SIPFENN_Krajewski2020_NovelMaterialsModel': 1}):

    formula = e['material']['reducedFormula']
    dH = e['properties']['FormationEnergy_SIPFENN_Krajewski2020_NovelMaterialsModel']


# In[5]:


# 
for e in MPDDcurated.find(
        {'material.system': {'$in': sysList},
         'material.anonymizedFormula': 'AB',
         'material.spaceGroupN': 225},
        {'material.reducedFormula': 1,
         'metadata.parentDatabaseID': 1,
         'properties.FormationEnergy_SIPFENN_Krajewski2020_NovelMaterialsModel': 1,
         'material.POSCAR': 1}
        ):

    formula = e['material']['reducedFormula']
    dH = e['properties']['FormationEnergy_SIPFENN_Krajewski2020_NovelMaterialsModel']
    POSCAR = e['material']['POSCAR']
    pid = e['metadata']['parentDatabaseID']
    mpddid = str(e['_id'])
    HM_FORM_0K= dH * 96491.5666 #Unit change from eV/atom to J/mol, 1 eV/atom=96491.5666 J/mol
    print(f'{formula} - {HM_FORM_0K} - {mpddid} - {pid}\n{POSCAR}\n\n')
    first=formula[0]+formula[1] # modify here in terms of character number element
    second=formula[2]
    dictionary ={
  "components":[first, second],
  "phases":[formula],
  "solver":{
	  "mode": "manual",
	  "sublattice_site_ratios": [0.5, 0.5],
	  "sublattice_configurations": [[first, second]] 
  },
  "conditions": {
	  "P": 101325,
	  "T": 0,
  },
  "output": "HM_FORM",
  "values":   [[[HM_FORM_0K]]],
  "reference": "MPDDid:"+mpddid
  #"comment":
            }
  
    with open(formula+"_HM_FORM_0K_"+'_'+mpddid+".json", "w") as outfile:
        json.dump(dictionary, outfile,indent=1)


# In[ ]:




