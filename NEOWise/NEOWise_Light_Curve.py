# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 17:28:32 2021

@author: antonfj
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Manually set Spitzer values
spitzer_data = pd.DataFrame([[7.12, 5.63, 53285.0]], columns=['w1mpro','w2mpro','mjd'])
print(spitzer_data)

# Import NEOWise data
allwise_data = pd.read_csv(r'/home/antonfj/L1551/AllWISE_results_L1551_IRS_5.csv', na_values=['null'])
neowise_data = pd.read_csv(r'/home/antonfj/L1551/neowise_catalog_L1551_IRS_5.csv')

# Get W1 and W2 fluxes
allwise_fluxes = pd.DataFrame(allwise_data, columns=['ra', 'dec','w1mpro','w2mpro','mjd','dist'])
neowise_fluxes = pd.DataFrame(neowise_data, columns=['ra','dec','w1mpro','w2mpro','mjd','dist'])

ra_mean = neowise_fluxes['ra'].mean()
dec_mean = neowise_fluxes['dec'].mean()
print(ra_mean)
print(dec_mean)

ra_diff = (neowise_fluxes['ra'] - ra_mean)*60*60
dec_diff = (neowise_fluxes['dec'] - dec_mean)*60*60
print(ra_diff.mean())
print(ra_diff.std())
print(dec_diff.mean())
print(dec_diff.std())
print(np.sqrt(ra_diff.std()**2 + dec_diff.std()**2))

pos_err = pd.DataFrame(np.sqrt(ra_diff**2 + dec_diff**2), columns=['pos_err'])
neowise_fluxes = pd.concat([neowise_fluxes, pos_err], axis=1)
print(neowise_fluxes)
neowise_fluxes = neowise_fluxes.loc[lambda neowise_fluxes: neowise_fluxes['pos_err'] < 0.7, :]
print(neowise_fluxes)

all_fluxes = pd.concat([spitzer_data, allwise_fluxes, neowise_fluxes])
neowise_fluxes.plot.scatter(x='mjd', y='w1mpro')
neowise_fluxes.plot.scatter(x='mjd', y='w2mpro')