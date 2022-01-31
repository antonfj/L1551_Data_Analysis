# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 17:28:32 2021

@author: antonfj
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans

# Manually set Spitzer values
spitzer_data = pd.DataFrame([[7.12, 5.63, 53285.0]], columns=['w1mpro','w2mpro','mjd'])
print(spitzer_data)

# Import AllWise and NEOWise data
allwise_data = pd.read_csv(r'/home/antonfj/L1551/L1551_Data_Analysis/NEOWise/AllWISE_results_L1551_IRS_5.csv', na_values=['null'])
neowise_data = pd.read_csv(r'/home/antonfj/L1551/L1551_Data_Analysis/NEOWise/neowise_catalog_L1551_IRS_5.csv')

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

# Remove measurements where pos. was more than 0.7 arcsec from the mean pos.
pos_err = pd.DataFrame(np.sqrt(ra_diff**2 + dec_diff**2), columns=['pos_err'])
neowise_fluxes = pd.concat([neowise_fluxes, pos_err], axis=1)
print(neowise_fluxes)
neowise_fluxes = neowise_fluxes.loc[lambda neowise_fluxes: neowise_fluxes['pos_err'] < 0.7, :]
print(neowise_fluxes)

all_fluxes = pd.concat([spitzer_data, allwise_fluxes, neowise_fluxes], ignore_index=True)

kmeans_m1 = KMeans(n_clusters=17).fit(all_fluxes[['w1mpro','mjd']])
centroids_m1 = kmeans_m1.cluster_centers_
epochs_m1 = pd.DataFrame(np.transpose(kmeans_m1.labels_), columns=['epoch'])
all_fluxes = pd.concat([all_fluxes, epochs_m1], axis=1)
print(all_fluxes)

#kmeans_m2 = KMeans(n_clusters=17).fit(all_fluxes[['w2mpro','mjd']])
#centroids_m2 = kmeans_m2.cluster_centers_

w1_means = np.array([])
w1_stds = np.array([])
mjd_w1_means = np.array([])
w2_means = np.array([])
w2_stds = np.array([])
mjd_w2_means = np.array([])

for i in range(0,17):
    epoch = all_fluxes.loc[all_fluxes['epoch'] == i]
    
    if len(epoch) == 1:
        w1_means = np.append(w1_means, epoch['w1mpro'])
        w1_stds = np.append(w1_stds, 0)
        mjd_w1_means = np.append(mjd_w1_means, epoch['mjd'])
    else:
        # Remove measurements outside 15th and 85th percentiles
        epoch = epoch.loc[(epoch['w1mpro'] > epoch['w1mpro'].quantile(.15))]    
        epoch = epoch.loc[(epoch['w1mpro'] < epoch['w1mpro'].quantile(.85))] 
        w1_mean = epoch['w1mpro'].mean()
        w1_std = epoch['w1mpro'].std()
        mjd_mean = epoch['mjd'].mean()
        w1_means = np.append(w1_means, w1_mean)
        w1_stds = np.append(w1_stds, w1_std)
        mjd_w1_means = np.append(mjd_w1_means, mjd_mean)

    # Repeat for w2 flux measurements
    epoch = all_fluxes.loc[all_fluxes['epoch'] == i]
    
    if len(epoch) == 1:
        w2_means = np.append(w2_means, epoch['w2mpro'])
        w2_stds = np.append(w2_stds, 0)
        mjd_w2_means = np.append(mjd_w2_means, epoch['mjd'])
    else:
        # Remove measurements outside 15th and 85th percentiles
        epoch = epoch.loc[(epoch['w2mpro'] > epoch['w2mpro'].quantile(.15))]    
        epoch = epoch.loc[(epoch['w2mpro'] < epoch['w2mpro'].quantile(.85))] 
        w2_mean = epoch['w2mpro'].mean()
        w2_std = epoch['w2mpro'].std()
        mjd_mean = epoch['mjd'].mean()
        w2_means = np.append(w2_means, w2_mean)
        w2_stds = np.append(w2_stds, w2_std)
        mjd_w2_means = np.append(mjd_w2_means, mjd_mean)

#all_fluxes.plot.scatter(x='mjd', y='w1mpro', ylim=(8.0, 7.0))
plt.figure(0)
plt.ylim(8.0, 7.0)
#plt.scatter(centroids_m1[:, 1], centroids_m1[:, 0], c='red', s=50, ylim=(8.0, 7.0))
plt.errorbar(mjd_w1_means, w1_means, yerr=w1_stds, color='blue', fmt='o')

#all_fluxes.plot.scatter(x='mjd', y='w2mpro', ylim=(6.0, 4.4))
plt.figure(1)
plt.ylim(6.0, 4.4)
plt.errorbar(mjd_w2_means, w2_means, yerr=w2_stds, color='blue', fmt='o')
