#!/usr/bin/env python
# coding: utf-8
"""
should be run from examples folder as current directory
downloads, processes W cascades of 2keV from cascadesDB
shows some basic plots from the processed data
runs a http server and opens the web-app with processed data loaded
"""

# In[ ]:


import sys
import shutil
import os
import json
import sqlite3
import pandas as pd
import seaborn as sns
from pandas import DataFrame
import pickle
import numba
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.neighbors import NearestNeighbors
import umap
import logging

pathToAnuvikar = "."
sys.path.append(pathToAnuvikar)
from pysrc.anuvikar_cdb_helper import getDefaultConfig, writeResultsToJSON
buildDirs = []
buildDirs.append(os.path.join(pathToAnuvikar, "_build"))
buildDirs.append(os.path.join(pathToAnuvikar, "build"))
libPaths = []
libPaths.append(os.path.join(buildDirs[0], "libanuvikar_shared.so"))
libPaths.append(os.path.join(buildDirs[0], "libanuvikar_shared.dylib"))
libPaths.append(os.path.join(buildDirs[1], "libanuvikar_shared.so"))
libPaths.append(os.path.join(buildDirs[1], "libanuvikar_shared.dylib"))
libPath = ''
for _libPath in libPaths:
  if (os.path.exists(_libPath)):
    libPath = _libPath
if libPath == '':
  print("Library not found at", libPath)
  print("This might be due to build errors in cmake.")
  print("If built successfully, edit this source and correct build directory & lib file (so / dylib / dll) path.")
  exit(1)

def dist(a, b):
    res = 0.0
    for x, y in zip(a, b):
        if (abs(x) > 1e-6):
            res += ((x - y)**2 * 1.0) / (1.0*x)
    return round(res, 4)

"""
Helper distance function for dimensionality reduction
"""
@numba.njit()
def chiSqr(x, y, startA, startB, endA, endB):  # brat_curtis
    numerator = 0.0
    denominator = 0.0
    for i, j in zip(range(startA, endA), range(startB, endB)):
        numerator += np.abs(x[i] - y[j])
        denominator += np.abs(x[i] + y[j])

    if denominator > 0.0:
        return float(numerator) / denominator
    else:
        return 0.0



"""
Distance function for dimensionality reduction
"""
@numba.njit()
def quad(x, y):
    l = x.shape[0]
    a = chiSqr(x, y, 0, 0, 36, 36)
    d = chiSqr(x, y, 36, 36, l, l)
    preA = chiSqr(x, y, 0, 1, 35, 36)
    postA = chiSqr(x, y, 1, 0, 36, 35)
    preD = chiSqr(x, y, 36, 37, l - 1, l)
    postD = chiSqr(x, y, 37, 36, l, l - 1)
    wA = 1.2
    wD = 0.9
    wAs = 0.4
    wDs = 0.25
    cA = (wAs * (preA + postA) + a) * wA / (2.0 * wAs + 1.0)
    cD = (wDs * (preD + postD) + d) * wD / (2.0 * wDs + 1.0)
    return (cA + cD) / (wA + wD)

def clusterClassData(data):
    feat = []
    tag = []
    for i, x in enumerate(data):
        for y in x['features']:
            #feat.append(x['features'][y]['angle'] + x['features'][y]['dist'])
            feat.append(x['features'][y]['angle'] + x['features'][y]['dist'])
            tag.append((x['id'], y, i))
    return (feat, tag)

def quadCustom(wA, wD):
    def quad(x, y):
        l = x.shape[0]
        a = chiSqr(x, y, 0, 0, 36, 36)
        d = chiSqr(x, y, 36, 36, l, l)
        preA = chiSqr(x, y, 0, 1, 35, 36)
        postA = chiSqr(x, y, 1, 0, 36, 35)
        preD = chiSqr(x, y, 36, 37, l - 1, l)
        postD = chiSqr(x, y, 37, 36, l, l - 1)
        wAs = 0.4
        wDs = 0.25
        cA = (wAs * (preA + postA) + a) * wA / (2.0 * wAs + 1.0)
        cD = (wDs * (preD + postD) + d) * wD / (2.0 * wDs + 1.0)
        return (cA + cD) / (wA + wD)
    return quad

# if old's nearest is from new then update else let it be
# add new ones with it.
# new dimensionality reduction :|
def cookNewComparison(oldFt, feat, tag):
  topsize = 5
  neigh = {}
  keys = ['angle', 'dist', 'all']
  quadAngle = quadCustom(1.0, 0.0)
  quadDist = quadCustom(0.0, 1.0)
  quadBoth = quad
  defaultK = topsize * 3 if topsize * 3 < len(feat) else len(feat) - 1
  neigh[keys[0]] = NearestNeighbors(n_neighbors = defaultK, metric=quadAngle)
  neigh[keys[1]] = NearestNeighbors(n_neighbors = defaultK, metric=quadDist)
  neigh[keys[2]] = NearestNeighbors(n_neighbors = defaultK, metric=quadBoth)
  dists = {}
  neighbours = {}
  allFeat = oldFt['feat'] + feat
  for key in neigh:
    if len(feat) == 0: continue
    neigh[key].fit(allFeat)
    dists[key], neighbours[key] = neigh[key].kneighbors()
  allTags = oldFt['tag'] + tag
  oldLen = len(oldFt['tag'])
  additions = {}
  refs = {}
  # add old
  for index, tagv in enumerate(oldFt['tag']):
    isUpdate = False
    vals = {}
    curRef = {}
    for key in neigh:
      vals[key] = []
      curRef[key] = []
      for x, y in  zip(dists[key][index][:topsize], neighbours[key][index][:topsize]):
        if y > oldLen: isUpdate = True
        vals[key].append((round(x, 2), allTags[y][0], allTags[y][1]))
        curRef[key].append((allTags[y][0], allTags[y][1], allTags[y][2] if (y > oldLen - 1) else -1))
    if isUpdate:
      additions[tagv] = vals
      refs[tagv] = curRef
  # add new
  for index, tagv in enumerate(tag):
    totalIndex = index + oldLen
    additions[tagv] = {}
    refs[tagv] = {}
    for key in neigh:
      additions[tagv][key] = [(round(x,2), allTags[y][0], allTags[y][1]) for x, y in zip(
                dists[key][totalIndex][:topsize], neighbours[key][totalIndex][:topsize])]
      refs[tagv][key] = [(allTags[y][0], allTags[y][1], allTags[y][2] if (y > oldLen - 1) else -1) for y in neighbours[key][totalIndex][:topsize]]
  return additions, refs, allFeat, allTags
  ##curLen = oldNN['size1'][index]
  ##lenDiff = [(abs(curLen - len(data[tag[x][0]]['clusters'][tag[x][1]])), i)
  ##             for i, x in enumerate(neighbours[key][index])]
  ##  lenDiff.sort()
  ##  cascade['clust_cmp_size'][cid][key] = [
  ##      (dists[key][index][x[1]], tag[neighbours[key][index][x[1]]][0], tag[neighbours[key][index][x[1]]][1]) for x in lenDiff[:topsize]]

def mergeCascadeDbs(dbNew, dataPath, dest, dbOld = None, oldFtPath=None):
  if (not os.path.exists(dataPath)): return (False, "Can not access: " + dataPath)
  dataFile = open(dataPath, 'r')
  data = json.load(dataFile)
  dataFile.close()
  feat, tag = clusterClassData(data)
  oldFt = {"feat":[], "tag":[], "reducer": None}
  rndSeed = 42
  #reducer = umap.UMAP(n_components=2, n_neighbors=6, min_dist=0.45, metric=quad, random_state=rndSeed)
  reducer = umap.UMAP(n_components=2, n_neighbors=6, min_dist=0.45, random_state=rndSeed)
  dims = None
  if (oldFtPath and os.path.exists(oldFtPath)):
    ftFile = open(oldFtPath, 'rb')
    oldFt = pickle.load(ftFile)
    reducer = oldFt['reducer']
    dims = reducer.transform(feat).tolist()
    ftFile.close()
  else:
    dims = reducer.fit_transform(feat).tolist()
  additions, addrefs, allFeat, allTags = cookNewComparison(oldFt, feat, tag)
  #reducer = umap.UMAP(n_components=2, n_neighbors=6, min_dist=0.45, metric=quad, random_state=rndSeed).fit(allFeat)
  reducer = umap.UMAP(n_components=2, n_neighbors=6, min_dist=0.45, random_state=rndSeed).fit(allFeat)
  #print(updates)
  #print(additions['all'])
  #print(allTags)
  #print(allTags.keys())
  saveNN(data, {'feat':allFeat, 'tag':allTags, "reducer":reducer}, len(feat), len(oldFt['feat']), dest)
  #updateComparison(dbNew, additions)
  #print(dbNew)
  if dbOld: shutil.copy(dbOld, dest)
  else: shutil.copy(dbNew, dest)
  con = sqlite3.connect(dest)
  cur = con.cursor()
  if dbOld: addToDb(dbNew, cur, con)
  newTags = set(tag)
  for key in additions:
    cascadeid = key[0]#data[key[0]]['id']
    name = key[1]
    valJson = json.dumps(additions[key])
    pairsJson = getComparisonPairs(data, addrefs[key], cur)
    if key in newTags:
      dim = dims[key[2]]
      cur.execute("UPDATE clusters set cmp= ?, cmpsize=?, cmppairs = ? where cascadeid = ? and name = ?", (valJson, valJson, pairsJson, cascadeid, name))
      #cur.execute("UPDATE clusters set cmp= ?, cmpsize=?, cmppairs = ?, hdbx=?, hdby=? where cascadeid = ? and name = ?", (valJson, valJson, pairsJson, round(dim[0], 2), round(dim[1], 2), cascadeid, name))
      #cur.execute("UPDATE clusters set cmp= ?, cmpsize=?, cmppairs = ? where cascadeid = ? and name = ?", (valJson, valJson, pairsJson, cascadeid, name))
    else:
      cur.execute("UPDATE clusters set cmp= ?, cmpsize=?, cmppairs = ? where cascadeid = ? and name = ?", (valJson, valJson, pairsJson, cascadeid, name))
  con.commit()
  for (t, z) in zip(tag, dims):
    cur.execute("UPDATE clusters set hdbx=?, hdby=? where cascadeid = ? and name = ?", (z[0], z[1], t[0], t[1]))
  con.commit()
  con.close()
  return (True, "")
  #store
  #updateComparison(db2Path, updates)
  #cpDb(db1Path, db2Path)

def addToDb(dbNew, cur, con):
  cur.execute("Attach ? as nu", (dbNew, ))
  cur.execute("BEGIN")
  cols = 'cascadeid, name, savimorph, size, coordtype, coords, hdbx, hdby, cmp, cmpsize, cmppairs, morphdesc, properties'
  cur.execute("INSERT INTO cascades select * from nu.cascades")
  q1 = "INSERT INTO clusters ("+ cols +") select "+ cols + " from nu.clusters"
  cur.execute(q1)
  con.commit()
  cur.execute("detach database nu")

def saveNN(data, neigh, newLen, oldLen, basepath):
  i = 0
  while i < newLen:
    x = neigh['tag'][oldLen + i]
    neigh['tag'][oldLen + i] = (x[0], x[1], data[x[2]]['clusterSizes'][x[1]])
    i += 1
  f = open(basepath +"_tree.pickle", "wb")
  pickler = pickle.Pickler(f)
  pickler.dump(neigh)
  f.close()

def cookCmpCascadeInfo(row):
  return {"id":row["id"], "substrate": row["substrate"], "energy": row["energy"], "temperature": row["temperature"], "potentialused": row["potentialUsed"], "author": row["author"]};

def cookCmpCascadeInfoFromDb(cascadeid, cur):
  cols = ["id", "substrate", "energy", "temperature", "potentialused", "author"]
  q = "Select " + ", ".join(cols)  + " from cascades WHERE id =?"
  row = cur.execute(q, (cascadeid, )).fetchone()
  res = {key:val for (key, val) in zip(cols, row)}
  return res

def getComparisonPairs(data, cmp, cur):
  pairs = {};
  for key in cmp:
    ar = cmp[key]
    for val in ar:
      if val[2] >= 0:
        pairs[''+str(val[0])+","+str(val[1])] = cookCmpCascadeInfo(data[val[2]])
      else:
        if cur: pairs[''+str(val[0])+","+str(val[1])] = cookCmpCascadeInfoFromDb(val[0], cur)
        else: print ("error") # TODO throw error
  return json.dumps(pairs)

if __name__ == "__main__":
  if len(sys.argv) < 3:
    print("please provide dir. with new data (having output from anuvikar validate script), destination db path, existing added db path (optional).")
  else:
    newPath = sys.argv[1]
    dest = sys.argv[2]
    dbNew = os.path.join(newPath, "anuvikar.db")
    data = os.path.join(newPath, "anuvikar.json")
    dbOld = None
    nnOld = None
    if len(sys.argv) > 3: 
      dbOld = sys.argv[3]
      nnOld = sys.argv[3] + "_tree.pickle"
    isSuccess, msg = mergeCascadeDbs(dbNew, data, dest, dbOld, nnOld)
    if (isSuccess):
      print("New database written to ", dest)
      print("Copy this file to Csaransh src/db/dev.csaransh.db to view with Csaransh")