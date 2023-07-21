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
import os
import json
import sqlite3
import pandas as pd
import seaborn as sns
from pandas import DataFrame
import xmltodict
import numba
import umap
from sklearn.neighbors import NearestNeighbors
import shutil
import numpy as np
import pickle

import logging

dbToJsonMap = {
  "id": "id",
  "cascadeid": "id",
  'ncell': 'ncell',
  'energy': 'energy',
  'latticeconst': 'latticeConst',
  'structure': 'structure',
  'temperature': 'temperature',
  'simulationtime': 'simulationTime',
  'infile': 'infile',
  'xyzfilepath': "xyzFilePath",
  'substrate': 'substrate',
  'simulationcode': 'simulationCode',
  'potentialused': 'potentialUsed',
  'author': 'author',
  'es': 'es',
  'tags': 'tags',
  'ndefects': 'n_defects',
  "nclusters": 'n_clusters',
  "maxclustersize": 'max_cluster_size',
  "maxclustersizei": 'max_cluster_size_I',
  "maxclustersizev": 'max_cluster_size_V',
  "incluster": 'in_cluster',
  "inclusteri": 'in_cluster_I',
  "inclusterv": 'in_cluster_V',
  "ndclusti": 'dclustI_count',
  "dclustsecimpact": 'dclust_sec_impact',
  "hullvol": 'hull_vol',
  "hulldensity": 'hull_density',
  'coords': 'coords',
  'eigencoords': 'eigen_coords',
  'dclustcoords': 'dclust_coords',
  'codefects': 'coDefects',
  'clusters': 'clusters',
  'clusterclasses': 'clusterClasses',
  'siavenu': 'siavenu',
  'simboxfoc': 'pka',
  'boxsize': 'boxSize'
};

dbToJsonMap2 = (
  ("id", "id"),
  ("cascadeid", "id"),
  ('ncell', 'ncell'),
  ('energy', 'energy'),
  ('latticeconst', 'latticeConst'),
  ('structure', 'structure'),
  ('temperature', 'temperature'),
  ('simulationtime', 'simulationTime'),
  ('infile', 'infile'),
  ('xyzfilepath', "xyzFilePath"),
  ('substrate', 'substrate'),
  ('simulationcode', 'simulationCode'),
  ('potentialused', 'potentialUsed'),
  ('author', 'author'),
  ('es', 'es'),
  ('tags', 'tags'),
  ('ndefects', 'n_defects'),
  ("nclusters", 'n_clusters'),
  ("maxclustersize", 'max_cluster_size'),
  ("maxclustersizei", 'max_cluster_size_I'),
  ("maxclustersizev", 'max_cluster_size_V'),
  ("incluster", 'in_cluster'),
  ("inclusteri", 'in_cluster_I'),
  ("inclusterv", 'in_cluster_V'),
  ("ndclustv", 'dclustV_count'),
  ("dclustsecimpact", 'dclust_sec_impact'),
  ("hullvol", 'hull_vol'),
  ("hulldensity", 'hull_density'),
  ('coords', 'coords'),
  ('eigencoords', 'eigen_coords'),
  ('dclustcoords', 'dclust_coords'),
  ('codefects', 'coDefects'),
  ('clusters', 'clusters'),
  ('clusterclasses', 'clusterClasses'),
  ('siavenu', 'siavenu'),
  ('simboxfoc', 'pka'),
  ('boxsize', 'boxSize'),
);

dbTypes = {
  "id": "string",
  "cascadeid": "string",
  'ncell': 'integer',
  'energy': 'integer',
  'boxsize': 'real',
  'latticeconst': 'real',
  'structure': 'string',
  'temperature': 'real',
  'simulationtime': 'integer',
  'infile': 'string',
  'xyzfilepath': "string",
  'substrate': 'string',
  'simulationcode': 'string',
  'potentialused': 'string',
  'author': 'string',
  'es': 'integer',
  'tags': 'text',
  'ndefects': 'integer',
  "nclusters": 'integer',
  "maxclustersize": 'integer',
  "maxclustersizei": 'integer',
  "maxclustersizev": 'integer',
  "incluster": 'integer',
  "inclusteri": 'integer',
  "inclusterv": 'integer',
  "ndclustv": 'integer',
  "dclustsecimpact": 'real',
  "hullvol": 'real',
  "hulldensity": 'real',
  'viewfields': 'text',
  'created_at': 'text'
};

def cookCascadesDbTuple(cascades):
  rows = []
  for cascade in cascades:
    row = []
    for val in dbToJsonMap2:
      #if val[0] == "id": continue
      if val[0] == "coords": break
      res = cascade[val[1]]
      #if (dbTypes[val[0]] == 'string' and type(res) != str): res = str(res)
      row.append(res)
    row.append(json.dumps({
      'coords': cascade['coords'],
      'savi': cascade['savi'] if 'savi' in cascade else {},
      'clusters': cascade['clusters'],
      'clustersizes': cascade['clusterSizes'],
      'clusterclasses': cascade['clusterClasses'] if 'clusterClasses' in cascade else {"savi":{}}, # TODO insert anyway
      'eigencoords': cascade['eigen_coords'],
      'dclustcoords': cascade['dclust_coords'],
      'siavenu': cascade['siavenu'] if 'siavenu' in cascade else [],
      'simboxfoc': cascade['pka'],
      'boxsize': cascade['boxSize']
    }))
    rows.append(tuple(row))
  columns = []
  for val in dbToJsonMap2:
    #if val[0] == "id": continue
    if val[0] == "coords": 
      columns.append('viewfields')
      break
    columns.append(val[0])
  return (rows, columns)

def addCascadesTable(cascades, cur):
  cur.execute('''create table cascades
               (id text Primary key, cascadeid text unique, ncell integer, energy integer, latticeconst real not null,
                structure string, temperature real, simulationtime real, infile string, xyzfilepath string not null,
                substrate sring, simulationcode string, potentialused string, author string,
                es integer, tags text, ndefects integer, nclusters integer,
                maxclustersize integer, maxclustersizei integer, maxclustersizev integer,
                incluster integer, inclusteri integer, inclusterv integer, ndclustv integer,
                dclustsecimpact integer, hullvol real, hulldensity real, viewfields text,
                created_at DATETIME DEFAULT CURRENT_TIMESTAMP)''')
  cascadesTuple, columns = cookCascadesDbTuple(cascades)
  if len(cascadesTuple) > 0:
    li = ["?"]*len(cascadesTuple[0])
    qStr = "("+ ",".join(li) + ")"
    #print(cascadesTuple[0])
    cStr = "(" + ",".join(columns) + ")"
    #print(cStr)
    #print(cascadesTuple[0])
    #cur.execute("INSERT INTO cascades" + cStr + " VALUES " + qStr, cascadesTuple[0]);
    cur.executemany("INSERT INTO cascades" + cStr + " VALUES " + qStr, cascadesTuple);

def getCoordType (row, cid):
  cid = str(cid);
  if len(cid) == 0 or cid not in row['savi']: return -1;
  if 'savi' in row and cid in row['savi'] and 'venu' in row['savi'][cid]: return 1;
  return 0;

def getClusterCoord(row, cid):
  c = [[],[],[]];
  if (len(cid) == 0 or not row or cid not in row['eigen_features']): return c
  for x in row['eigen_features'][cid]['coords']:
    c[0].append(x[0]);
    c[1].append(x[1]);
    c[2].append(x[2]);
  return c

def getClusterLineCoord(row, cid):
  lines = [];
  linesT = [];
  pointsI = [[], [], [], []];
  pointsV = [[], [], [], []];
  cid = str(cid)
  if cid:
    for x in row['savi'][cid]['venu']['linesT']:
      c = [[],[],[],[], [-1.0, -1.0]]
      c[4] = x['orient']
      for y in x['main']:
        curCoord = row['coords'][y]
        c[0].append(curCoord[0]);
        c[1].append(curCoord[1]);
        c[2].append(curCoord[2]);
        c[3].append(str(x['orient']));
      linesT.append(c);
    for x in row['savi'][cid]['venu']['lines']:
      c = [[],[],[],[], [-1.0, -1.0]];
      c2 = [[],[],[],[], [-1.0, -1.0]];
      c[4] = x['orient']
      c2[4] = x['orient']
      for y in x['main']:
        curCoord = row['coords'][y]
        c[0].append(curCoord[0]);
        c[1].append(curCoord[1]);
        c[2].append(curCoord[2]);
        c[3].append(str(x['orient']));
      for y in x['sub']:
        curCoord = row['coords'][y]
        c2[0].append(curCoord[0]);
        c2[1].append(curCoord[1]);
        c2[2].append(curCoord[2]);
        c2[3].append(str(x['orient']));
      lines.append({'main':c, 'sub':c2});
    for x in row['savi'][cid]['venu']['pointsI']:
      curCoord = row['coords'][x]
      pointsI[0].append(curCoord[0]);
      pointsI[1].append(curCoord[1]);
      pointsI[2].append(curCoord[2]);
      pointsI[3].append(x);
    for x in row['savi'][cid]['venu']['pointsV']:
      curCoord = row['coords'][x]
      pointsV[0].append(curCoord[0]);
      pointsV[1].append(curCoord[1]);
      pointsV[2].append(curCoord[2]);
      pointsV[3].append(x);
  return {"lines":lines, "linesT":linesT, "pointsI":pointsI, "pointsV":pointsV};

def cookClustersDbTuple(cascades):
  rows = []
  for cascade in cascades:
    for clusterName in cascade['clusters']:
      row = [cascade['id'], clusterName]
      curClusterClass = cascade['clusterClasses']['savi'][clusterName]
      row.append(cascade['clusterSizes'][clusterName])
      row.append(curClusterClass['morph'])
      coordType = getCoordType(cascade, clusterName);
      coords = getClusterLineCoord(cascade, clusterName) if (coordType == 1) else getClusterCoord(cascade, clusterName);
      row.append(coordType)
      row.append(json.dumps(coords))
      # coords and coordtype
      row.append(curClusterClass['hdbpoint'][0])
      row.append(curClusterClass['hdbpoint'][1])
      rows.append(tuple(row))
  columns = ["cascadeid", "name", "size", "savimorph", "coordtype", "coords", "hdbx", "hdby"]
  return (rows, columns)

def addClustersTable(cascades, cur):
  cur.execute('''create table clusters
               (id integer Primary Key, 
                cascadeid text not null, name integer, savimorph text, size integer, coordtype integer,
                coords text, hdbx real, hdby real, cmp text, cmpsize text, cmppairs text, morphdesc text,
                properties text, created_at DATETIME DEFAUlT CURRENT_TIMESTAMP, 
                foreign key (cascadeid) references cascades (id), unique(cascadeid, name))''')
  clusterTuples, columns = cookClustersDbTuple(cascades)
  if len(clusterTuples) > 0:
    li = ["?"]*len(clusterTuples[0])
    qStr = "("+ ",".join(li) + ")"
    cStr = "(" + ",".join(columns) + ")"
    cur.executemany("INSERT INTO clusters" + cStr + " VALUES " + qStr, clusterTuples);


def writeMlResultsToSqliteDb(cascades, config, isOverwrite=True):
  con = None
  try: 
    con = sqlite3.connect(config['outputDbPath'])
  except sqlite3.Error as e:
    print(e)
    logging.error("Error in saving db file: " + str(e))
    if con: con.close()
    return False
  cur = con.cursor()
  addCascadesTable(cascades, cur)
  addClustersTable(cascades, cur)
  cur.execute("CREATE UNIQUE INDEX 'cascades_xyzfilepath_unique' on 'cascades' ('xyzfilepath')")
  cur.execute("CREATE UNIQUE INDEX 'clusters_cascadeid_name_unique' on 'clusters' ('cascadeid', 'name')")
  con.commit()
  con.close()
  return True

def xmlFileToDict(fname):
    f = open(fname, 'r')
    xmlStr = f.read()
    di = xmltodict.parse(xmlStr)
    return di

def writeResultsToJSON(res, config):
    f = open(config['outputJSONFilePath'], "w")
    json.dump(res, f)
    f.close()

def  summarizeLog(config):
    pass

## == Add functions


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