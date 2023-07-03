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
pathToAnuvikar = "."
sys.path.append(pathToAnuvikar)
from pysrc.anuvikar_cdb_helper import getDefaultConfig, writeResultsToJSON
buildDirs = []
libPaths = []

from pysrc.anuvikar_cdb_helper import processXyzFilesInDirGivenMetaFile, _unzipFile
from pysrc.anuvikar_ml import validateForCdb
import logging

def stageEkaCpp(metaFilePath, xyzDir, config, prefixLen):
    isSuccess, cascades = processXyzFilesInDirGivenMetaFile(metaFilePath, xyzDir, config, isRemoveDirectoryFromNames=prefixLen)
    msg = ''
    if not isSuccess: msg = cascades
    return [isSuccess, msg, cascades]

def stageDwiMl(cascades):
  return validateForCdb(cascades)

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

def validateArchive(srcDir, extractionDir, metaFiles, overwriteJson, overwriteDb):
  config = getDefaultConfig("warning", "error") # check bottom cell for various keys / options to configure
  config['logFilePath'] = os.path.join(extractionDir, "avi_cpp.log")
  config['outputJSONFilePath'] = os.path.join(extractionDir, "anuvikar.json")
  config['outputDbPath'] = os.path.join(extractionDir, "anuvikar.db")
  config['anuvikarLib'] = libPath
  cascades = []
  if not(os.path.exists(config['outputJSONFilePath'])) or overwriteJson:
    for metaFilePath in metaFiles:
      print("parsing and finding archive for", metaFilePath)
      meta = xmlFileToDict(metaFilePath)
      archiveName = meta['cdbml']['cdbrecord']['data']['archive_name']
      archivePath = os.path.join(srcDir, archiveName)
      #archiveName = os.path.basename(archivePath)
      if not os.path.exists(archivePath):
        logging.error("Archive missing from the source directory: " + archivePath)
        print("Archive missing from the source directory: " + archivePath + " for " + metaFilePath)
        continue
      xyzDir = _unzipFile(extractionDir, archivePath, archiveName)
      if len(xyzDir) == 0: 
        logging.error("Error in unzipping the file: " + archivePath + " for " + metaFilePath)
        print("Error in unzipping the file: " + archivePath)
        continue
      print("processing: ", archiveName)
      isSuccess, msg, curCascades = stageEkaCpp(metaFilePath, xyzDir, config, len(extractionDir))
      if not isSuccess: 
        logging.error("Error in cpp processing files in : " + xyzDir + " for " + metaFilePath + ": " + msg)
        print("Error in cpp processing of file: " + metaFilePath + ": " + msg)
      else:
        cascades += curCascades
        #return [isSuccess, msg]
    writeResultsToJSON(cascades, config)
    print("Analysis results written to JSON: " + config['outputJSONFilePath'])
  else:
    print("Using existing processed json file. To carry out fresh analysis delete/move anuvikar.json from extraction dir.")
    logging.warning("Processed json file 'anuvikar.json' already exists. Using it for further processing.")
    f = open(config['outputJSONFilePath'], 'r')
    cascades = json.load(f)
    f.close()
  cascades = stageDwiMl(cascades)
  config['outputJSONFilePath'] =  os.path.join(extractionDir, "anuvikar_ml.json")
  writeResultsToJSON(cascades, config)
  print("Json with Ml results written: " + config['outputJSONFilePath'])
  if not(os.path.exists(config['outputDbPath'])) or overwriteDb:
    if os.path.exists(config['outputDbPath']):
      os.remove(config['outputDbPath'])
    writeMlResultsToSqliteDb(cascades, config)
    print("Analysis results written to db: " + config['outputDbPath'])
  print("Please go through logs written to: ", extractionDir, ": avi_cpp.log and avi_py.log.")
  print("Run add script to generate db that can be viewed with Csaransh.")
  summarizeLog(config)
  # saved cascades.json, cascades.db, log-summary.md, log.txt
  return [True, ""]

def writeResultsToJSON(res, config):
    f = open(config['outputJSONFilePath'], "w")
    json.dump(res, f)
    f.close()

def  summarizeLog(config):
    pass

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("please provide srcdir, extractiondir, metaFiles")
        exit(1)
    srcDir = sys.argv[1]
    extractionDir = sys.argv[2]
    if not(os.path.exists(srcDir)):
      print("Source directory path is not accessible.")
      exit(1)
    if not(os.path.exists(extractionDir)):
      print("Extraction directory path is not accessible.")
      exit(1)
    #overwriteJson = True
    overwriteJson = False
    overwriteDb = True
    buildDirs.append(os.path.join(pathToAnuvikar, "_build"))
    buildDirs.append(os.path.join(pathToAnuvikar, "build"))
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
    #overwriteJson = (sys.argv[3] == "1")
    #overwriteDb = (sys.argv[4] == "1")
    metaFiles = [x for x in sys.argv[3:]]
    logging.basicConfig(filename=os.path.join(extractionDir, "avi_py.log"), level=logging.WARNING, format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.warning("Atho Anuvikar validation")
    validateArchive(srcDir, extractionDir, metaFiles, overwriteJson, overwriteDb)
    logging.warning("Iti")