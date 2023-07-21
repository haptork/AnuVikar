import json
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics.cluster import adjusted_rand_score
from sklearn.metrics.cluster import adjusted_mutual_info_score
from sklearn.neighbors import NearestNeighbors
from scipy.spatial import ConvexHull
import math
import sys
import os

#import umap
import logging

from pysrc.venuFeat import getPointDefectLines, linesForCascade, makeLatticeGroups
from pysrc.savi import addFullComponentInfo

"""
Fits and applies PCA transformation to the given points.
Applies same transformation to other points given as second argument list.
"""


def findEigen(points, others):
    pca = PCA(n_components=3)
    pca.fit(points)
    var = pca.explained_variance_ratio_.tolist()
    eigen_coords = pca.transform(points)
    other_trans = []
    for x in others:
        temp = pca.transform(x)
        temp = [[round(x[0], 4), round(x[1], 4), round(x[2], 4)]
                for x in temp.tolist()]
        other_trans.append(temp)
    eigen_coords = [[round(x[0], 4), round(x[1], 4), round(x[2], 4)]
                    for x in eigen_coords.tolist()]
    var = [round(x, 4) for x in var]
    return eigen_coords, var, other_trans


"""
  Applies PCA transformation to the whole cascade coordinates and PKA point.
"""


def transformCascade(coords, point):
    all_coords = np.asarray(coords)
    point = np.asarray([point])
    if (len(coords) == 0):
        return(all_coords, point, [0, 0, 0])
    all_coords = all_coords[:, 0:3]
    trans_coords, var, trans_others = findEigen(all_coords, [point])
    trans_point = trans_others[0]
    return (trans_coords, trans_point, var)


"""
Applies DBSCAN to the vacancy and interstitial points separately
"""


def densityClustering(coords, trans_coords, thresh):
    vac_coords = []
    int_coords = []
    for i, x in enumerate(coords):
        if x[5] == 0:
            continue  # annihilated / psedo defect
        if x[3] == 0:
            vac_coords.append(trans_coords[i])
        else:
            int_coords.append(trans_coords[i])
    if len(vac_coords) == 0:
        labels_vac = []
    else:
        dbscan_vac = DBSCAN(eps=thresh, min_samples=3).fit(vac_coords)
        labels_vac = dbscan_vac.labels_
    dbscan_int = DBSCAN(eps=thresh, min_samples=3).fit(int_coords)
    labels_int = dbscan_int.labels_
    return (labels_vac, labels_int)

def addEigenAndSubcascades(data):
    for fdata in data:
        sys.stdout.write('\rto ' + fdata["infile"] + " "*10)
        sys.stdout.flush()
        # eigen coords and variance for cascade
        pka = fdata['pka']
        eigen_coords, eigen_pka, var = transformCascade(fdata['coords'], pka)
        fdata['eigen_coords'] = eigen_coords
        fdata['eigen_pka'] = eigen_pka
        fdata['eigen_var'] = var
        # eigen coords and variance for each cluster
        eigen_features = {}
        for x in fdata['clusters']:
            if len(fdata['clusters'][x]) >= 3:
                c = [fdata['coords'][y][:3] for y in fdata['clusters'][x]]
                ec, ev, _ = findEigen(c, [])
                eigen_features[x] = {'coords': ec, 'var': ev}
            else:
                eigen_features[x] = {'coords': [fdata['coords'][y][:3] for y in fdata['clusters'][x]], 'var': 1.0}
        fdata['eigen_features'] = eigen_features
        if len(fdata['coords']) == 0:
            fdata['eigen_coords'] = []
            fdata['eigen_pka'] = [round(x[0], 4)
                                  for x in fdata['eigen_pka'].tolist()]
            fdata['dclust_coords'] = {}
            fdata['dclustI_count'] = 0
            fdata['dclustV_count'] = 0
            fdata['dclust_sec_impact'] = 0
            continue
        # density clustering
        bias = 3.0
        nn4 = fdata['latticeConst'] * (3**0.5)
        labelsV, labelsI = densityClustering(
            fdata['coords'], eigen_coords, bias * nn4)
        indicesV = [i for i, x in enumerate(
            fdata['coords']) if x[3] == 0 and x[5] == 1]
        subsv = {}
        for i, x in enumerate(labelsV):
            if x == -1:
                continue
            x = str(x) + "v"
            if not x in subsv:
                subsv[x] = []
            subsv[x].append(indicesV[i])
        indicesI = [i for i, x in enumerate(
            fdata['coords']) if x[3] == 1 and x[5] == 1]
        subsi = {}
        for i, x in enumerate(labelsI):
            if x == -1:
                continue
            x = str(x) + "i"
            if not x in subsi:
                subsi[x] = []
            subsi[x].append(indicesI[i])
        # sorting subsv by length
        lenV = []
        totalV = 0
        for x in subsv:
            lenV.append([len(subsv[x]), x])
            totalV += len(subsv[x])
        lenV.sort(reverse=True)
        # criterion to filter only substantial subcascades
        dclustNamesV = [x[1] for i, x in enumerate(lenV) if (
            x[0] / (totalV/len(lenV)) > 0.55 and x[0] > 4) or (i < 2 and x[0] > 4)]
        fdata['dclust_coords'] = {}
        doLimit = len(dclustNamesV) > 1
        for x in subsv:
            if x in dclustNamesV and (not doLimit or len(subsv[x]) > 4):
                fdata['dclust_coords'][x] = subsv[x]
        fdata['dclustI_count'] = len(subsi)
        fdata['dclustV_count'] = len(fdata['dclust_coords'])
        if fdata['dclustV_count'] == 0 and fdata['n_defects'] > 0: fdata['dclustV_count'] = 1
        dclust_len = [0, 0]
        if (len(lenV) > 0):
            dclust_len[0] = lenV[0][0]
        if (len(lenV) > 1):
            dclust_len[1] = lenV[1][0]
        fdata['dclust_sec_impact'] = 0
        if (fdata['dclustV_count']>1 and dclust_len[0] > 0):
            fdata['dclust_sec_impact'] = dclust_len[1] * 100 / dclust_len[0]

# chiSqr distance criterion for cluster comparison

def dist(a, b):
    res = 0.0
    for x, y in zip(a, b):
        if (abs(x) > 1e-6):
            res += ((x - y)**2 * 1.0) / (1.0*x)
    return round(res, 4)

"""
Helper distance function for dimensionality reduction
"""
#@numba.njit()
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
#@numba.njit()
def quad1(x, y):
    l = x.shape[0]
    a = chiSqr(x, y, 0, 0, 36, 36)
    d = chiSqr(x, y, 36, 36, l, l)
    preA = chiSqr(x, y, 0, 1, 35, 36)
    postA = chiSqr(x, y, 1, 0, 36, 35)
    preD = chiSqr(x, y, 36, 37, l - 1, l)
    postD = chiSqr(x, y, 37, 36, l, l - 1)
    return (0.5 * (preA + postA) + a + 0.9 * (0.1 * (preD + postD) + d)) / (2.0 + 1.2*0.9)


"""
Distance function for dimensionality reduction
"""
#@numba.njit()
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
            feat.append(x['features'][y]['angle'] + x['features'][y]['dist'])
            tag.append((i, y))
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


"""
Adds top 5 matching clusters for each cluster in the data
"""

def addClusterCmp(data, feat, tag):
    topsize = 5
    neigh = {}
    keys = ['angle', 'dist', 'all']
    #quadAngle = quadCustom(1.0, 0.0)
    #quadDist = quadCustom(0.0, 1.0)
    #quadBoth = quad
    defaultK = topsize * 3 if topsize * 3 < len(feat) else len(feat) - 1
    neigh[keys[0]] = NearestNeighbors(defaultK)
    neigh[keys[1]] = neigh[keys[0]] 
    neigh[keys[2]] = neigh[keys[0]] 
    dists = {}
    neighbours = {}
    for key in neigh:
        if len(feat) == 0: continue
        neigh[key].fit(feat)
        dists[key], neighbours[key] = neigh[key].kneighbors()
    for index, (cascadeIndex, cid) in enumerate(tag):
        cascade = data[cascadeIndex]
        if not 'clust_cmp' in cascade:
            cascade['clust_cmp'] = {}
            cascade['clust_cmp_size'] = {}
            cascade['clust_cmp'][cid] = {}
            cascade['clust_cmp_size'][cid] = {}
        elif not cid in cascade['clust_cmp']:
            cascade['clust_cmp'][cid] = {}
            cascade['clust_cmp_size'][cid] = {}
        for key in neigh:
            if key not in dists:
              cascade['clust_cmp'][cid][key] = []
              cascade['clust_cmp_size'][cid][key] = []
            cascade['clust_cmp'][cid][key] = [(x, tag[y][0], tag[y][1]) for x, y in zip(
                dists[key][index][:topsize], neighbours[key][index][:topsize])]
            curLen = len(cascade['clusters'][cid])
            lenDiff = [(abs(curLen - len(data[tag[x][0]]['clusters'][tag[x][1]])), i)
                       for i, x in enumerate(neighbours[key][index])]
            lenDiff.sort()
            cascade['clust_cmp_size'][cid][key] = [
                (dists[key][index][x[1]], tag[neighbours[key][index][x[1]]][0], tag[neighbours[key][index][x[1]]][1]) for x in lenDiff[:topsize]]

def clusterClasses(data, feat, tag):
    classesData = [{}] # first item reserved for line-comp
    show_dim = []
    if len(feat) == 0: return False
    rndSeed = 42
    #reduced_dim = umap.UMAP(n_components=2, n_neighbors=6, min_dist=0.40,
    #                        metric=quad, random_state=rndSeed).fit_transform(feat).tolist()
    reduced_dim = TSNE(n_components=2, random_state=rndSeed).fit_transform(feat).tolist()
    preId = -1
    for (tcas, tcid), dim in zip(tag,reduced_dim):
        cascade = data[tcas]
        if cascade['clusterSizes'][tcid] < 2: 
        #if cascade['clusterSizes'][tcid] < 2 or len(cascade['clusters'][tcid]) > 3000: 
          compLabel = "7-?"
          if cascade['clusterSizes'][tcid] < 0: compLabel = "9-v"
          #else: print("-- ignoring", cascade['id'], cascade['xyzFilePath'], tcid, cascade['clusterSizes'][tcid], cascade['infile'], len(cascade['clusters'][tcid]), cascade['n_defects'])
          if not 'clusterClasses' in cascade: cascade['clusterClasses'] = {}
          if not 'savi' in cascade['clusterClasses']: cascade['clusterClasses']['savi'] = {}
          cascade['clusterClasses']['savi'][tcid] = {"morph":compLabel}
        else:
          if tcas != preId:
            preId = tcas
            triads, pairs = makeLatticeGroups(cascade)
            if 'siavenu' not in cascade:
              cascade['siavenu'] = getPointDefectLines(cascade, triads, pairs)
          #print("adding savi info for", cascade['id'], cascade['xyzFilePath'], tcid, cascade['clusterSizes'][tcid], cascade['infile'], len(cascade['clusters'][tcid]))
          isSuccess = addFullComponentInfo(cascade, tcid, triads, pairs)
          if not isSuccess:
            compLabel = "7-?"
            if not 'clusterClasses' in cascade: cascade['clusterClasses'] = {}
            if not 'savi' in cascade['clusterClasses']: cascade['clusterClasses']['savi'] = {}
            cascade['clusterClasses']['savi'][tcid] = {"morph":compLabel}
          #print(cascade['clusterClasses']['savi'][tcid]['morph'])
        cascade['clusterClasses']['savi'][tcid]['hdbpoint'] = dim
    for cascade in data:
      if 'siavenu' not in cascade:
        triads, pairs = makeLatticeGroups(cascade)
        cascade['siavenu'] = getPointDefectLines(cascade, triads, pairs)
      if 'savi' not in cascade:
        cascade['savi'] = {}
    return True

def addHull(data):
    for cascade in data:
        li = [x for x in cascade['coords'] if x[3] == 0 and x[5] == 1]
        li = [x[0:3] for x in li]
        cascade['hull_vol'] = 0.0
        cascade['hull_area'] = 0.0
        cascade['hull_density'] = 0.0
        cascade['hull_vertices'] = []
        cascade['hull_simplices'] = []
        cascade['hull_nvertices'] = 0
        cascade['hull_nsimplices'] = 0
        #cascade['hull_neigh'] = []
        if (len(li) < 4):
            continue
        try:
            hull = ConvexHull(li)
            cascade['hull_vol'] = round(hull.volume, 2)
            cascade['hull_area'] = round(hull.area, 2)
            cascade['hull_density'] = (cascade['n_defects'] / hull.volume)
            cascade['hull_vertices'] = hull.vertices.tolist()
            cascade['hull_simplices'] = hull.simplices.tolist()
            cascade['hull_nvertices'] = len(hull.vertices)
            cascade['hull_nsimplices'] = hull.nsimplex
        except:
            pass
            #cascade['hull_neigh'] = hull.neighbors.tolist()

def checkIds(cascades):
    uniqueIdCount = len(set([cascade['id'] for cascade in cascades]))
    if uniqueIdCount == len(cascades): return
    print("Duplicate Ids found. Resetting all ids...")
    logging.warning("Duplicate Ids found. Resetting all ids...")
    for i in range(len(cascades)):
        cascades[i]['id'] = i + 1

# top level user functions
# -------------------------

def getBasicLines(cascade, cid):
    return linesForCascade(cascade, cid)

def validateForCdb(cascades, isInit=True, isAddClusterComparison=False, isAddClassification=True):
    if isInit:
      checkIds(cascades)
      print("Adding coordinates in eigen dimensions and convex hulls for cascades...")
      addEigenAndSubcascades(cascades)
      addHull(cascades)
      print('finished.')
    if isAddClusterComparison or isAddClassification:
      feat, tag = clusterClassData(cascades)
    if isAddClusterComparison:
      print("Adding cluster comparison...")
      addClusterCmp(cascades, feat, tag)
      print('finished.')
      sys.stdout.flush()
    if isAddClassification:
      print("Defect morphology identification & classification...")
      isSuccess = clusterClasses(cascades, feat, tag)
      if (isSuccess): 
          print("finished.")
      else: 
          logging.error("Finished with error!")
          print("finished with error!")
      sys.stdout.flush()
    return cascades