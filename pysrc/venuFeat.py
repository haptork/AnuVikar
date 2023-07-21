#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import math
import seaborn as sns
import matplotlib.pyplot as plt
from .geo import Point, Line, use_degrees, Plane, abs2
import json
use_degrees()
angular_unit = 180.0/math.pi
debug = False
import numpy as np
from sklearn.neighbors import NearestNeighbors
from collections import Counter

# In[1098]:

def getLinesData(data):
    finalCluster = []
    attrs = []
    tags = []
    for i, cascade in enumerate(data):
        sorted_clusters = sorted(cascade['clusters'].keys())
        #print(i, cascade['id'])
        triads, pairs = makeLatticeGroups(cascade)
        for cid in sorted_clusters:
            if cascade['clusterSizes'][cid] < 0: continue
            lines, pointLineMap, freeIs = getTriadLines(cascade, cid, triads)
            linesT, pointLineMapT = getThreshLines(cascade, cid, len(lines), pairs)
            extraVs = getExtraFreeVs(cascade, cid, pointLineMap, pointLineMapT)
            freeVs, extraVsMap = addCollinearToTriads(cascade, cid, lines, pointLineMap, linesT, pointLineMapT, extraVs)
            attrs.append(cookAttrs(cascade, cid, lines + linesT, freeIs, extraVs))
            allVs = list(freeVs) + list(extraVs)
            finalCluster.append({'lines':lines, 'linesT':linesT, 'pointsI':freeIs, 'pointsV':allVs})
            tags.append((i, cid))
    return finalCluster, tags, attrs

def getLinesDataOnly(data):
    finalCluster = []
    tags = []
    for i, cascade in enumerate(data):
        sorted_clusters = sorted(cascade['clusters'].keys())
        #print(i, cascade['id'])
        triads, pairs = makeLatticeGroups(cascade)
        for cid in sorted_clusters:
            if cascade['clusterSizes'][cid] < 0: continue
            lines, pointLineMap, freeIs = getTriadLines(cascade, cid, triads)
            linesT, pointLineMapT = getThreshLines(cascade, cid, len(lines), pairs)
            extraVs = getExtraFreeVs(cascade, cid, pointLineMap, pointLineMapT)
            freeVs, extraVsMap = addCollinearToTriads(cascade, cid, lines, pointLineMap, linesT, pointLineMapT, extraVs)
            allVs = list(freeVs) + list(extraVs)
            finalCluster.append({'lines':lines, 'linesT':linesT, 'pointsI':freeIs, 'pointsV':allVs})   
            tags.append((i, cid))
    return finalCluster, tags

def linesForCascade(cascade, cid):
    finalCluster = []
    attrs = []
    tags = []
    if cascade['clusterSizes'][cid] < 0: return None
    triads, pairs = makeLatticeGroups(cascade)
    lines, pointLineMap, freeIs = getTriadLines(cascade, cid, triads)
    for line in lines:
      line['forceAlign'] = forceAlign(line['eq'])
    #linesT, pointLineMapT = getThreshLines(cascade, cid, len(lines), pairs)
    #extraVs = getExtraFreeVs(cascade, cid, pointLineMap, pointLineMapT)
    #freeVs, extraVsMap = addCollinearToTriads(cascade, cid, lines, pointLineMap, linesT, pointLineMapT, extraVs)
    #attrs.append(cookAttrs(cascade, cid, lines + linesT, freeIs, extraVs))
    #allVs = list(freeVs) + list(extraVs)
    return lines#, triads, pairs

def lineFeatsForCluster(cascade, cid, triads, pairs):
    if cascade['clusterSizes'][cid] < 0: return None
    #triads, pairs = makeLatticeGroups(cascade)
    lines, pointLineMap, freeIs = getTriadLines(cascade, cid, triads)
    for line in lines:
      line['forceAlign'] = forceAlign(line['eq'])
    linesT, pointLineMapT = getThreshLines(cascade, cid, len(lines), pairs)
    extraVs = getExtraFreeVs(cascade, cid, pointLineMap, pointLineMapT)
    freeVs, extraVsMap = addCollinearToTriads(cascade, cid, lines, pointLineMap, linesT, pointLineMapT, extraVs)
    allVs = list(freeVs) + list(extraVs)
    return {'allLines':lines+linesT, 'lines':lines, 'linesT':linesT, 'pointsI':freeIs, 'pointsV':allVs}

def getAllAttrs(finalCluster, tags, cascades):
    attrs = []    
    for x, (tcas, tcid) in zip(finalCluster, tags):
        #print(tcas, tcid)
        attrs.append(cookAttrs(cascades[tcas], tcid, x['lines'] + x['linesT'], x['pointsI'], x['pointsV']))
    return attrs

def getAllHistAttrs(attrs, tags, cascades):
    maxAttrs = cookMaxAttrs(attrs)
    histAttrs = []
    for x, (tcas, _) in zip(attrs, tags):
        nnDistBins = cookNNBins(cascades[tcas]['latticeConst']) ## Make it take max for last bin edge
        histAttrs.append(cookHistAttrs(x, maxAttrs, nnDistBins))
    return histAttrs

def addLinesData(data, finalClusters, tags):
    for lineC, (tcas, tcid) in zip(finalClusters, tags):
        di = {'lines':[], 'linesT':[], 'pointsI':lineC['pointsI'], 'pointsV': lineC['pointsV']}
        #print(lineC)
        for line in lineC['lines']:
            if 'parent' in line or 'del' in line: continue
            sub = []
            if line['subLine']: sub = [x[1] for x in line['subLine']['points']]
            di['lines'].append({'main': [x[1] for x in line['mainPoints']], 'sub': sub, 'orient': line['forceAlign']['dir']})
        for line in lineC['linesT']:
            if 'parent' in line or 'del' in line: continue
            di['linesT'].append({'main': [x[1] for x in line['mainPoints']], 'orient': line['forceAlign']['dir']})
        data[tcas]['features'][tcid]['lines'] = di

def getLineFeat(lineC, clusterLi):
    di = {'lines':[], 'linesT':[], 'pointsI':[clusterLi[x] for x in lineC['pointsI']], 'pointsV': [clusterLi[x] for x in lineC['pointsV']]}
    for line in lineC['lines']:
        if 'parent' in line or 'del' in line: continue
        sub = []
        if line['subLine']: sub = [clusterLi[x[1]] for x in line['subLine']['points']]
        di['lines'].append({'main': [clusterLi[x[1]] for x in line['mainPoints']], 'sub': sub, 'orient': line['forceAlign']['dir']})
    for line in lineC['linesT']:
        if 'parent' in line or 'del' in line: continue
        di['linesT'].append({'main': [clusterLi[x[1]] for x in line['mainPoints']], 'orient': line['forceAlign']['dir']})
    return di

# In[1235]:



stdPlanes = [Plane(Point(0.0,0.0,0.0), Point(0.0, 2.0, 5.0), Point(0.0, 4.0, 9.0)),
             Plane(Point(0.0,0.0,0.0), Point(1.0, 0.0, 4.0), Point(4.0, 0.0, 5.0)),
             Plane(Point(0.0,0.0,0.0), Point(1.0, 2.0, 0.0), Point(3.0, 4.0, 0.0))]
def getSymmetricPlaneAngle(line):
    res = []
    for x in stdPlanes:
        res.append(x.angle_to(line))
    return res

def getSymmetricPlaneAngleOld(line):
    res = [0.0, 0.0, 0.0]
    for i, x in enumerate(line.t):
        j = i + 1 if i < 2 else 0
        angle = 0.0
        if abs(line.t[j]) > 0.0: angle = abs(angular_unit*math.atan(line.t[i] / line.t[j]))
        if angle > 45: angle = 90 - angle
        res[i] = angle
    return res

def forceAlign(line):
    angleTol = 18
    milErrTol = 0.20
    milErrTotTol = 0.35
    subP = (line.r - line.r2)
    divP = max(abs(subP))
    if divP == 0.0:
      return {"type": (0,0,0), 'type1': (0,0,0), 'dir': (0.0, 0.0), 'eq':line, 'angles': (0.0,0.0), 'err': (False, 0.0, 0.0, 0.0)}
    milBasic = subP/divP
    milBasicAbs = abs(milBasic)
    # method 1
    milErrSum = 0.0
    milErrMax = 0.0    
    milErrList = []
    milMult = 1
    perfect = False
    for milMult in range(1, 3):
        curMilAbs = milBasicAbs * milMult      
        milErrList = abs(curMilAbs - np.rint(curMilAbs))
        milErrSum = sum(milErrList)
        milErrMax = max(milErrList)
        #print(milMult, milErrList)        
        if milErrMax < (milErrTol / milMult) and milErrSum < (milErrTotTol / milMult):
            perfect = True
            break
    #print(milErrList)
    if not perfect:
        milErrOne = sum(abs(milBasicAbs - np.rint(milBasicAbs)))
        if milErrOne < milErrSum * 2:
            milErr = milErrOne
            milMult = 1
    #print(milErr, milErrList)
    milFamOne = tuple(sorted([abs(int(x)) for x in np.rint(milBasic)]))
    milDir = np.rint(milBasic * milMult)
    milFam = tuple(sorted([abs(int(x)) for x in milDir]))
    p1 = Point(line.r)
    pMil = Point(p1.r - milDir)
    if abs2(p1.r - pMil.r) > 0.0001: 
      lineN = Line(p1, pMil)
      angles = [0.0, 0.0, 0.0]
      angErr = line.angle_to(lineN)
    else:
      #print("Error in force align")
      lineN = line
      angles = 0.0
      angErr = 0.0
    rawDir = tuple(map(lambda x: round(x,2), sorted(abs(milBasic))[:-1]))
    return {"type": milFam, 'type1': milFamOne, 'dir': rawDir, 'eq':lineN, 'angles': angles, 'err': (perfect, milErrSum, milErrMax, angErr)}

# In[1104]:

# In[1106]:


"""
- in all neighbouring points, check if in a line and add line only if angle is less
- 
"""

def lineDist(a, b):
    angleBetween = a['eq'].angle_to(b['eq'])
    return (angleBetween, a['eq'].distance_to(Point(b['eq'].r)), a['eq'].distance_to(Point(b['eq'].r2)))

def addLine(line, otherLine, dist):
    if 'children' not in line: 
        line['children'] = {'lines':[]}
    line['children']['lines'].append((dist, otherLine['id']))
    if 'parent' not in otherLine:
        otherLine['parent'] = []
    otherLine['parent'].append((dist, line['id']))

def addForceAlignAndSort(cascade, cid, lines, linesT):
    sortPointsHelper = lambda a: sum(cascade['eigen_features'][cid]['coords'][a[1]])
    for line in lines:
        if 'parent' in line: continue
        line['mainPoints'].sort(key=sortPointsHelper)        
        if line['subLine']:
            line['forceAlign'] = forceAlign(line['subLine']['eq'])
            line['subLine']['forceAlign'] = line['forceAlign']
            line['subLine']['points'].sort(key=sortPointsHelper)
        else:
            line['forceAlign'] = forceAlign(line['eq'])
    for line in linesT:
        if 'parent' in line: continue
        line['mainPoints'].sort(key=sortPointsHelper)        
        line['forceAlign'] = forceAlign(line['eq'])
    # TODO what if a line's neighbours are added while that line is not added?
    # possibly add Vs that don't add to subline to freeVs and points to the line?


"""
TODO
Finds if there is a free vacancy,
if other prob-vacancies around it and that vacancy forms a line that is either at an angle > eps or
at a distance if yes creates a sub line.
"""

def findSubLine(cascade, cid, line, vacsOrig):
    # TODO, continue with probP vacs with dist > threshold are more than 2
    #distTol = 0.20
    vacs = list(reversed(sorted(vacsOrig)))
    #if vacs[1][0] < distTol: return None
    mainPIndex = vacs[0][1]
    #mainP = cascade['eigen_features'][cid]['coords'][mainPIndex]
    mainPReal = cascade['coords'][cascade['clusters'][cid][mainPIndex]][:3]
    otherP = vacs[1][1]
    otherPReal = cascade['coords'][cascade['clusters'][cid][otherP]][:3]
    err = vacs[0][0]#line['vacErr'][1]
    err += vacs[1][0]
    if abs2(np.array(mainPReal) - np.array(otherPReal)) < 0.00001: return None
    # print("Error in def. subline")
    # print(mainPReal, otherPReal)
    subLine = Line(Point(mainPReal), Point(otherPReal))
    subLinePs = list()    
    subLinePs.append((0.0, mainPIndex))
    subLinePs.append((0.0, otherP))
    #if line['id'] == 11: print (vacs, subLinePs)
    for mainDist, index in vacs[2:]:
        curPoint = cascade['coords'][cascade['clusters'][cid][index]][:3]
        curDist = subLine.distance_to(Point(curPoint))
        #if line['id'] == 13: print (index, mainDist, curDist)
        if (curDist < mainDist):
            subLinePs.append((curDist, index))
            err += (mainDist - curDist)
    angles = getSymmetricPlaneAngle(subLine)
    angleErr = line['eq'].angle_to(subLine)
    if angleErr > 15 and len(subLinePs) <= 2:
        return None
    #if line['id'] == 13: print(subLinePs)
    return {'eq': subLine, 'points': subLinePs, 'angles':angles, 'err': [angleErr, err]}
# %%

def addCollinearToThresh(cascade, cid, lines, coordIndexMap, linesT, coordIndexMapT, nn, extraVs, extraVsMap):
    angleTol = 18
    distTol = [0.8, cascade['latticeConst']]
    # different criterion for there are only vac and interstitial. intersitial
    # can be slightly shifted in whatever angle so we can not consider that two
    # t-line as one line based on distance from the line subtended
    for line in linesT:
        if 'parent' in line: continue
        pointsToLookNeighs = [cascade['eigen_features'][cid]['coords'][x] for x in line['points']]
        IgnoreNeighs = set(line['points'])
        while pointsToLookNeighs:
            dist, neigh = nn.radius_neighbors(pointsToLookNeighs)
            curNeigh = set()
            for x in neigh:
                for y in x:
                    if y in IgnoreNeighs: continue
                    if y in coordIndexMap: continue
                    curNeigh.add(y)
                    if y in coordIndexMapT: 
                        IgnoreNeighs.update(linesT[coordIndexMapT[y] - len(lines)]['points'])
                    else:
                        IgnoreNeighs.add(y)
            nextPointsToLook = []
            #if (line['id'] == 10): print(curNeigh)
            for neighIndex in curNeigh:
                neighIndex = int(neighIndex)
                if neighIndex in coordIndexMapT:
                    otherLine = linesT[coordIndexMapT[neighIndex]- len(lines)]
                    #print(neighIndex, line['id'])
                    #dist = lineDist(line, otherLine)
                    angleBetween = line['eq'].angle_to(otherLine['eq'])
                    dist = line['eq'].distance_to(Point(otherLine['eq'].r))
                    vacDist = Point(line['eq'].r).distance_to(Point(otherLine['eq'].r))
                    if angleBetween < angleTol and dist < distTol[0] and vacDist < distTol[1]:
                        addLine(line, otherLine, vacDist)
                        nextPointsToLook.append(cascade['eigen_features'][cid]['coords'][neighIndex]) # TODO check index to add or something else
                        for y in otherLine['points']:
                            coordIndexMapT[y] = line['id']
                            coordIndexMapT.pop(y)# = line['id']                           
                        # add to error                                              
                elif neighIndex in extraVs:
                    coordNeighIndex = cascade['coords'][cascade['clusters'][cid][neighIndex]][:3]
                    dist = line['eq'].distance_to(Point(coordNeighIndex))
                    if neighIndex not in extraVsMap or extraVsMap[neighIndex][0] > dist:
                        if neighIndex in extraVsMap: 
                            preLineIndex = extraVsMap[neighIndex][1]
                            if preLineIndex < len(lines):
                                removeExtraV(lines[preLineIndex], neighIndex)
                            else:
                                removeExtraV(linesT[preLineIndex - len(lines)], neighIndex)
                        addExtraV(line, neighIndex, dist)
                        extraVsMap[neighIndex] = (dist, line['id'])
            pointsToLookNeighs = nextPointsToLook
        line['subLine'] = None
        curLinePs = list()
        curLinePs.append((0.0, line['points'][0]))
        curLinePs.append((0.0, line['points'][1]))                          
        if 'children' in line:
            for lineChIndex in line['children']['lines']:
                lineCh = linesT[lineChIndex[1] - len(lines)]
                for x in lineCh['points']:
                    curLinePs.append((lineChIndex[0], x))
        line['freeVs'] = list()
        line['mainPoints'] = curLinePs
        #line['mainPoints'].sort(key=sortPointsHelper)        
        #line['forceAlign'] = forceAlign(line['eq'])
    # TODO what if a line's neighbours are added while that line is not added?
    # possibly add Vs that don't add to subline to freeVs and points to the line?

def addCollinearToTriadsAgain(cascade, cid, lines, coordIndexMap, linesT, coordIndexMapT):
    latConst = cascade['latticeConst']
    distTolErr = 0.10
    distTol = [0.30 + distTolErr, 0.50 + distTolErr, 0.70 + (distTolErr * 2)]#latConst * 0.3
    nn = NearestNeighbors(n_neighbors=5, radius=latConst * 2.0)
    tIs = []
    tIsMap = {}
    for line in linesT:
        if 'parent' in line: continue
        tIsMap[len(tIs)] = line['points'][1]
        tIs.append(cascade['eigen_features'][cid]['coords'][line['points'][1]])
    if len(tIs) == 0: return
    nn.fit(tIs)
    allLines = lines + linesT
    for line in allLines:
        if 'parent' in line or len(line['mainPoints']) < 3: continue
        if len(line['mainPoints']) <= 3 and not line['subLine']: continue
        pointsToLookNeighs = [cascade['eigen_features'][cid]['coords'][x[1]] for x in line['mainPoints']]
        IgnoreNeighs = set([x[1] for x in line['mainPoints']])
        #print(line['id'], IgnoreNeighs)
        line['oldEq'] = line['eq']
        line['eq'] = Line([Point(cascade['coords'][cascade['clusters'][cid][x[1]]][:3]) for x in line['mainPoints']])
        while pointsToLookNeighs:
            dist, neigh = nn.radius_neighbors(pointsToLookNeighs)
            curNeigh = set()
            for x in neigh:
                for y in x:
                    #print(y, tIsMap[int(y)])
                    whatever = tIsMap[int(y)]
                    if whatever in IgnoreNeighs: continue
                    curNeigh.add(whatever)
                    IgnoreNeighs.add(whatever)
            nextPointsToLook = []
            #if (line['id'] == 10): (curNeigh)
            for neighIndex in curNeigh:
                if neighIndex in coordIndexMapT:
                    otherLine = linesT[coordIndexMapT[neighIndex]- len(lines)]
                    #(neighIndex, line['id'])
                    interstitialDist = line['eq'].distance_to(Point(otherLine['eq'].r2))
                    vacDist = line['subLine']['eq'].distance_to(Point(otherLine['eq'].r)) if line['subLine'] else line['eq'].distance_to(Point(otherLine['eq'].r))
                    if interstitialDist < distTol[1] and vacDist < distTol[1]: # interstitial inline
                        nextPointsToLook.append(cascade['eigen_features'][cid]['coords'][neighIndex])
                        #print(line['id'], otherLine['id'])
                        addLine(line, otherLine, interstitialDist)
                        cMap = coordIndexMap if line['id'] < len(lines) else coordIndexMapT
                        for y in otherLine['points']:
                            cMap[y] = line['id']
                            if cMap == coordIndexMap: coordIndexMapT.pop(y)
                        line['mainPoints'].append((interstitialDist, otherLine['points'][1]))
                        if line['subLine']:
                            line['subLine']['points'].append((vacDist, otherLine['points'][0]))
                        else:
                            line['mainPoints'].append((vacDist, otherLine['points'][0]))
                            # add to error                                              
            pointsToLookNeighs = nextPointsToLook
    return

def addExtraV(line, neighIndex, dist):
    distTol = 0.70
    if 'extraVs' not in line: line['extraVs'] = set()
    line['extraVs'].add(neighIndex)
    #line['extraVs'][neighIndex] = (dist < distTol, dist)

def removeExtraV(line, neighIndex):
    line['extraVs'].remove(neighIndex)
    #line['extraVs'].pop(neighIndex)
    
def calcSize(lines, linesT):
    for line in lines:
        if 'parent' in line: continue
        size = 1
        if 'children' in line:
            for x in line['children']['lines']:
                if x[1] < len(lines):
                    size += 1
        if 'extraVs' in line:
            line['extraVs'] = list(line['extraVs'])
            size -= len(line['extraVs'])
        line['defectCount'] = size
    for line in linesT:
        if 'parent' in line: continue
        size = 0
        if 'extraVs' in line:
            line['extraVs'] = list(line['extraVs'])
            size -= len(line['extraVs'])
        line['defectCount'] = size

def addCollinearToTriads(cascade, cid, lines, coordIndexMap, linesT, coordIndexMapT, extraVs):
    latConst = cascade['latticeConst']
    angleTol = [15, 27]
    distTolErr = 0.10
    distTol = [0.30 + distTolErr, 0.50 + distTolErr, 0.70 + (distTolErr * 2)]#latConst * 0.3
    nn = NearestNeighbors(n_neighbors=5, radius=latConst * 2.0)
    nn.fit(cascade['eigen_features'][cid]['coords'])
    allFreeVs = []
    allLines = lines + linesT
    extraVsMap = {}
    for line in lines:
        if 'parent' in line: continue
        pointsToLookNeighs = [cascade['eigen_features'][cid]['coords'][x] for x in line['points']]
        IgnoreNeighs = set(line['points'])
        subLineVacs = [(line['vacErr'][1] + 1.0, line['points'][1])]
        isCheckSubLine = line['isVac'] != 2
        probableLines = {}
        while pointsToLookNeighs:
            dist, neigh = nn.radius_neighbors(pointsToLookNeighs)
            curNeigh = set()
            for x in neigh:
                for y in x:
                    if y in IgnoreNeighs: continue
                    curNeigh.add(y)
                    if y in coordIndexMap: 
                        IgnoreNeighs.update(lines[coordIndexMap[y]]['points'])
                        # not ignoring probable points which might be farther
                    elif y in coordIndexMapT: 
                        IgnoreNeighs.update(linesT[coordIndexMapT[y] - len(lines)]['points'])
                        # not ignoring probable points which might be farther
                    else:
                        IgnoreNeighs.add(y)
            nextPointsToLook = []
            #if (line['id'] == 10): (curNeigh)
            for neighIndex in curNeigh:
                neighIndex = int(neighIndex)
                if neighIndex in coordIndexMap:
                    otherLine = lines[coordIndexMap[neighIndex]]
                    if 'children' in otherLine: continue
                    #(neighIndex, line['id'])
                    dist = lineDist(line, otherLine)
                    if (dist[0] < angleTol[0] and (dist[1] < distTol[0] or dist[2] < distTol[0])):
                        dist = dist[1] + dist[2]
                        addLine(line, otherLine, dist)
                        nextPointsToLook.append(cascade['eigen_features'][cid]['coords'][neighIndex]) # TODO check index to add or something else
                        for y in otherLine['points']:
                            coordIndexMap[y] = line['id']
                        if not isCheckSubLine: isCheckSubLine = otherLine['isVac'] != 2
                        subLineVacs.append((dist, otherLine['points'][1]))
                elif neighIndex in coordIndexMapT:
                    otherLine = linesT[coordIndexMapT[neighIndex]- len(lines)]
                    #(neighIndex, line['id'])
                    dist = lineDist(line, otherLine)
                    if (dist[0] < angleTol[1] and (dist[1] < distTol[1] and dist[2] < distTol[1])):
                        #if (line['id'] == 1 and otherLine['id'] == 36): print("first", otherLine['id'], dist)
                        totalDist = dist[1] + dist[2]
                        addLine(line, otherLine, totalDist)
                        nextPointsToLook.append(cascade['eigen_features'][cid]['coords'][neighIndex]) # TODO check index to add or something else
                        for y in otherLine['points']:
                            coordIndexMap[y] = line['id']
                            coordIndexMapT.pop(y)# = line['id']                           
                        subLineVacs.append((dist[0], otherLine['points'][0]))
                        if dist[0] > distTol[0]: isCheckSubLine = True
                    else:
                        interstitialDist = line['eq'].distance_to(Point(otherLine['eq'].r2))
                        #if line['id'] == 1 and otherLine['id'] == 36:  print("check-2", otherLine['eq'].r2, interstitialDist)
                        if interstitialDist < distTol[2]: # interstitial inline
                            #if line['id'] == 1 and otherLine['id'] == 36:  print("second", otherLine['eq'].r2, interstitialDist)
                            probableLines[otherLine['points'][0]] = (otherLine['id'], interstitialDist)
                            nextPointsToLook.append(cascade['eigen_features'][cid]['coords'][otherLine['points'][1]])
                            isCheckSubLine = True
                            vacDist = line['eq'].distance_to(Point(otherLine['eq'].r))
                            subLineVacs.append((vacDist, otherLine['points'][0]))
                        # add to error                                              
                elif neighIndex in extraVs:
                    coordNeighIndex = cascade['coords'][cascade['clusters'][cid][neighIndex]][:3]
                    dist = line['eq'].distance_to(Point(coordNeighIndex))
                    if neighIndex not in extraVsMap or extraVsMap[neighIndex][0] > dist:
                        if neighIndex in extraVsMap: removeExtraV(lines[extraVsMap[neighIndex][1]], neighIndex)
                        addExtraV(line, neighIndex, dist)
                        extraVsMap[neighIndex] = (dist, line['id'])
            pointsToLookNeighs = nextPointsToLook
        subLinePs = set()
        freeVs = []
        if isCheckSubLine and len(subLineVacs) > 1:
            subLine = findSubLine(cascade, cid, line, subLineVacs)
            #if line['id'] == 13:  print(subLine)
            if subLine:
                #subLine['forceAlign'] = forceAlign(subLine['eq'])
                #subLine['points'].sort(key=sortPointsHelper)
                for subP in subLine['points']:
                    subLinePs.add(subP[1])
                    if subP[1] not in probableLines: continue
                    otherLineI = probableLines[subP[1]]
                    addLine(line, allLines[otherLineI[0]], otherLineI[1])
                    for y in allLines[otherLineI[0]]['points']:
                        coordIndexMap[y] = line['id']
                        if y in coordIndexMapT: coordIndexMapT.pop(y)
            line['subLine'] = subLine
        else:
            line['subLine'] = None
        curLinePs = list()
        curLinePs.append((0.0, line['points'][0]))
        curLinePs.append((0.0, line['points'][2]))                          
        if line['subLine'] == None or not line['points'][1] in subLinePs:
            if line['isVac'] == 0:
              freeVs.append(line['points'][1])
            else:
              curLinePs.append((line['vacErr'][1], line['points'][1]))
        if 'children' in line:
            for lineChIndex in line['children']['lines']:
                lineCh = allLines[lineChIndex[1]]
                for i, x in enumerate(lineCh['points']):
                    if line['subLine'] and x in subLinePs: continue
                    if i == 1 and 'isVac' in lineCh and lineCh['isVac'] == 0: 
                        freeVs.append(x)
                        continue
                    curLinePs.append((lineChIndex[0], x))
        allFreeVs += freeVs
        line['freeVs'] = freeVs
        line['mainPoints'] = curLinePs
        #line['mainPoints'].sort(key=sortPointsHelper)        
        #line['forceAlign'] = forceAlign(line['eq'])
    addCollinearToThresh(cascade, cid, lines, coordIndexMap, linesT, coordIndexMapT, nn, extraVs, extraVsMap)
    addCollinearToTriadsAgain(cascade, cid, lines, coordIndexMap, linesT, coordIndexMapT)
    addForceAlignAndSort(cascade, cid, lines, linesT)
    calcSize(lines, linesT)
    return allFreeVs, extraVsMap
    # TODO what if a line's neighbours are added while that line is not added?
    # possibly add Vs that don't add to subline to freeVs and points to the line?

def getThreshLines(cascade, cid, startingIndex, pairs):
    lines = []
    coordIndexMap = {}    
    if cascade['clusterSizes'][cid] < 0: return lines, coordIndexMap    
    for cCoordIndex, coordIndex in enumerate(cascade['clusters'][cid]):
        if cascade['coords'][coordIndex][5] == 1: continue
        if coordIndex not in pairs: continue
        otherIndex = pairs[coordIndex]
        cCoordIndex2 = cCoordIndex - 1
        #if cascade['coords'][coordIndex][3] == 0 or cascade['coords'][coordIndex][5] == 0: continue
        #coordIndex = cascade['clusters'][cid][cCoordIndex]
        p1 = Point(cascade['coords'][coordIndex][:3])
        p2 = Point(cascade['coords'][otherIndex][:3])
        line = Line(p1, p2)
        angles = getSymmetricPlaneAngle(line)
        includedPoints = [cCoordIndex2, cCoordIndex]
        coordIndexMap[cCoordIndex] = len(lines) + startingIndex
        coordIndexMap[cCoordIndex2] = len(lines) + startingIndex
        curLine = {'id': len(lines) + startingIndex, 'eq':line, 'points':includedPoints, 'angles':angles}
        lines.append(curLine)
    return lines, coordIndexMap


# In[1105]:


# not using currently
def getSymmetricCosineAngle(line):
    res = [0.0, 0.0, 0.0]    
    for i, x in enumerate(line.t):
        angle = abs(angular_unit*math.acos(x))
        if angle > 90: angle = 90 - angle
        if angle > 45: angle = abs(90 - angle)
        res[i] = angle
    return res

"""
triads represent dumbbells while pairs represent annihilated lattice-site and atom pair that may
be collinear to a dumbbell forming a crowdion.
Key is an atom, value is either a lattice site or a lattice site and annihilated atom in case of a dumbbell.
"""
def makeLatticeGroups(cascade):
    triads = {}
    pairs = {}
    pre = 0
    for cur in cascade['coDefects']:
        if (cur - pre) == 2:
            pairs[pre + 1] = pre
        elif (cur - pre) > 2:
            triads[pre + 2] = (pre, pre + 1)
        pre = cur
    return triads, pairs

def getPointDefectLines(cascade, triads, pairs):
    pointDefectLines = []
    for coordIndex, coord in enumerate(cascade['coords']):
      if coord[4] != -1: continue
      if coordIndex in triads:
        vacIndex = triads[coordIndex][0]
        secAtomIndex =  triads[coordIndex][1]
        secIndex = secAtomIndex
        points = [coordIndex, vacIndex, secAtomIndex]
      elif coordIndex in pairs:
        vacIndex = pairs[coordIndex]
        secIndex = vacIndex
        points = [coordIndex, vacIndex]
      else:
        if cascade['coords'][coordIndex][3] == 1 and cascade['coords'][coordIndex][5] == 1:
            pass # error shouldn't happen by definition of triads and -1 tag TODO log warning.
        continue
      p1 = Point(coord[:3])
      p2 = Point(cascade['coords'][secIndex][:3])
      line = Line(p1, p2)
      fa = forceAlign(line)
      pointDefectLines.append((points, fa['dir']))
    return pointDefectLines


def getTriadLines(cascade, cid, triads):
    lines = []
    coordIndexMap = {}
    probablePoints = {}
    if cascade['clusterSizes'][cid] < 0: return lines, coordIndexMap, probablePoints
    dumbbellPairSurviveIs = triads
    #dumbbellPairSurviveIs = set([triads[key][1] for key in triads])
    #dumbbellPairIs = set([triads[key][0] for x in triads]) + dumbbellPairSurviveIs
    cvmap = {x: i for i, x in enumerate(cascade['clusters'][cid]) if cascade['coords'][x][3] == 0}
    freeIs = []
    for cCoordIndex, coordIndex in enumerate(cascade['clusters'][cid]):
        #if cascade['coords'][coordIndex][3] == 0 or cascade['coords'][coordIndex][5] == 0: continue
        if coordIndex not in dumbbellPairSurviveIs: 
            if cascade['coords'][coordIndex][3] == 1 and cascade['coords'][coordIndex][5] == 1:
                freeIs.append(cCoordIndex)
            continue
        p1 = Point(cascade['coords'][coordIndex][:3])
        p2 = Point(cascade['coords'][triads[coordIndex][1]][:3])
        # TODO do we change perfect line or is it okay to store the wrong line Eq
        #line = Line(p1, Point(p1.r - p3)) if perfect else Line(p1, p2)
        line = Line(p1, p2)
        angles = getSymmetricPlaneAngle(line)     
        vacIndex = triads[coordIndex][0]
        if vacIndex not in cvmap:
          continue  
        vacClusterIndex = cvmap[vacIndex]
        latticeSite = Point(cascade['coords'][vacIndex][:3])      
        line1 = Line(p1, latticeSite)        
        line2 = Line(p2, latticeSite)
        vacAngle = abs(line1.angle_to(line2))
        includedPoints = [cCoordIndex, vacClusterIndex, cCoordIndex - 1]
        for pindex in includedPoints:
            coordIndexMap[pindex] = len(lines)
        isVac = 2
        if vacAngle < 18:
            #coordIndexMap[vacClusterIndex] = len(lines)
            isVac = 2
        elif vacAngle < 25:
            isVac = 1
        else:
            isVac = 0
        vacDist = line.distance_to(latticeSite)
        curLine = {'id': len(lines), 'eq':line, 'points':includedPoints, 'angles':angles, 'isVac': isVac, 'vacErr': (vacAngle, vacDist)}
        lines.append(curLine) #({[line1, line2]})
    return lines, coordIndexMap, freeIs#, probablePoints, freeVs

def getExtraFreeVs(cascade, cid, coordIndexMap, coordIndexMapT):
    extraVs = set()
    for cCoordIndex, coordIndex in enumerate(cascade['clusters'][cid]):
        if cCoordIndex not in coordIndexMap and cCoordIndex not in coordIndexMapT:
            if cascade['coords'][coordIndex][3] == 0: extraVs.add(cCoordIndex)
    return extraVs
