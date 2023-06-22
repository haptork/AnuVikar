import sys
import os
pathToCsaranshPP = "/home/utkarsh/code/AnuVikar/" # change if required
sys.path.append(pathToCsaranshPP)
from pysrc.anuvikar_cdb_helper import getDefaultConfig, getDefaultInfos, mpiProcessXyzFilesInDirGivenInfo
buildDir = os.path.join(pathToCsaranshPP, "_build")
libPath = os.path.join(buildDir, "libanuvikar_shared.so")  # path to anuvikar library
if (not os.path.exists(buildDir) or not os.path.exists(libPath)):
    print("Library not found. Might be due to build errors in cmake.")
    print("If built successfully, specify correct build directory & lib file (so / dlib / dlib) above.")
    
import plotly as py
import plotly.express as px
import json
cascades = None
debug = False
data = None

def isInterstitial(coord):
    return 1 == coord[3]

def isSurviving(coord):
    return 1 == coord[5]

config = getDefaultConfig()
config['logFilePath'] = "local-log.txt"
config['anuvikarLib'] = libPath
info, extraInfo = getDefaultInfos()
extraInfo['substrate'] = "Fe"
extraInfo['isPKAGiven'] = False
info['structure'] = "bcc"
info['isIgnoreBoundary'] = False
info['xyzFileType'] = "Generic"
info['xyzColumnStart'] = 2
info["latticeConst"] = 2.86
#info["latticeConst"] = 3.18#652
extraInfo['energy'] = 20
#info['originType'] = 1
info['extraColumnStart'] = -2
extraInfo["author"] = "Andrea"
#xyzDir = "../T300K_50keV_data/"
xyzDir = "/home/utkarsh/code/aviml/parallelize/006-fpos-20/"
#xyzDir = "/media/utkarsh/data/andreas/3-Derlet-Dudarev/fpos"
cascades = mpiProcessXyzFilesInDirGivenInfo(xyzDir, info, extraInfo, config, prefix=[], suffix=["z"], outputJsonFilename="out.json")
