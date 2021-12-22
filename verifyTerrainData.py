#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 05, 2021

'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import struct



jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

blockX = params["blockX"]
blockY = params["blockY"]

pointsPerBlockX = 6000
pointsPerBlockY = 6000

totalPointX = ( pointsPerBlockX - 1 ) * blockX + 1;
totalPointY = ( pointsPerBlockY - 1 ) * blockY + 1;

totalTerrain = np.zeros( [totalPointY, totalPointX] )

totalTerrainFile = open( "./TerrainDir/totalTerrain.bin", "rb" )
#totalTerrainFile = open( "./TerrainDir/gradTotalTerrain.bin", "rb" )

data = np.fromfile( totalTerrainFile, dtype='float32', count = totalPointX * totalPointY )

totalTerrain = np.reshape( data, ( totalPointY, totalPointX ) )

plt.pcolor( totalTerrain[::10, ::10], cmap = "seismic" )
plt.colorbar( )


