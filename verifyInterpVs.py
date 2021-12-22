#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 12, 2021
15:32
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt



jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 2

k = 8

outputPath = params["out"]

lonFileName = params["out"] + "/lon"
latFileName = params["out"] + "/lat"
interpVsFileName = params["out"] + "/InterpVs%d"%k
#interpVsFileName = params["out"] + "/Vs"

FAST_AXIS = params["FAST_AXIS"]
'''
for Z in range( grid.PZ ):
	for Y in range( grid.PY ):
		for X in range( grid.PX ):
			print( "nx = %d, ny = %d, nz = %d\n" % ( grid.nx[X], grid.ny[Y], grid.nz[Z] ) )
'''

sliceX = params["sliceX"] - grid.frontNX
sliceY = params["sliceY"] - grid.frontNY
sliceZ = params["sliceZ"] - grid.frontNZ



if FAST_AXIS == 'Z':
	lon = np.zeros( 	[grid.NX, grid.NY] )
	lat = np.zeros( 	[grid.NX, grid.NY] )
	interpVs = np.zeros( [grid.NX, grid.NY] )
else:
	lon = np.zeros( [grid.NY, grid.NX] )
	lat = np.zeros( [grid.NY, grid.NX] )
	interpVs = np.zeros( [grid.NY, grid.NX] )


mpiZ = grid.PZ - 1
for mpiY in range( grid.PY ):
	for mpiX in range( grid.PX ):
		lonFile = open( "%s_mpi_%d_%d_%d.bin" % ( lonFileName, mpiX, mpiY, mpiZ ), "rb" )
		latFile = open( "%s_mpi_%d_%d_%d.bin" % ( latFileName, mpiX, mpiY, mpiZ ), "rb" )
		interpVsFile = open( "%s_mpi_%d_%d_%d.bin" % ( interpVsFileName, mpiX, mpiY, mpiZ ), "rb" )
		#interpVsFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( interpVsFileName, mpiX, mpiY, mpiZ ), "rb" )
		
		ny = grid.ny[mpiY]
		nx = grid.nx[mpiX]

		print( "ny = %d, nx = %d" % ( ny, nx ) )
		lon_ = np.fromfile( lonFile, dtype='float32', count = ny * nx )
		lat_ = np.fromfile( latFile, dtype='float32', count = ny * nx )
		interpVs_ = np.fromfile( interpVsFile, dtype='float32', count = ny * nx )

		J  = grid.frontNY[mpiY]
		J_ = grid.frontNY[mpiY] + ny
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx
		if FAST_AXIS == 'Z':
			lon		[I:I_, J:J_] = np.reshape( lon_, 		(nx, ny) )
			lat		[I:I_, J:J_] = np.reshape( lat_, 		(nx, ny) )
			interpVs [I:I_, J:J_] = np.reshape( interpVs_, 	(nx, ny) )
		else:
			lon[J:J_, I:I_] = np.reshape( lon_, (ny, nx) )
			lat[J:J_, I:I_] = np.reshape( lat_, (ny, nx) )
			interpVs[J:J_, I:I_] = np.reshape( interpVs_, (ny, nx) )

nPML = params["nPML"]
NX = grid.NX
NY = grid.NY

np.save( "lon.npy", lon[nPML:NY - nPML, nPML:NX-nPML] )
np.save( "lat.npy", lat[nPML:NY - nPML, nPML:NX-nPML] )

dpi = 300
fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1000 # 1km = 1000m
data = np.zeros( [grid.NY, grid.NX] )
#plt.pcolor( lon[::sample, ::sample], lat[::sample, ::sample], interpVs[::sample, ::sample], cmap = "jet" )
plt.contourf( lon[::sample, ::sample], lat[::sample, ::sample], interpVs[::sample, ::sample], 80, cmap = "jet" )
plt.colorbar( )
plt.axis( "image" )
plt.savefig( "interpVs.png" )

#plt.plot( dataZ[grid.NZ - 1, :] )
