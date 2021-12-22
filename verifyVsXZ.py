#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt





jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )

sample = 1

outputPath = params["out"]
fileNameX = params["out"] + "/coordX"
fileNameY = params["out"] + "/coordY"
fileNameZ = params["out"] + "/coordZ"


fileName = params["out"] + "/Vs"


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

for mpiSliceX in range( grid.PX ):
	if sliceX[mpiSliceX] >= 0 and sliceX[mpiSliceX] < grid.nx[mpiSliceX]:
		break

for mpiSliceY in range( grid.PY ):
	if sliceY[mpiSliceY] >= 0 and sliceY[mpiSliceY] < grid.ny[mpiSliceY]:
		break

for mpiSliceZ in range( grid.PZ ):
	if sliceZ[mpiSliceZ] >= 0 and sliceZ[mpiSliceZ] < grid.nz[mpiSliceZ]:
		break



if FAST_AXIS == 'Z':
	dataX = np.zeros( [grid.NX, grid.NZ] )
	dataZ = np.zeros( [grid.NX, grid.NZ] )
	data  = np.zeros( [grid.NX, grid.NZ] )
else:
	dataX = np.zeros( [grid.NZ, grid.NX] )
	dataZ = np.zeros( [grid.NZ, grid.NX] )
	data  = np.zeros( [grid.NZ, grid.NX] )


mpiY = mpiSliceY
for mpiZ in range( grid.PZ ):
	for mpiX in range( grid.PX ):
		fileX = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_Y_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		nx = grid.nx[mpiX]
		nz = grid.nz[mpiZ]
		print( "nx = %d, nz = %d" % ( nx, nz ) )
		datax = np.fromfile( fileX, dtype='float32', count = nx * nz )
		print( np.shape( datax ) )
		dataz = np.fromfile( fileZ, dtype='float32', count = nx * nz )
		data_  = np.fromfile( file, dtype='float32', count = nx * nz )
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx
		K  = grid.frontNZ[mpiZ]
		K_ = grid.frontNZ[mpiZ] + nz

		if FAST_AXIS == 'Z':
			dataX[I:I_, K:K_] = np.reshape( datax, ( nx, nz ) )
			dataZ[I:I_, K:K_] = np.reshape( dataz, ( nx, nz ) )
			data [I:I_, K:K_] = np.reshape( data_, ( nx, nz ) )
		else:
			dataX[K:K_, I:I_] = np.reshape( datax, ( nz, nx ) )
			dataZ[K:K_, I:I_] = np.reshape( dataz, ( nz, nx ) )
			data [K:K_, I:I_] = np.reshape( data_, ( nz, nx ) )


dpi = 300
fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
plt.pcolormesh( dataX[::sample, ::sample] // unit, dataZ[::sample, ::sample] // unit, data[::sample, ::sample] // unit, cmap = "jet" )
plt.colorbar( )
plt.axis( "image" )
plt.savefig( "VsXZ.png" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )


