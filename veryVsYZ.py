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
	dataY = np.zeros( [grid.NY, grid.NZ] )
	dataZ = np.zeros( [grid.NY, grid.NZ] )
	data  = np.zeros( [grid.NY, grid.NZ] )
else:
	dataY = np.zeros( [grid.NZ, grid.NY] )
	dataZ = np.zeros( [grid.NZ, grid.NY] )
	data  = np.zeros( [grid.NZ, grid.NY] )

mpiX = mpiSliceX
for mpiZ in range( grid.PZ ):
	for mpiY in range( grid.PY ):
		fileY = open( "%s_X_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
		fileZ = open( "%s_X_mpi_%d_%d_%d.bin" % ( fileNameZ, mpiX, mpiY, mpiZ ), "rb" )
		file  = open( "%s_X_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		ny = grid.ny[mpiY]
		nz = grid.nz[mpiZ]
		print( "ny = %d, nz = %d" % ( ny, nz ) )
		datay = np.fromfile( fileY, dtype='float32', count = ny * nz )
		dataz = np.fromfile( fileZ, dtype='float32', count = ny * nz )
		data_ = np.fromfile( file , dtype='float32', count = ny * nz )
		
		J  = grid.frontNY[mpiY]
		J_ = grid.frontNY[mpiY] + ny
		K  = grid.frontNZ[mpiZ]
		K_ = grid.frontNZ[mpiZ] + nz

		if FAST_AXIS == 'Z':
			dataY[J:J_, K:K_] = np.reshape( datay, ( ny, nz ) )
			dataZ[J:J_, K:K_] = np.reshape( dataz, ( ny, nz ) )
			data [J:J_, K:K_] = np.reshape( data_, ( ny, nz ) )
		else:
			dataY[K:K_, J:J_] = np.reshape( datay, ( nz, ny ) )
			dataZ[K:K_, J:J_] = np.reshape( dataz, ( nz, ny ) )
			data [K:K_, J:J_] = np.reshape( data_, ( nz, ny ) )


dpi = 300
fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
plt.pcolor( dataY[::sample, ::sample] // unit, dataZ[::sample, ::sample] // unit, data [::sample, ::sample] // unit, cmap = "jet" )
plt.colorbar( )
plt.axis( "image" )
plt.savefig( "VsYZ.png" )
#plt.plot( dataZ[grid.NZ - 1, :] )

#print( grid )


