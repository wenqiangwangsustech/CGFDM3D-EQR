#!/usr/bin/env python

'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np
from pyscripts.GRID import GRID
import matplotlib.pyplot as plt
import  sys


FreeSurf = 1

it = 800

var = 'Vx' #Vy Vz



if len( sys.argv ) > 1:
	it = int( sys.argv[1] )
	var = str( sys.argv[2] )
	FreeSurf = int( sys.argv[3] )


jsonsFile = open( "params.json" )
params = json.load( jsonsFile )
grid = GRID( params )
FAST_AXIS = params["FAST_AXIS"]

sample = 1

outputPath = params["out"]
fileNameX = params["out"] + "/coordX"
fileNameY = params["out"] + "/coordY"


if FreeSurf == 1:
	fileName = params["out"] + "/FreeSurf%s_%d" % ( var, it )
	fileName2 = "output_homo_original"+ "/FreeSurf%s_%d" % ( var, it )
    #fileName = params["out"] + "/PGV"
	varname = "./png/compFreeSurf%s_%d" % ( var, it )
else:
	fileName = params["out"] + "/%s_%d" % ( var, it )
	fileName2 = "output_homo_original"+ "/%s_%d" % ( var, it )
	varname = "./png/comp%s_%d" % ( var, it )


print( "Draw " + fileName  )



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
	dataX = np.zeros( [grid.NX, grid.NY] )
	dataY = np.zeros( [grid.NX, grid.NY] )
	data  = np.zeros( [grid.NX, grid.NY] )
	data2  = np.zeros( [grid.NX, grid.NY] )
else:
	dataX = np.zeros( [grid.NY, grid.NX] )
	dataY = np.zeros( [grid.NY, grid.NX] )
	data  = np.zeros( [grid.NY, grid.NX] )
	data2  = np.zeros( [grid.NY, grid.NX] )





for mpiY in range( grid.PY ):
	for mpiX in range( grid.PX ):
		mpiZ = mpiSliceZ
		XFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameX, mpiX, mpiY, mpiZ ), "rb" )
		YFile = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileNameY, mpiX, mpiY, mpiZ ), "rb" )
		if FreeSurf:
			mpiZ = grid.PZ - 1
		else:
			mpiZ = mpiSliceZ
		File  = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName , mpiX, mpiY, mpiZ ), "rb" )
		File2 = open( "%s_Z_mpi_%d_%d_%d.bin" % ( fileName2 , mpiX, mpiY, mpiZ ), "rb" )
		
		ny = grid.ny[mpiY]
		nx = grid.nx[mpiX]

		print( "ny = %d, nx = %d" % ( nx, ny ) )
		datax = np.fromfile( XFile, dtype='float32', count = ny * nx )
		datay = np.fromfile( YFile, dtype='float32', count = ny * nx )
		data_ = np.fromfile(  File, dtype='float32', count = ny * nx )
		data_2 = np.fromfile(  File2, dtype='float32', count = ny * nx )

		J  = grid.frontNY[mpiY]
		J_ = grid.frontNY[mpiY] + ny
		I  = grid.frontNX[mpiX]
		I_ = grid.frontNX[mpiX] + nx

		if FAST_AXIS == 'Z':
			dataX[I:I_, J:J_] = np.reshape( datax, (nx, ny) )
			dataY[I:I_, J:J_] = np.reshape( datay, (nx, ny) )
			data [I:I_, J:J_] = np.reshape( data_, (nx, ny) )
			data2 [I:I_, J:J_] = np.reshape( data_2, (nx, ny) )
		else:
			dataX[J:J_, I:I_] = np.reshape( datax, (ny, nx) )
			dataY[J:J_, I:I_] = np.reshape( datay, (ny, nx) )
			data [J:J_, I:I_] = np.reshape( data_, (ny, nx) )
			data2 [J:J_, I:I_] = np.reshape( data_2, (ny, nx) )


dpi = 300
#fig = plt.figure( dpi = dpi, figsize = ( 1920 // dpi, 1080 // dpi ) )
unit = 1 # 1km = 1000m
vm = np.max( np.abs( data - data2 ) )
#plt.contourf( dataX, dataY, data, 10 )
#plt.pcolormesh( data[::sample, ::sample] // unit, cmap = "seismic" )
plt.pcolormesh( dataX, dataY, np.abs( data - data2 ), vmax = vm, vmin = 0.0, cmap = "seismic" ) #origin= "lower" )
#plt.pcolormesh( dataX, dataY, data, cmap = "jet" ) #origin= "lower" )
plt.colorbar( )
plt.axis( "image" )
plt.title( fileName )

#plt.show( )
plt.savefig( varname + ".png" )
#plt.plot( dataY[grid.NZ - 1, :] )

#print( grid )


