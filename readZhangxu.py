#!/usr/bin/python
'''
Wenqiang Wang@SUSTech 

Date: 2021/12/21

'''


import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import json as js
import struct
from pyproj import Proj
from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator

def gaussFunc( t, t0, a ):
    source = np.exp(-( t - t0 )**2 / ( a ** 2 ) ) / ( np.sqrt( np.pi ) * a )
    return source



#right-hand rule
def rotateX( xTheta, x, y, z ):
	newX = x
	newY = x * 0 + y * np.cos( xTheta ) - z * np.sin( xTheta )
	newZ = x * 0 + y * np.sin( xTheta ) + z * np.cos( xTheta ) 
	return newX, newY, newZ

def rotateY( yTheta, x, y, z ):
	newX = x * np.cos( yTheta ) + z * np.sin( yTheta )
	newY = y
	newZ =-x * np.sin( yTheta ) + z * np.cos( yTheta )
	return newX, newY, newZ

def rotateZ( zTheta, x, y, z ):
	newX = x * np.cos( zTheta ) - y * np.sin( zTheta )
	newY = x * np.sin( zTheta ) + y * np.cos( zTheta )
	newZ = z
	return newX, newY, newZ
	

def localElementMesh( Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, n1, n2 ):
	xmesh = np.zeros( [n2, n1] ) 
	ymesh = np.zeros( [n2, n1] ) 
	zmesh = np.zeros( [n2, n1] ) 

	xmesh[:, 0] = np.linspace( Ax, Cx, n2 )
	ymesh[:, 0] = np.linspace( Ay, Cy, n2 )
	zmesh[:, 0] = np.linspace( Az, Cz, n2 )

	xmesh[:, -1] = np.linspace( Bx, Dx, n2 )
	ymesh[:, -1] = np.linspace( By, Dy, n2 )
	zmesh[:, -1] = np.linspace( Bz, Dz, n2 )

	for j in range( n2 ):
		xmesh[j, :] = np.linspace( xmesh[j,0], xmesh[j,-1], n1 ) 
		ymesh[j, :] = np.linspace( ymesh[j,0], ymesh[j,-1], n1 ) 
		zmesh[j, :] = np.linspace( zmesh[j,0], zmesh[j,-1], n1 )

	return xmesh, ymesh, zmesh 
	#return xmesh[:-1, :-1], ymesh[:-1, :-1], zmesh[:-1, :-1] 

def thickenMesh( nx, ny, thickenX, thickenY, X, Y, Z, XMesh, YMesh, ZMesh ):
	NX = ( nx - 1 ) * thickenX + 1
	NY = ( ny - 1 ) * thickenY + 1
	for j in range( 1, ny ):
		for i in range( 1, nx ):
			J_ = ( j - 1 ) * thickenY
			J  = j * thickenY
			I_ = ( i - 1 ) * thickenX
			I  = i * thickenX
			index = i - 1 + ( j - 1 ) * nx
			Ax, Ay, Az = X[index], Y[index], Z[index] 
			index = i + ( j - 1 ) * nx
			Bx, By, Bz = X[index], Y[index], Z[index] 
			index = i - 1 + j * nx
			Cx, Cy, Cz = X[index], Y[index], Z[index]
			index = i + j * nx
			Dx, Dy, Dz = X[index], Y[index], Z[index]
			n1 = thickenX
			n2 = thickenY

			XCoord, YCoord, ZCoord = localElementMesh( Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, Dx, Dy, Dz, n1 + 1, n2 + 1 )
			for jj in range( J_, J + 1 ):
				for ii in range( I_, I + 1 ):
					index = ii + jj * NX
					XMesh[index], YMesh[index], ZMesh[index] = XCoord[jj - J_, ii - I_], YCoord[jj - J_, ii - I_], ZCoord[jj - J_, ii - I_], 


def dealWithBoundaryValue( NX, NY, DataMesh ):
	i = NX - 1
	for j in range( NY ):
		index  = i + j * NX
		index2 = i - 1 + j * NX
		DataMesh[index] = DataMesh[index2]

	j = NY - 1
	for i in range( NX ):
		index  = i + j * NX
		index2 = i + ( j - 1 ) * NX
		DataMesh[index] = DataMesh[index2]


def dealNanValue( NX, NY, DataMesh ):
	StoreNan = np.isnan( DataMesh )
	Nnan = 0
	#print( StoreNan )
	
	for i in range( NX * NY ):
		if StoreNan[i] == True:
			Nnan += 1
			DataMesh[i] = 0.0
	#print( "The total Nan value number is %d\n" % Nnan )
		      


NT = 1000 


lengthPoint = 14
depthPoint =  8


thickenX = 20
thickenZ = 20

zy_data = np.loadtxt("ZhangxuModel.txt")
NT = len( zy_data[0,10:] )

zy_lat         = np.zeros( [ lengthPoint * depthPoint ] )
zy_lon         = np.zeros( [ lengthPoint * depthPoint ] )
zy_depth       = np.zeros( [ lengthPoint * depthPoint ] )
zy_size_strike = np.zeros( [ lengthPoint * depthPoint ] )
zy_size_dip    = np.zeros( [ lengthPoint * depthPoint ] )
zy_strike      = np.zeros( [ lengthPoint * depthPoint ] )
zy_dip         = np.zeros( [ lengthPoint * depthPoint ] )
zy_rake        = np.zeros( [ lengthPoint * depthPoint ] )
zy_slip        = np.zeros( [ lengthPoint * depthPoint ] )
zy_sample_rate = np.zeros( [ lengthPoint * depthPoint ] )
zy_slip_rate   = np.zeros( [lengthPoint * depthPoint, NT] )


for j in range( depthPoint ):
	for i in range( lengthPoint ):
		index = i + j * lengthPoint 
		pos =  j + i * depthPoint 
		zy_lat         [index] = zy_data[pos,0]
		zy_lon         [index] = zy_data[pos,1]
		zy_depth       [index] = zy_data[pos,2] * (-1000)
		zy_size_strike [index] = zy_data[pos,3] * 1000
		zy_size_dip    [index] = zy_data[pos,4] * 1000
		zy_strike      [index] = zy_data[pos,5]
		zy_dip         [index] = zy_data[pos,6]
		zy_rake        [index] = zy_data[pos,7]
		zy_slip        [index] = zy_data[pos,8]
		zy_sample_rate [index] = zy_data[pos,9]
		zy_slip_rate   [index] = zy_data[pos,10:]



strike = zy_strike[0]
dip = zy_dip[0]


src_npts = len(zy_data)
dt = 1 / zy_sample_rate[0]
lon0 = ( np.max( zy_lon ) + np.min( zy_lon ) )* 0.5
lat0 = ( np.max( zy_lat ) + np.min( zy_lat ) )* 0.5

proj = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84" )
XCenter, YCenter = proj( zy_lon, zy_lat )
ZCenter = zy_depth 


fig = plt.figure( 1 )
ax = Axes3D( fig )
ax.scatter( XCenter[:10], YCenter[:10], ZCenter[:10] )


nx = lengthPoint + 1
nz = depthPoint + 1
npts = nz * nx

lon  = np.zeros( npts )
lat  = np.zeros( npts )
depth = np.zeros( npts )

'''
X = np.zeros( npts )
Y = np.zeros( npts )
Z =  np.zeros( npts )


slip = np.zeros( npts  )
slipRate = np.zeros( [ npts, NT]  )
rake = np.zeros( npts  )
'''

xTheta = 0#- 45 * np.pi / 180
yTheta = -dip * np.pi / 180#45 * np.pi / 180
#yTheta = 0#45 * np.pi / 180
zTheta = strike * np.pi / 180

XCenter, YCenter, ZCenter = rotateZ( zTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateX( xTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateY( yTheta, XCenter, YCenter, ZCenter )

fig = plt.figure( 2)
ax = Axes3D( fig )
ax.scatter( XCenter, YCenter, ZCenter )






'''
disXmove = 0.5 * zy_size_strike[0]#( XCenter[lengthPoint] - XCenter[0] ) * 0.5
disYmove = 0.5 * zy_size_dip[0]#( YCenter[1] - YCenter[0] ) * 0.5




for j in range( nz - 1 ):
	for i in range( nx - 1 ):
		index = i + j * nx
		pos = i + j * lengthPoint
		X[index] = XCenter[pos] - disXmove
		Y[index] = YCenter[pos] - disYmove
fig = plt.figure( 6 )
ax = Axes3D( fig )
ax.scatter( X, Y, Z )		


j = nz - 1
for i in range( 0, nx - 1 ):
	index = i + j * nx
	pos = i + ( j - 1 ) * nx
	X[index] = X[pos] + 2 * disXmove
	Y[index] = Y[pos] #disYmove  


i = nx - 1
for j in range( 0, nz ):
	index = i + j * nx
	pos = i - 1 + j * nx
	X[index] = X[pos]
	Y[index] = Y[pos] + 2 * disYmove


fig = plt.figure( 15 )
ax = Axes3D( fig )
ax.scatter( X, Y, Z )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), ZCenter )
Z = interp( X, Y )


interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), zy_slip )
slip = interp( X, Y )

for it in range( NT ):
	interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), zy_slip_rate[:, it] )
	slipRate[:, it] = interp( X, Y )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), zy_rake )
rake = interp( X, Y )
'''


NX = ( nx - 1 ) * thickenX + 1
NY = ( nz - 1 ) * thickenZ + 1
NPTS = NX * NY
XMesh = np.zeros( NPTS )
YMesh = np.zeros( NPTS )
ZMesh = np.zeros( NPTS )

XMesh2D = np.zeros( [NY, NX] )
YMesh2D = np.zeros( [NY, NX] )
DataMesh2D = np.zeros( [NY, NX] )



SlipMesh = np.zeros( NPTS )
AreaMesh = np.zeros( NPTS )
rakeMesh = np.zeros( NPTS )
SlipRateMesh = np.zeros( [NPTS, NT] )
IntegaralSlipMesh = np.zeros( [NPTS] )
RakeTimeMesh = np.zeros( [NPTS, NT] )
StrikeMesh = np.zeros( NPTS ) + strike
DipMesh = np.zeros( NPTS  ) + dip


#thickenMesh( nx, nz, thickenX, thickenZ, X, Y, Z, XMesh, YMesh, ZMesh )

Xmin = np.min( XCenter )
Xmax = np.max( XCenter )

Ymin = np.min( YCenter )
Ymax = np.max( YCenter )

XDH = ( Xmax - Xmin ) / ( NX - 1 )
YDH = ( Ymax - Ymin ) / ( NY - 1 )


for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XMesh[index] = Xmin + XDH * i
		YMesh[index] = Ymin + YDH * j
		XMesh2D[j, i] = Xmin + XDH * i
		YMesh2D[j, i] = Ymin + YDH * j

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), ZCenter )
ZMesh = interp( XMesh, YMesh )


'''
for j in range( 1, NY ):
	for i in range( 1, NX ):
		index2 = i + j * NX
		index  = i - 1 + ( j - 1 ) * NX
		x = np.abs( XMesh[index2] - XMesh[index] )
		y = np.abs( YMesh[index2] - YMesh[index] )
		AreaMesh[index] = x * y

i = NX - 1
for j in range( NY ):
	index  = i + j * NX
	index2 = i - 1 + j * NX
	AreaMesh[index] = AreaMesh[index2]



j = NY - 1
for i in range( NX ):
	index  = i + j * NX
	index2 = i + ( j - 1 ) * NX
	AreaMesh[index] = AreaMesh[index2]




i = NX - 1
for j in range( NY ):
	index  = i + j * NX
	index2 = i - 1 + j * NX
	SlipMesh[index] = SlipMesh[index2]



j = NY - 1
for i in range( NX ):
	index  = i + j * NX
	index2 = i + ( j - 1 ) * NX
	SlipMesh[index] = SlipMesh[index2]
'''



AreaMesh = np.zeros( NPTS ) + XDH * YDH


interp = CloughTocher2DInterpolator( list(zip(XCenter, YCenter)), zy_slip)
SlipMesh = interp( XMesh, YMesh )
dealNanValue( NX, NY, SlipMesh )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), zy_rake )
rakeMesh = interp( XMesh, YMesh )




for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		DataMesh2D[j, i] = SlipMesh[index]
        
        
fig = plt.figure( 3 )
plt.pcolormesh( XMesh2D, YMesh2D, DataMesh2D, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )



for i in range( NY * NX ):
	for it in range( NT ):
		RakeTimeMesh[i, it] = rakeMesh[i]

for it in range( 0, NT ):
	interp = CloughTocher2DInterpolator(list(zip(XCenter, YCenter)), zy_slip_rate[:,it])
	SlipRateMesh[:, it] = interp( XMesh, YMesh)
	dealNanValue( NX, NY, SlipRateMesh[:, it] )



for i in range( NY * NX ):
	for it in range( 0, NT ):
		if SlipRateMesh[i, it] < 0. :
			SlipRateMesh[i, it] = 0.0

'''
if NORM == 1:
	svmax = np.max( np.abs( SlipRateMesh ) )
	for i in range( NPTS ):    
		SlipRateMesh[i, :] = SlipRateMesh[i, :] / svmax * SlipMesh[i] 
'''

mu = 3e10
SlipRateMesh = SlipRateMesh / ( zy_size_dip[0] * zy_size_strike[0] ) / mu
IntegaralSlipMesh = IntegaralSlipMesh / ( zy_size_dip[0] * zy_size_strike[0] ) / mu


for i in range( NY * NX ):
	sum0 = 0.
	for it in range( 1, NT ):
		sum0 += ( SlipRateMesh[i, it - 1] + SlipRateMesh[i, it - 1] ) * dt * 0.5
	IntegaralSlipMesh[i] = sum0 
	





XMesh, YMesh, ZMesh = rotateY( -yTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateX( -xTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateZ( -zTheta, XMesh, YMesh, ZMesh )



#fig = plt.figure( 14 )
#ax = Axes3D( fig )
#ax.scatter( XMesh, YMesh, ZMesh )

lonMesh = np.zeros( NPTS )
latMesh = np.zeros( NPTS )

proj1 = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84")
lonMesh, latMesh = proj1( XMesh, YMesh, inverse = "true" )


jsonsFile = open( "params.json" )
params = js.load( jsonsFile )

sourceFileName = params["sourceDir"] + "/zhangxu_source.bin"
sourceFile = open( sourceFileName, "wb" )

value = struct.pack( "i", NPTS )


print( "NPTS = %d\n" % NPTS  )
sourceFile.write( value )

print( "NT = %d\n" % NT  )
value = struct.pack( "i", NT )
sourceFile.write( value )

print( "dt = %f\n" % dt  )
value = struct.pack( "f", dt )
sourceFile.write( value )

srcInfo = open( "source/zhangxuTxt.txt", "w" )

for i in range( NPTS ):
	srcInfo.write( "Lon: %f, Lat: %f, Z: %f\n" % ( lonMesh[i], latMesh[i], ZMesh[i]  )  )
	value = struct.pack( "f",  lonMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  latMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  ZMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  AreaMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  StrikeMesh[i] )
	sourceFile.write( value )
	value = struct.pack( "f",  DipMesh[i] )
	sourceFile.write( value )
	tvalue = struct.pack( "f" * NT,  *( RakeTimeMesh[i, :] ) )
	sourceFile.write( tvalue )
	tvalue = struct.pack( "f" * NT,  *( SlipRateMesh[i, :] ) )
	sourceFile.write( tvalue )

sourceFile.close( )
srcInfo.close( )

XCoord = np.zeros( [NY, NX] )
YCoord = np.zeros( [NY, NX] )
ZCoord = np.zeros( [NY, NX] )

data = np.zeros( [NY, NX] )
dataSlip = np.zeros( [NY, NX] )
dataArea = np.zeros( [NY, NX] )
dataRake = np.zeros( [NY, NX] )
dataSlipRate = np.zeros( [NY, NX, NT] )

for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XCoord[j, i] = lonMesh[index]
		YCoord[j, i] = latMesh[index]
		ZCoord[j, i] = ZMesh[index]
		dataSlip[j, i] = SlipMesh[index]
		#dataSlip[j, i] = IntegaralSlipMesh[index]
		dataArea[j, i] = AreaMesh[index]
		dataRake[j, i] = rakeMesh[index]
		dataSlipRate[j, i, :] = SlipRateMesh[index, :]
		

fig = plt.figure( 4)
#plt.pcolormesh( XCoord, YCoord, dataSlip, cmap = "seismic" )
plt.pcolormesh( XCoord, YCoord, dataSlip, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )


'''
fig = plt.figure( 12 )
plt.pcolormesh( XCoord, YCoord, dataSlip, cmap = "jet" )
#plt.pcolormesh( LonCoord, LatCoord, dataSlipRate, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

plt.savefig( "fault.png" )


fig = plt.figure( 4 )
ax = Axes3D( fig )
ax.scatter( X, Y, slip )
XCoord = np.zeros( [NY, NX]  )
YCoord = np.zeros( [NY, NX]  )
ZCoord = np.zeros( [NY, NX]  )

data = np.zeros( [NY, NX] )
dataSlip = np.zeros( [NY, NX] )
dataArea = np.zeros( [NY, NX] )
dataRake = np.zeros( [NY, NX] )
dataSlipRate = np.zeros( [NY, NX, NT] )

for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XCoord[j, i] = XMesh[index]
		YCoord[j, i] = YMesh[index]
		ZCoord[j, i] = ZMesh[index]
		data[j, i] = SlipMesh[index]
		dataSlip[j, i] = SlipMesh[index]
		dataSlip[j, i] = IntegaralSlipMesh[index]
		dataArea[j, i] = AreaMesh[index]
		dataRake[j, i] = rakeMesh[index]
		dataSlipRate[j, i, :] = SlipRateMesh[index, :]
		

fig = plt.figure( )
plt.pcolormesh( LonCoord, LatCoord, dataSlip, cmap = "seismic" )
#plt.pcolormesh( LonCoord, LatCoord, dataSlipRate, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

plt.savefig( "fault.png" )
'''



'''
plt.ion( )
plt.figure( )
for it in range( 0, NT, 10 ):
	plt.clf( )
	vm = np.max( dataSlipRate[:, :, it] )
	#vm = 0.1
	plt.pcolormesh( dataSlipRate[:, :, it], vmin = - vm, vmax = vm, cmap = "seismic")
	plt.xlabel( "x" )
	plt.ylabel( "y" )
	plt.colorbar( )
	plt.axis( "image" )
	plt.show( )
	plt.pause( 0.002 )
'''
#fig = plt.figure( 1 )
#plt.pcolormesh( Slip, cmap = "jet" )
#plt.colorbar( )
#plt.axis( "image" )









#fig = plt.figure( 0 )
#plt.gca().set_box_aspect( ( 250, 130, 4)  )




#fig = plt.figure( 1 )
#plt.scatter( XCoord, YCoord, XCoord )
#plt.pcolormesh( XCoord, YCoord, data, cmap = "jet" )
#plt.pcolormesh( XCoord, YCoord, abs( dataRake ),  cmap = "jet" )

#plt.contourf( XCoord, YCoord, abs( dataArea ), label = 15,  cmap = "jet" )

#plt.colorbar( )
#plt.pcolormesh( XCoord[:-1, :-1], YCoord[:-1, :-1], XCoord[:-1, :-1], cmap = "jet" )
#plt.axis( "equal" )

#plt.scatter( X, Y, slip )


