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


def sourceTimeFunc( r, nt ):
	ts = np.linspace( 0, r, nt )
	s = ( 1 - np.cos( 2 * np.pi * ts / r ) )/ r
	return s


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



NT = 1000
strike = 289
dip = 84


lengthPoint = 27
depthPoint = 12


thickenX = 10
thickenZ = 10

Dx = 2000
Dy = 2000

src_npts = lengthPoint * depthPoint

usgsData = np.loadtxt( "basic_inversion.param" )

rakeCenter 	  = np.zeros( src_npts )

latRectCenter	= np.zeros( src_npts )
lonRectCenter	= np.zeros( src_npts )
depthRectCenter = np.zeros( src_npts )

timeFieldCenter	= np.zeros( src_npts )
halfDurationCenter = np.zeros( src_npts ) 
slipCenter = np.zeros( src_npts )


for j in range( src_npts ):
	lonRectCenter  [j]    = usgsData[j, 1]
	latRectCenter  [j]    = usgsData[j, 0]
	depthRectCenter[j]    = usgsData[j, 2] * 1000
	slipCenter[j]		  = usgsData[j, 3] * 0.01
	rakeCenter[j]	      = usgsData[j, 4]
	timeFieldCenter[j]    = usgsData[j, 7] 
	halfDurationCenter[j] = usgsData[j, 8] * 0.5

depthRectCenter = - depthRectCenter


fig = plt.figure( 1 )
ax = Axes3D( fig )
ax.scatter( lonRectCenter, latRectCenter, depthRectCenter ) 

if ( src_npts != lengthPoint * depthPoint ):
	print( "src_npts = %d, lengthPoint = %d, depthPoint = %d" %( src_npts, lengthPoint, depthPoint ) )
	print( "The lengthPoint and depthPoint you set is wrong"  )
	sys.exit( )


nx = lengthPoint
nz = depthPoint
npts = nz * nx



X = np.zeros( npts )
Y = np.zeros( npts )
Z = np.zeros( npts )

lon  = np.zeros( npts )
lat  = np.zeros( npts )
depth = np.zeros( npts )



slip = np.zeros( npts )
timeField = np.zeros( npts )
halfDuration = np.zeros( npts )
Area = np.zeros( npts )
rake = np.zeros( npts  )

lon0 = ( np.max( lonRectCenter ) + np.min( lonRectCenter ) )* 0.5 
lat0 = ( np.max( latRectCenter ) + np.min( latRectCenter ) )* 0.5 
print( "Longitude and latitude center: lon0 = %f, lat = %f" % ( lon0, lat0 ) )
proj = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84" )

XCenter, YCenter = proj( lonRectCenter, latRectCenter )
ZCenter = depthRectCenter

xTheta = 0#- 45 * np.pi / 180
yTheta = -dip * np.pi / 180#45 * np.pi / 180
#yTheta = 0#45 * np.pi / 180
zTheta = strike * np.pi / 180

XCenter, YCenter, ZCenter = rotateZ( zTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateX( xTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateY( yTheta, XCenter, YCenter, ZCenter )


disXmove = 0.5 * Dx#( XCenter[lengthPoint] - XCenter[0] ) * 0.5
disYmove = 0.5 * Dy#YCenter[1] - YCenter[0] ) * 0.5


'''
for j in range( nz - 1 ):
	for i in range( nx - 1 ):
		index = i + j * nx
		pos = i + j * lengthPoint
		X[index] = XCenter[pos] - disXmove
		Y[index] = YCenter[pos] - disYmove
		


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
'''
X = XCenter
Y = YCenter


interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), slipCenter )
slip = interp( X, Y )
interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), halfDurationCenter )
halfDuration = interp( X, Y )
interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), timeFieldCenter )
timeField = interp( X, Y )
interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), rakeCenter )
rake = interp( X, Y )

slip = slipCenter
halfDuration = halfDurationCenter
timeField = timeFieldCenter
rake = rakeCenter



M0 = 0.0
miu = 3e10
sumSlip = np.sum( slip )
M0 = sumSlip * Dx * Dy  * miu
Mw = 2/3 * np.log10( M0 ) - 6.06
print( "The M0 is %e"% M0  )
print( "The Mw is %e"% Mw  )


fig = plt.figure( 2 )
ax = Axes3D( fig )
ax.scatter( X, Y, Z ) 


NX = ( nx - 1 ) * thickenX + 1
NY = ( nz - 1 ) * thickenZ + 1
NPTS = NX * NY
XMesh = np.zeros( NPTS )
YMesh = np.zeros( NPTS )
ZMesh = np.zeros( NPTS )

SlipMesh = np.zeros( NPTS )
HalfDurationMesh = np.zeros( NPTS )
TimeFieldMesh = np.zeros( NPTS )
AreaMesh = np.zeros( NPTS )
rakeMesh = np.zeros( NPTS )
SlipRateMesh = np.zeros( [NPTS, NT] )
RakeTimeMesh = np.zeros( [NPTS, NT] )
StrikeMesh = np.zeros( NPTS ) + strike
DipMesh = np.zeros( NPTS  ) + dip




XCoord = np.zeros( [NY, NX]  )
YCoord = np.zeros( [NY, NX]  )
ZCoord = np.zeros( [NY, NX]  )

LonCoord = np.zeros( [NY, NX]  )
LatCoord = np.zeros( [NY, NX]  )

data = np.zeros( [NY, NX] )
dataTime = np.zeros( [NY, NX] )
dataHalfDuration = np.zeros( [NY, NX] )
dataSlip = np.zeros( [NY, NX] )
dataArea = np.zeros( [NY, NX] )
dataRake = np.zeros( [NY, NX] )
dataSlipRate = np.zeros( [NY, NX, NT] )



thickenMesh( nx, nz, thickenX, thickenZ, X, Y, Z, XMesh, YMesh, ZMesh )




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

interp = CloughTocher2DInterpolator( list(zip(X, Y)), slip)
SlipMesh = interp( XMesh, YMesh )

for i in range( NPTS ):
	if SlipMesh[i] < 0:
		SlipMesh[i] = 0.0

interp = LinearNDInterpolator( list(zip(X, Y)), halfDuration)
HalfDurationMesh = interp( XMesh, YMesh )
interp = LinearNDInterpolator( list(zip(X, Y)), timeField)
TimeFieldMesh = interp( XMesh, YMesh )
interp = LinearNDInterpolator( list(zip(X, Y)), rake )
rakeMesh = interp( XMesh, YMesh )

print( "AreaMesh = %f" % AreaMesh[0]  )


M0 = 0.0
miu = 3e10
#sumSlipMesh = np.sum( SlipMesh )
for i in range( NPTS ):
	M0 += SlipMesh[i] * AreaMesh[i] * miu
#M0 = sumSlipMesh * AreaMesh[0] * miu
Mw = 2/3 * np.log10( M0 ) - 6.06
print( "The M0 is %e"% M0  )
print( "The Mw is %e"% Mw  )

frontT = np.max( HalfDurationMesh ) * 4
backT = np.max( HalfDurationMesh ) * 4
dt = ( np.max( TimeFieldMesh ) + frontT + backT ) / NT


t = np.linspace( 0, NT * dt, NT )
tmp = np.zeros( NT ) + 1.0

plt.figure( 15 )
for i in range( NPTS ):
	source = np.zeros( NT )
	t0 = TimeFieldMesh[i] + frontT
	'''
	source = gaussFunc( t, t0, HalfDurationMesh[i] )#
	'''
	r = HalfDurationMesh[i] * 2
	startT = int(t0 / dt)
	endT = int((t0 + r)/dt )
	nt = endT - startT
	soooo = sourceTimeFunc( r, nt )
	source[startT:startT+nt] = soooo
	if r < 1e-5:
		source[:] = 0.0
	SlipRateMesh[i, :] = source * SlipMesh[i]
	#plt.plot( SlipRateMesh[i, :] )
	RakeTimeMesh[i, :] = rakeMesh[i] * tmp


for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XCoord[j, i] = XMesh[index]
		YCoord[j, i] = YMesh[index]
		ZCoord[j, i] = ZMesh[index]
		data[j, i] = SlipMesh[index]
		dataTime[j, i] = TimeFieldMesh[index]
		dataHalfDuration[j, i] = HalfDurationMesh[index]
		dataSlip[j, i] = SlipMesh[index]
		dataArea[j, i] = AreaMesh[index]
		dataRake[j, i] = rakeMesh[index]
		dataSlipRate[j, i, :] = SlipRateMesh[index, :]



XMesh, YMesh, ZMesh = rotateY( -yTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateX( -xTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateZ( -zTheta, XMesh, YMesh, ZMesh )

lonMesh = np.zeros( NPTS )
latMesh = np.zeros( NPTS )

proj1 = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84")
lonMesh, latMesh = proj1( XMesh, YMesh, inverse = "true" )



jsonsFile = open( "params.json" )
params = js.load( jsonsFile )

sourceFileName = params["sourceDir"] + "/USGSsource.bin"
sourceFile = open( sourceFileName, "wb" )

value = struct.pack( "i", NPTS )
sourceFile.write( value )

value = struct.pack( "i", NT )
sourceFile.write( value )

value = struct.pack( "f", dt )
sourceFile.write( value )

for i in range( NPTS ):
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



for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		LonCoord[j, i] = lonMesh[index]
		LatCoord[j, i] = latMesh[index]
		

fig = plt.figure( 3 )
plt.pcolormesh( LonCoord, LatCoord, dataSlip, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

plt.savefig( "fault.png" )


for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		LonCoord[j, i] = lonMesh[index]
		LatCoord[j, i] = latMesh[index]
		

fig = plt.figure( 4 )
plt.pcolormesh( XCoord, YCoord, dataSlip, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

plt.savefig( "fault.png" )


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


