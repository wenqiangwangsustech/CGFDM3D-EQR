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


def triangleFunc( t, t0, a, dt ):
	h = 1 / a
	startT = int( ( t0 - a ) / dt )
	centerT = int( t0 / dt )
	endT   = int( ( t0 + a ) / dt )
	source = np.zeros( len( t ) )
	source[startT : centerT] = h / a * ( t[startT:centerT] - t[startT] )# - h * h * np.abs( t[startT : endT] - t0 ) 
	source[centerT : endT ] = h / a * ( t[endT] - t[centerT:endT] )#h - h * h * np.abs( t[startT : endT] - t0 ) 
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



NT = 1000

lengthPoint = 17
depthPoint = 9

thickenX = 20
thickenZ = 20

file1 = open( "./wangweimindata/slipsub.dat", "r" )
file2 = open( "./wangweimindata/slipsub.par", "r" )

src_npts = lengthPoint * depthPoint

wwm_lon = np.zeros( src_npts )
wwm_lat = np.zeros( src_npts )
wwm_depth = np.zeros( src_npts )
wwm_slip = np.zeros( src_npts )
wwm_dip = np.zeros( src_npts )
wwm_rake = np.zeros( src_npts )
wwm_strike  = np.zeros( src_npts )
wwm_timeField = np.zeros( src_npts )
wwm_halfDuration = np.zeros( src_npts )

wwm_Dx = np.zeros( src_npts )
wwm_Dy = np.zeros( src_npts )



lines1 = [ ]
numLine = 0
while True:
	line = file1.readline( )
	if not line or line.isspace( ):
		break	
	numLine += 1
	lines1.append( line )


lines2 = [ ]
numLine = 0
while True:
	line = file2.readline( )
	if not line or line.isspace( ):
		break	
	numLine += 1
	lines2.append( line )


for i in range( src_npts ):
	strLine = lines1[i].split( )
	wwm_lon[i]       	= float( strLine[0] )
	wwm_lat[i]       	= float( strLine[1] )
	wwm_depth[i]     	= - np.abs( float( strLine[2] ) ) * 1000
	wwm_slip[i]      	= float( strLine[3] )*0.01
	wwm_rake[i]      	= float( strLine[4] )
	wwm_timeField[i] 	= float( strLine[5] )
	wwm_halfDuration[i] = float( strLine[6] ) / 1.628




for i in range( src_npts ):
	strLine = lines2[i].split( )
	wwm_Dx[i]   = float( strLine[3] ) * 1000
	wwm_Dy[i]   = float( strLine[4] ) * 1000
	wwm_dip[i] = float( strLine[5] )
	wwm_strike[i]  = float( strLine[6] )

wwm2_lon			= np.zeros( src_npts )
wwm2_lat			= np.zeros( src_npts )
wwm2_depth			= np.zeros( src_npts )
wwm2_slip			= np.zeros( src_npts )
wwm2_dip			= np.zeros( src_npts )
wwm2_rake			= np.zeros( src_npts )
wwm2_strike			= np.zeros( src_npts )
wwm2_timeField		= np.zeros( src_npts )
wwm2_halfDuration	= np.zeros( src_npts )
wwm2_Dx				= np.zeros( src_npts )
wwm2_Dy				= np.zeros( src_npts )


for j in range( depthPoint ):
	for i in range( lengthPoint ):
		index = i + j * lengthPoint
		pos =  j + i * depthPoint
		wwm2_lon			[index] = wwm_lon			[pos] 
		wwm2_lat			[index] = wwm_lat			[pos] 
		wwm2_depth			[index] = wwm_depth			[pos] 
		wwm2_slip			[index] = wwm_slip			[pos] 
		wwm2_dip			[index] = wwm_dip			[pos] 
		wwm2_rake			[index] = wwm_rake			[pos] 
		wwm2_strike			[index] = wwm_strike		[pos] 
		wwm2_timeField		[index] = wwm_timeField		[pos] 
		wwm2_halfDuration	[index] = wwm_halfDuration	[pos] 
		wwm2_Dx				[index] = wwm_Dx			[pos] 
		wwm2_Dy				[index] = wwm_Dy			[pos] 







strike = wwm2_strike[0]
dip = wwm2_dip[0]


lon0 = ( np.max( wwm2_lon ) + np.min( wwm2_lon ) )* 0.5
lat0 = ( np.max( wwm2_lat ) + np.min( wwm2_lat ) )* 0.5

proj = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84" )
XCenter, YCenter = proj( wwm2_lon, wwm2_lat )
ZCenter = wwm2_depth 




nx = lengthPoint + 1
nz = depthPoint + 1
npts = nz * nx

lon  = np.zeros( npts )
lat  = np.zeros( npts )
depth = np.zeros( npts )


X = np.zeros( npts )
Y = np.zeros( npts )
Z =  np.zeros( npts )

slip = np.zeros( npts  )
slipRate = np.zeros( [ npts, NT]  )
rake = np.zeros( npts  )
half_duration = np.zeros( npts  )
timeField = np.zeros( npts )

xTheta = 0#- 45 * np.pi / 180
yTheta = -dip * np.pi / 180#45 * np.pi / 180
#yTheta = 0#45 * np.pi / 180
zTheta = strike * np.pi / 180

XCenter, YCenter, ZCenter = rotateZ( zTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateX( xTheta, XCenter, YCenter, ZCenter )
XCenter, YCenter, ZCenter = rotateY( yTheta, XCenter, YCenter, ZCenter )



disXmove = 0.5 * wwm2_Dx[0]#( XCenter[lengthPoint] - XCenter[0] ) * 0.5
disYmove = 0.5 * wwm2_Dy[0]#( YCenter[1] - YCenter[0] ) * 0.5



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


interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), ZCenter )
Z = interp( X, Y )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), wwm2_slip )
slip = interp( X, Y )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), wwm2_rake )
rake = interp( X, Y )


interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), wwm2_halfDuration )
halfDuration = interp( X, Y )

interp = NearestNDInterpolator( list(zip(XCenter, YCenter)), wwm2_timeField )
timeField = interp( X, Y )




NX = ( nx - 1 ) * thickenX + 1
NY = ( nz - 1 ) * thickenZ + 1
NPTS = NX * NY
XMesh = np.zeros( NPTS )
YMesh = np.zeros( NPTS )
ZMesh = np.zeros( NPTS )

SlipMesh = np.zeros( NPTS )
AreaMesh = np.zeros( NPTS )
rakeMesh = np.zeros( NPTS )
SlipRateMesh = np.zeros( [NPTS, NT] )
IntegaralSlipMesh = np.zeros( [NPTS] )
RakeTimeMesh = np.zeros( [NPTS, NT] )
StrikeMesh = np.zeros( NPTS ) + strike
DipMesh = np.zeros( NPTS  ) + dip
HalfDurationMesh = np.zeros( NPTS )
TimeFieldMesh = np.zeros( NPTS )


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
for i in range( NY * NX ):
	if SlipMesh[i] < 0. :
		SlipMesh[i] = 0.0

interp = LinearNDInterpolator( list(zip(X, Y)), halfDuration)
HalfDurationMesh = interp( XMesh, YMesh )
interp = CloughTocher2DInterpolator( list(zip(X, Y)), timeField)
TimeFieldMesh = interp( XMesh, YMesh )
interp = NearestNDInterpolator( list(zip(X, Y)), rake )
rakeMesh = interp( XMesh, YMesh )

frontT = np.max( HalfDurationMesh ) * 4
backT = np.max( HalfDurationMesh ) * 4
dt = ( np.max( TimeFieldMesh ) + frontT + backT ) / NT


t = np.linspace( 0, NT * dt, NT )
tmp = np.zeros( NT ) + 1.0



for i in range( NPTS ):
	t0 = TimeFieldMesh[i] + frontT
	source = gaussFunc( t, t0, HalfDurationMesh[i] )#
	#source = triangleFunc( t, t0, HalfDurationMesh[i], dt )#
	SlipRateMesh[i, :] = source * SlipMesh[i]
	RakeTimeMesh[i, :] = rakeMesh[i] * tmp



for i in range( NY * NX ):
	for it in range( NT ):
		RakeTimeMesh[i, it] = rakeMesh[i]

M = 0.0
miu = 3e10
for i in range( NPTS ):
	M += SlipMesh[i] * AreaMesh[i] * miu


Mw = 2/3 * np.log10( M ) - 6.06
print( "M = %e, Mw = %f" % ( M, Mw ) )


for i in range( NY * NX ):
	sum0 = 0.
	for it in range( 1, NT ):
		sum0 += ( SlipRateMesh[i, it - 1] + SlipRateMesh[i, it - 1] ) * dt * 0.5
	IntegaralSlipMesh[i] = sum0 
	

XMesh, YMesh, ZMesh = rotateY( -yTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateX( -xTheta, XMesh, YMesh, ZMesh )
XMesh, YMesh, ZMesh = rotateZ( -zTheta, XMesh, YMesh, ZMesh )



lonMesh = np.zeros( NPTS )
latMesh = np.zeros( NPTS )

proj1 = Proj( proj='aeqd', lat_0 = lat0, lon_0 = lon0, ellps = "WGS84")
lonMesh, latMesh = proj1( XMesh, YMesh, inverse = "true" )



jsonsFile = open( "params.json" )
params = js.load( jsonsFile )

sourceFileName = params["sourceDir"] + "/wangweimin_source.bin"
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



'''
XCoord = np.zeros( [NY, NX] )
YCoord = np.zeros( [NY, NX] )
ZCoord = np.zeros( [NY, NX] )
dataHalfDuration = np.zeros( [NY, NX] )

for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XCoord[j, i] = XMesh[index]
		YCoord[j, i] = YMesh[index]
		ZCoord[j, i] = ZMesh[index]
		dataHalfDuration[j, i] = HalfDurationMesh[index]
		

fig = plt.figure( 12 )
plt.pcolormesh( XCoord, YCoord, dataHalfDuration, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

data = np.zeros( [NY, NX] )
dataSlip = np.zeros( [NY, NX] )
dataArea = np.zeros( [NY, NX] )
dataRake = np.zeros( [NY, NX] )
dataStrike = np.zeros( [NY, NX] )
dataDip = np.zeros( [NY, NX] )
dataSlipRate = np.zeros( [NY, NX, NT] )

for j in range( NY ):
	for i in range( NX ):
		index = i + j * NX
		XCoord[j, i] = lonMesh[index]
		YCoord[j, i] = latMesh[index]
		ZCoord[j, i] = ZMesh[index]
		dataSlip[j, i] = SlipMesh[index]
		#dataSlip[j, i] = IntegaralSlipMesh[index]
		#dataArea[j, i] = AreaMesh[index]
		#dataStrike[j, i] = StrikeMesh[index]
		#dataRake[j, i] = rakeMesh[index]
		#dataSlipRate[j, i, :] = SlipRateMesh[index, :]
		

fig = plt.figure( 13)
plt.pcolormesh( XCoord, YCoord, dataSlip, cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )

fig = plt.figure( 16)
plt.pcolormesh( np.flipud( dataSlip ), cmap = "seismic" )
plt.colorbar( )
plt.axis( "equal" )
'''
