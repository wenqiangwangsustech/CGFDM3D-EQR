#!/usr/bin/env python
'''
Wenqiang Wang@SUSTech 

Date: 2023/3/13

'''


import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import struct
import json as js
from pyproj import Proj

import math

import scipy as sp
import scipy.ndimage

import matplotlib.animation as animation

from scipy.interpolate import NearestNDInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import CloughTocher2DInterpolator

def gaussFunc( t, t0, a ):
    source = np.exp(-( t - t0 )**2 / ( a ** 2 ) ) / ( np.sqrt( np.pi ) * a )
    return source

def trifunc(t, t0, sig, dt):
    '''
    Triangle source time function
    '''
    stf = np.zeros(len(t))
    startT = int(( t0 - sig ) / dt)  # index
    endT   = int(( t0 + sig ) / dt)
    stf[startT : endT] = 1.0/sig - (1.0/sig**2) * np.abs(t[startT : endT] - t0)
    return stf


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
	


def lonlatdepth2XYZ( lon, lat, depth, lonC, latC ):
	#lonC = ( np.max( lon ) + np.min( lon ) ) / 2
	#latC = ( np.max( lat ) + np.min( lat ) ) / 2
	print( "lonC = %f, latC = %f"%( lonC, latC ) )

	proj = Proj( proj='aeqd', lat_0 = latC, lon_0 = lonC, ellps = "WGS84" )
	X, Y = proj( lon, lat )
	Z = - depth

	return X, Y, Z

def triangle_area_cross(x1, y1, x2, y2, x3, y3):
	return 0.5 * np.abs(x1 * y2 + x2 * y3 + x3 * y1 - x2 * y1 - x3 * y2 - x1 * y3)


def MergeFault( X, Y, Z, strike, dip, duration, mu, s_slip, d_slip, s_slip_total, d_slip_total, timeField, nx, ny, thicken, NF, Num ):
	averStrike = np.average( strike )
	averDip = np.average( dip )
	sig = 3
	xTheta = 0#- 45 * np.pi / 180
	yTheta = -averDip * np.pi / 180#45 * np.pi / 180
	#yTheta = 0#45 * np.pi / 180
	zTheta = averStrike * np.pi / 180

	slip = np.zeros_like( s_slip )

	X, Y, Z = rotateZ( zTheta, X, Y, Z )
	X, Y, Z = rotateX( xTheta, X, Y, Z )
	X, Y, Z = rotateY( yTheta, X, Y, Z )

	Xmax = np.max( X )
	#print( "Xmax = %f" % Xmax )
	Ymax = np.max( Y )# + 1000
	#print( "Ymax = %f" % Ymax )
	#if NF == 1:
	#	Ymax += 1500

	Xmin = np.min( X )
	#if NF == 1:
	#	Xmin = -500
	#	Xmax -= 1000
	#print( "Xmin = %f" % Xmin )
	Ymin = np.min( Y )

	xcoord = np.zeros( [ny*nx] )
	ycoord = np.zeros( [ny*nx] )
	zcoord = np.zeros( [ny*nx] )

	strikeInterp = np.zeros( [ny*nx] )
	dipInterp    = np.zeros( [ny*nx] )
	durationInterp = np.zeros( [ny*nx] )
	muInterp 	 = np.zeros( [ny*nx] )

	rake_totalInterp = np.zeros( [ny*nx] )
	s_slip_totalInterp = np.zeros( [ny*nx] )
	d_slip_totalInterp = np.zeros( [ny*nx] )
	s_slipInterp = np.zeros( [N, ny*nx] )
	d_slipInterp = np.zeros( [N, ny*nx] )
	slipInterp   = np.zeros( [N, ny*nx] )	
	timeFieldInterp = np.zeros( [N, ny*nx] )


	x = np.linspace( Xmin, Xmax, nx )
	y = np.linspace( Ymin, Ymax, ny )
	print( "============================="  )
	print( "x1 - x0 = %f" % ( x[1] - x[0] ) )
	print( "============================="  )

	print( "============================="  )
	print( "y1 - y0 = %f" % ( y[1] - y[0] ) )
	print( "============================="  )

	for j in range( ny ):
		for i in range( nx ):
			index = i + j * nx
			xcoord[index] = x[i]	
			ycoord[index] = y[j]	

	interp = NearestNDInterpolator( list(zip(X, Y)), Z)
	zcoord = interp( xcoord, ycoord )
	zcoord = scipy.ndimage.filters.gaussian_filter( zcoord, sig )
	
	interp = NearestNDInterpolator( list(zip(X, Y)), strike)
	strikeInterp = interp( xcoord, ycoord )
	interp = NearestNDInterpolator( list(zip(X, Y)), dip)
	dipInterp    = interp( xcoord, ycoord )
	interp = NearestNDInterpolator( list(zip(X, Y)), duration)
	durationInterp = interp( xcoord, ycoord )
	interp = NearestNDInterpolator( list(zip(X, Y)), mu)
	muInterp 	 = interp( xcoord, ycoord )

	rake_total = np.arctan2( d_slip_total, s_slip_total ) * 180 / np.pi

	interp = NearestNDInterpolator( list(zip(X, Y)), rake_total)
	rake_totalInterp = interp( xcoord, ycoord )

	interp = NearestNDInterpolator( list(zip(X, Y)), s_slip_total)
	s_slip_totalInterp = interp( xcoord, ycoord )
	interp = NearestNDInterpolator( list(zip(X, Y)), d_slip_total)
	d_slip_totalInterp = interp( xcoord, ycoord )

	#print( "===========================" )
	#print( np.shape( s_slip[0,:] ) )
	#print( "***************************" )
	for n in range( Num ):
		slip[n, :] = np.sqrt( s_slip[n, :]**2 + d_slip[n, :]**2 )
		interp = NearestNDInterpolator( list(zip(X, Y)), s_slip[n,:])
		s_slipInterp[n,:] = interp( xcoord, ycoord )
		interp = NearestNDInterpolator( list(zip(X, Y)), d_slip[n,:])
		d_slipInterp[n,:]  = interp( xcoord, ycoord )
		interp = NearestNDInterpolator( list(zip(X, Y)), slip[n,:])
		slipInterp[n,:]  = interp( xcoord, ycoord )
		interp = NearestNDInterpolator( list(zip(X, Y)), timeField[n,:])
		timeFieldInterp[n,:]  = interp( xcoord, ycoord )


	NX = thicken[0] * nx
	NY = thicken[1] * ny


	XCoord = np.zeros( [NY*NX] )
	YCoord = np.zeros( [NY*NX] )
	ZCoord = np.zeros( [NY*NX] )

	Strike = np.zeros( [NY*NX] )
	Dip    = np.zeros( [NY*NX] )
	Duration = np.zeros( [NY*NX] )
	Mu 	 = np.zeros( [NY*NX] )
	
	Rake_total = np.zeros( [NY*NX] )
	S_slip_total = np.zeros( [NY*NX] )
	D_slip_total = np.zeros( [NY*NX] )
	
	S_slip = np.zeros( [N, NY*NX] )
	D_slip = np.zeros( [N, NY*NX] )
	Slip   = np.zeros( [N, NY*NX] )	
	TimeField = np.zeros( [N, NY*NX] )

	x = np.linspace( Xmin, Xmax, NX )
	y = np.linspace( Ymin, Ymax, NY )


	print( "*****************************"  )
	print( "X1 - X0 = %f" % ( x[1] - x[0] ) )
	print( "*****************************"  )

	print( "*****************************"  )
	print( "Y1 - Y0 = %f" % ( y[1] - y[0] ) )
	print( "*****************************"  )

	for j in range( NY ):
		for i in range( NX ):
			index = i + j * NX
			XCoord[index] = x[i]	
			YCoord[index] = y[j]	


	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), zcoord)
	ZCoord =  interp(XCoord, YCoord)

	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), strikeInterp)
	Strike = interp( XCoord, YCoord )
	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), dipInterp)
	Dip = interp( XCoord, YCoord )
	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), durationInterp)
	Duration = interp( XCoord, YCoord )
	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), muInterp)
	Mu = interp( XCoord, YCoord )
	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), rake_totalInterp)
	Rake_total = interp( XCoord, YCoord )
	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), s_slip_totalInterp)
	S_slip_total = interp( XCoord, YCoord )
	interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), d_slip_totalInterp)
	D_slip_total = interp( XCoord, YCoord )
	for n in range( Num ):
		interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), s_slipInterp[n,:])
		S_slip[n,:] = interp( XCoord, YCoord )
		interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), d_slipInterp[n,:])
		D_slip[n,:]  = interp( XCoord, YCoord )
		interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), slipInterp[n,:])
		Slip[n,:]  = interp( XCoord, YCoord )
		interp = CloughTocher2DInterpolator( list(zip(xcoord, ycoord)), timeFieldInterp[n,:])
		TimeField[n,:]  = interp( XCoord, YCoord )
	
	Area = np.zeros( [NY*NX] )

	for j in range( NY - 1 ):
		for i in range( NX - 1 ):
			index = i + j * NX
			indexA = ( i + 0 ) + ( j + 0 )* NX #A
			indexB = ( i + 1 ) + ( j + 0 )* NX #B
			indexC = ( i + 1 ) + ( j + 1 )* NX #C
			indexD = ( i + 0 ) + ( j + 1 )* NX #D
			'''
			A===B
			|   |
			|   |
			D===C
			'''
			xA = XCoord[indexA]
			xB = XCoord[indexB]
			xC = XCoord[indexC]
			xD = XCoord[indexD]

			yA = YCoord[indexA]
			yB = YCoord[indexB]
			yC = YCoord[indexC]
			yD = YCoord[indexD]
			x1 = xA; x2 = xB; x3 = xC; 
			y1 = yA; y2 = yB; y3 = yC; 
			Area[index] = triangle_area_cross(x1, y1, x2, y2, x3, y3)
			x1 = xA; x2 = xD; x3 = xC; 
			y1 = yA; y2 = yD; y3 = yC; 
			Area[index] += triangle_area_cross(x1, y1, x2, y2, x3, y3)
	
	for j in range( NY ):
		i = NX - 1; index = i + j * NX; 
		i = NX - 2; index2 = i + j * NX; 
		Area[index] = Area[index2]

	for i in range( NX ):
		j = NY - 1; index = i + j * NX; 
		j = NY - 2; index2 = i + j * NX; 
		Area[index] = Area[index2]


	XCoord, YCoord, ZCoord = rotateY( -yTheta, XCoord, YCoord, ZCoord )
	XCoord, YCoord, ZCoord = rotateX( -xTheta, XCoord, YCoord, ZCoord )
	XCoord, YCoord, ZCoord = rotateZ( -zTheta, XCoord, YCoord, ZCoord )

	return XCoord, YCoord, ZCoord, Strike, Dip, Duration, Mu, S_slip, D_slip, Slip, TimeField, S_slip_total, D_slip_total, Rake_total, Area

'''
	XMesh = np.zeros( [NY, NX] )
	YMesh = np.zeros( [NY, NX] )
	ZMesh = np.zeros( [NY, NX] )
	SlipMesh = np.zeros( [NY, NX] )

	for j in range( NY ):
		for i in range( NX ):
			index = i + j * NX
			XMesh[j, i] = XCoord[index]	
			YMesh[j, i] = YCoord[index]	
			ZMesh[j, i] = ZCoord[index]

	XMesh, YMesh, ZMesh = rotateY( -yTheta, XMesh, YMesh, ZMesh )
	XMesh, YMesh, ZMesh = rotateX( -xTheta, XMesh, YMesh, ZMesh )
	XMesh, YMesh, ZMesh = rotateZ( -zTheta, XMesh, YMesh, ZMesh )
'''

def plotGIF( plt, it, fig, XCoord, YCoord, ZCoord, Rate ):
	#for it in range( 0,NT,10 ):
	plt.clf()
	ax = Axes3D( fig )
	handle = ax.scatter( XCoord, YCoord, ZCoord, s = 10, c = Rate[:,it], vmin = 0.0, vmax = 0.2, cmap  = "jet" )
	#ax.scatter( XTemp1, YTemp1, 0.0, s = 10, c = 'k' )
	#ax.scatter( XTemp2, YTemp2, 0.0, s = 10, c = 'k' )
	ax.set_box_aspect((4,4,1)) # 设定坐标轴的长宽高比例
	ax.grid(False) # 不要网格
	plt.colorbar(handle, shrink = 0.5 )
	plt.pause( 0.002 )

def plotData( plt, fig, XCoord, YCoord, ZCoord, data ):
	ax = Axes3D( fig )
	handle = ax.scatter( XCoord, YCoord, ZCoord, s = 10, c = data, cmap  = "jet" )
	ax.set_box_aspect((4,4,1)) # 设定坐标轴的长宽高比例
	ax.grid(False) # 不要网格
	plt.colorbar(handle, shrink = 0.5 )
	plt.pause( 0.002 )

def writeData( sourceFileName, Lon, Lat, Z, Area, Strike, Dip, Rake, Rate ):
	sourceFile = open(sourceFileName, "wb")
	value = struct.pack( "i", NPTS )
	sourceFile.write( value )
	value = struct.pack( "i", NT )
	sourceFile.write( value )
	value = struct.pack( "f", dt )
	sourceFile.write( value )

	for i in range( NPTS ):
		value = struct.pack( "f",  Lon[i] )
		sourceFile.write( value )
		value = struct.pack( "f",  Lat[i] )
		sourceFile.write( value )
		value = struct.pack( "f",  Z[i] )
		sourceFile.write( value )
		value = struct.pack( "f",  Area[i] )
		sourceFile.write( value )
		value = struct.pack( "f",  Strike[i] )
		sourceFile.write( value )
		value = struct.pack( "f",  Dip[i] )
		sourceFile.write( value )
		tvalue = struct.pack( "f" * NT,  *( Rake[i, :] ) )
		sourceFile.write( tvalue )
		tvalue = struct.pack( "f" * NT,  *( Rate[i, :] ) )
		sourceFile.write( tvalue )

	sourceFile.close()



NT = 3000
dt = 0.01
TMAX = NT * dt

nx1 = 9
ny1 = 33

nx2 = 9
ny2 = 8

thicken =[10, 10]

fileTotal = open( "Luding.inv_total", "r" )
file = open( "Luding_staticMulti_str_168_dip_71_rake_-4.0036.inv", "r" )
linesTotal = [ ]
lines = [ ]
numLine = 0
while True:
	line = file.readline( )
	lineTotal = fileTotal.readline( )
	if not line or line.isspace( ):
		break	
	numLine += 1
	lines.append( line )
	linesTotal.append( lineTotal )

src_npts = numLine - 1#399
#src_npts = src_npts5 // N
src_npts = 369
N = 6

NF = 2  #Number of Fault

lon 			= np.zeros( [src_npts] ) 
lat 			= np.zeros( [src_npts] ) 
depth 			= np.zeros( [src_npts] )
strike 			= np.zeros( [src_npts] )
dip				= np.zeros( [src_npts] )
duration		= np.zeros( [src_npts] )
mu  			= np.zeros( [src_npts] )

s_slip_total    = np.zeros( [src_npts] )
d_slip_total    = np.zeros( [src_npts] )

s_slip			= np.zeros( [N, src_npts] )
d_slip			= np.zeros( [N, src_npts] )
timeField		= np.zeros( [N, src_npts] )

for i in range( src_npts ):
	line = lines[i+1].split( )
	lineTotal = linesTotal[i+1].split( )
	lon 		 [i] = float( line[1] )
	lat 		 [i] = float( line[2] )
	depth 		 [i] = float( line[3] ) * 1000
	strike		 [i] = float( line[4] ) 
	dip			 [i] = float( line[5] ) 
	duration	 [i] = float( line[7] ) / 1.628
	mu		     [i] = float( line[13] )
	s_slip_total [i] = float( lineTotal[8] )
	d_slip_total [i] = float( lineTotal[9] )



for n in range( N ):
	for i in range( src_npts ):
		index = i + n * src_npts
		line = lines[index + 1].split( )
		s_slip		 [n,i] = float( line[8] )
		d_slip		 [n,i] = float( line[9] )
		timeField	 [n,i] = float( line[12] )

lonC = ( np.max( lon ) + np.min( lon ) ) / 2
latC = ( np.max( lat ) + np.min( lat ) ) / 2


startNum1 = 0 
endNum1 = 297
X1, Y1, Z1 = lonlatdepth2XYZ( lon[startNum1:endNum1], lat[startNum1:endNum1], depth[startNum1:endNum1], lonC, latC )
strike1 = strike[startNum1:endNum1]
dip1 = dip[startNum1:endNum1] 
duration1 = duration[startNum1:endNum1] 
mu1 = mu[startNum1:endNum1] 
s_slip_total1 = s_slip_total[startNum1:endNum1] 
d_slip_total1 = d_slip_total[startNum1:endNum1] 

startNum2 = 297
endNum2 = 369
X2, Y2, Z2 = lonlatdepth2XYZ( lon[startNum2:endNum2], lat[startNum2:endNum2], depth[startNum2:endNum2], lonC, latC )
strike2 = strike[startNum2:endNum2]
dip2 = dip[startNum2:endNum2] 
duration2 = duration[startNum2:endNum2] 
mu2 = mu[startNum2:endNum2] 
s_slip_total2 = s_slip_total[startNum2:endNum2] 
d_slip_total2 = d_slip_total[startNum2:endNum2] 


s_slip1			= np.zeros( [N, endNum1 - startNum1] )
d_slip1			= np.zeros( [N, endNum1 - startNum1] )
timeField1		= np.zeros( [N, endNum1 - startNum1] )

s_slip2			= np.zeros( [N, endNum2 - startNum2] )
d_slip2			= np.zeros( [N, endNum2 - startNum2] )
timeField2		= np.zeros( [N, endNum2 - startNum2] )

for n in range( N ):
	s_slip1[n,:] = s_slip[n, startNum1:endNum1]
	d_slip1[n,:] = d_slip[n, startNum1:endNum1]
	s_slip2[n,:] = s_slip[n, startNum2:endNum2]
	d_slip2[n,:] = d_slip[n, startNum2:endNum2]
	timeField1[n,:] = timeField[n, startNum1:endNum1] 	
	timeField2[n,:] = timeField[n, startNum2:endNum2] 	


averStrike = np.average( strike )
averDip = np.average( dip )


# fig = plt.figure( 1 )
# ax = Axes3D( fig )
# ax.scatter( X1, Y1, Z1, s = 10 )
# ax.scatter( X2, Y2, Z2, s = 10 )


X, Y, Z = lonlatdepth2XYZ( lon, lat, depth, lonC, latC )
# fig = plt.figure( 2 )
# ax = Axes3D( fig )
# ax.scatter( X, Y, Z, s = 10 )



XCoord1, YCoord1, ZCoord1, Strike1, Dip1, Duration1, Mu1, S_slip1, D_slip1, Slip1, TimeField1, S_slip_total1, D_slip_total1, Rake_total1, Area1 = \
MergeFault( X1, Y1, Z1, strike1, dip1, duration1, mu1, s_slip1, d_slip1, s_slip_total1, d_slip_total1, timeField1, nx1, ny1, thicken, 1, N )


XCoord2, YCoord2, ZCoord2, Strike2, Dip2, Duration2, Mu2, S_slip2, D_slip2, Slip2, TimeField2, S_slip_total2, D_slip_total2, Rake_total2, Area2 = \
MergeFault( X2, Y2, Z2, strike2, dip2, duration2, mu2, s_slip2, d_slip2, s_slip_total2, d_slip_total2, timeField2, nx2, ny2, thicken, 2, N )



NX1 = nx1 * thicken[0]
NY1 = ny1 * thicken[1]
      
NX2 = nx2 * thicken[0]
NY2 = ny2 * thicken[1]

NPTS = NX1 * NY1 + NX2 * NY2

XCoord = np.zeros( [NPTS] )
YCoord = np.zeros( [NPTS] )
ZCoord = np.zeros( [NPTS] )

Strike = np.zeros( [NPTS] )
Dip    = np.zeros( [NPTS] )
Duration = np.zeros( [NPTS] )
Mu = np.zeros( [NPTS] )
S_slip_total = np.zeros( [NPTS] )
D_slip_total = np.zeros( [NPTS] )
Slip_total = np.zeros( [NPTS] )
Rake_total = np.zeros( [NPTS] )
Area  = np.zeros( [NPTS] )

S_slip = np.zeros( [N, NPTS] )
D_slip = np.zeros( [N, NPTS] )
Slip = np.zeros( [N, NPTS] )
TimeField = np.zeros( [N, NPTS] )


for j in range( NY1 ):
	for i in range( NX1 ):
		index = i + j * NX1
		index1 = i + j * NX1
		XCoord[index]    = XCoord1[index1] 	
		YCoord[index]    = YCoord1[index1] 	
		ZCoord[index]    = ZCoord1[index1] 
		Strike[index] 	 = Strike1[index1] 	 
		Dip[index]    	 = Dip1[index1]    	 
		Duration[index]  = Duration1[index1]  
		Mu[index] 		 = Mu1[index1] 		 
		S_slip_total[index] = S_slip_total1[index1]
		D_slip_total[index] = D_slip_total1[index1]
		Rake_total[index] = Rake_total1[index1]
		Area[index] = Area1[index1]
		for n in range( N ): 		 
			S_slip[n, index]    = S_slip1[n, index1]    
			D_slip[n, index]    = D_slip1[n, index1]    
			Slip[n, index]      = Slip1[n, index1]      
			TimeField[n, index] = TimeField1[n, index1]

for j in range( NY2 ):
	for i in range( NX2 ):
		index2 = i + j * NX2
		index = index2 + NX1 * NY1
		XCoord[index]    = XCoord2[index2] 	
		YCoord[index]    = YCoord2[index2] 	
		ZCoord[index]    = ZCoord2[index2] 
		Strike[index] 	 = Strike2[index2] 	 
		Dip[index]    	 = Dip2[index2]    	 
		Duration[index]  = Duration2[index2]  
		Mu[index] 		 = Mu2[index2]
		S_slip_total[index] = S_slip_total2[index2]
		D_slip_total[index] = D_slip_total2[index2]
		Rake_total[index] = Rake_total2[index2]
		Area[index] = Area2[index2]
		for n in range( N ): 		 
			S_slip[n, index]    = S_slip2[n, index2]    
			D_slip[n, index]    = D_slip2[n, index2]    
			Slip[n, index]      = Slip2[n, index2]      
			TimeField[n, index] = TimeField2[n, index2] 


Slip_total = np.sqrt( S_slip_total ** 2 + D_slip_total ** 2 )


# for i in range( NPTS ):
# 	dis = np.abs( XCoord[i] ) ** 2 + np.abs( YCoord[i] ) ** 2  + np.abs( ZCoord[i] ) ** 2 
# 	if dis < 0.1:
# 		print( "==================Index = %d==================" % i )
	
# 	if np.isnan( Slip[0, i] ) == True:
# 		print( "The interpolatation may be false since the Index = %d is nan.\n" % i )



#Rake = np.zeros( [N, NPTS, NT] )
#RateN = [N, NPTS, NT]
Rate = np.zeros( [NPTS, NT] )
Rake = np.zeros( [NPTS, NT] )
t = np.linspace( 0, TMAX, NT )


for i in range( NPTS ):
	#print( np.shape( Rate  ) )
	Rake[i] = Rake_total[i] * np.ones( [NT] )
	for n in range( N ):
		#Rake[n, i, :] = np.ones( [NT] )
		halfDura = Duration[i]/2
		t0 = TimeField[n, i] + halfDura + 1 
		#rate = Slip[n,i] * trifunc( t, t0, halfDura, dt)
		rate = Slip[n,i] * gaussFunc( t, t0, halfDura )
		Rate[i] = Rate[i] + rate


'''
lonlatTrace1 = np.loadtxt("lat_lonF1_smooth.txt")
lonlatTrace2 = np.loadtxt("lat_lonEAF_smooth.txt")

XTemp1, YTemp1, ZTemp1 = lonlatdepth2XYZ( lonlatTrace1[:,0], lonlatTrace1[:,1], 0, lonC, latC )
XTemp2, YTemp2, ZTemp2 = lonlatdepth2XYZ( lonlatTrace2[:,0], lonlatTrace2[:,1], 0, lonC, latC )
'''


fig = plt.figure( 6, figsize = (10, 8) )

#plotData( plt, fig, XCoord/1000, YCoord/1000, ZCoord/1000, Slip_total )
#plotGIF( plt, it, fig, XCoord/1000, YCoord/1000, ZCoord/1000, Rate )


M0 = np.sum( Mu * Area * Slip_total )
Mw = 2.0/3.0 * np.log10(M0) - 6.06
print( "Mw = %f" % Mw )


proj1 = Proj(proj='aeqd', lat_0 = latC, lon_0 = lonC, ellps = "WGS84")
LonCoord, LatCoord = proj1( XCoord, YCoord, inverse = "true")



jsonsFile = open( "params.json" )
params = js.load( jsonsFile )

sourceFileName = params["sourceDir"] + "/LudingXialei.bin"
writeData( sourceFileName, LonCoord, LatCoord, ZCoord, Area, Strike, Dip, Rake, Rate )

#rake = np.arctan2( d_slip_total, s_slip_total )  * 180 / np.pi
#plt.plot( rake )
def update( it ):
	print( "iterator: %d" % it )
	plotGIF( plt, it, fig, XCoord/1000, YCoord/1000, ZCoord/1000, Rate )
	return fig
	
itStart = 0
itEnd = NT
IT_SKIP = 50
anim = animation.FuncAnimation( fig, update, frames = range( itStart, itEnd, IT_SKIP ), interval = 100, repeat = False )
anim.save( "rate" + ".gif", writer = 'pillow', fps = 100 )
'''
'''
