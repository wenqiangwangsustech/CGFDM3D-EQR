#!/usr/bin/env python
'''
Wenqiang Wang@SUSTech 

Date: 2021/05/29

'''


import numpy as np 
import matplotlib.pyplot as plt
import math

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

if __name__ == "__main__":
	x = 1
	y = 1
	z = 1
	
	xTheta = 0#- 45 * np.pi / 180
	yTheta = 0#45 * np.pi / 180
	zTheta = 45 * np.pi / 180
	
	x, y, z = rotateX( xTheta, x, y, z )
	x, y, z = rotateY( yTheta, x, y, z )
	#print( "x = %f, y = %f, z = %f" % ( x, y, z ) )
	
	x, y, z = rotateZ( zTheta, x, y, z )
	print( "x = %f, y = %f, z = %f" % ( x, y, z ) )



#right-hand rule




