#/usr/bin/env python
'''
Author: Wenqiang Wang @ SUSTech on Sep 11, 2021
15:40
'''
import json
import numpy as np

class GRID:
	def __init__( self, params, HALO = 3 ):
		resX = 0 	
		resY = 0 
		resZ = 0 	
			
		self.PX = params["PX"] 
		self.PY = params["PY"] 
		self.PZ = params["PZ"] 

		self._NX_ = params["NX"] + 2 * HALO 
		self._NY_ = params["NY"] + 2 * HALO 
		self._NZ_ = params["NZ"] + 2 * HALO 

		self._NX = params["NX"] + HALO 
		self._NY = params["NY"] + HALO 
		self._NZ = params["NZ"] + HALO 

		self.NX = params["NX"] 
		self.NY = params["NY"] 
		self.NZ = params["NZ"] 

		resX = self.NX % self.PX 
		resY = self.NY % self.PY 
		resZ = self.NZ % self.PZ

		self.nx = np.zeros( self.PX, dtype = 'int32' ) + self.NX // self.PX 
		self.ny = np.zeros( self.PY, dtype = 'int32' ) + self.NY // self.PY 
		self.nz = np.zeros( self.PZ, dtype = 'int32' ) + self.NZ // self.PZ 
		
		self.frontNX = np.zeros( self.PX, dtype = 'int32' )
		self.frontNY = np.zeros( self.PY, dtype = 'int32' )
		self.frontNZ = np.zeros( self.PZ, dtype = 'int32' )
		


		for mpiX in range( self.PX ):
			if ( mpiX < resX ):
				self.nx[mpiX] += 1 
				self.frontNX[mpiX] = mpiX * self.nx[mpiX] 
			else:
				self.frontNX[mpiX] = resX * ( self.nx[mpiX] + 1 ) + ( mpiX - resX ) * self.nx[mpiX] 


		for mpiY in range( self.PY ):
			if ( mpiY < resY ):
				self.ny[mpiY] += 1 
				self.frontNY[mpiY] = mpiY * self.ny[mpiY] 
			else:
				self.frontNY[mpiY] = resY * ( self.ny[mpiY] + 1 ) + ( mpiY - resY ) * self.ny[mpiY] 

		for mpiZ in range( self.PZ ):
			if ( mpiZ < resZ ):
				self.nz[mpiZ] += 1 
				self.frontNZ[mpiZ] = mpiZ * self.nz[mpiZ] 
			else:
				self.frontNZ[mpiZ] = resZ * ( self.nz[mpiZ] + 1 ) + ( mpiZ - resZ ) * self.nz[mpiZ] 


		self._frontNX = self.frontNX + HALO 
		self._frontNY = self.frontNY + HALO 
		self._frontNZ = self.frontNZ + HALO 


		self._nx = self.nx + HALO 
		self._ny = self.ny + HALO 
		self._nz = self.nz + HALO 

		self._nx_ = self.nx + 2 * HALO 
		self._ny_ = self.ny + 2 * HALO 
		self._nz_ = self.nz + 2 * HALO 

		self.originalX = params["centerX"] 
		self.originalY = params["centerY"] 

		self._originalX = self.originalX + HALO 
		self._originalY = self.originalY + HALO 

		self.halo = HALO 

		self.DH = params["DH"] 

def main( ):
	jsonsFile = open( "params.json" )
	params = json.load( jsonsFile )
	grid = GRID( params )
	print( grid )


if __name__ == '__main__':
	main()




