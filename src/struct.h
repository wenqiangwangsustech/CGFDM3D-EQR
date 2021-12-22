/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:struct.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-04
*   Discription:
*
================================================================*/
#ifndef __STRUCT__
#define __STRUCT__
typedef struct PARAMS
{
	
    double TMAX;
    
	double DT;
    double DH;
    
	int NX;
    int NY;
    int NZ;

    int PX;
    int PY;
    int PZ;

	int centerX;
	int centerY;

	double centerLatitude ; 
	double centerLongitude;
    
	int sourceX; 
	int sourceY;
	int sourceZ;

	int IT_SKIP;
	
	int sliceX; 
	int sliceY;
	int sliceZ;
	int sliceFreeSurf;


	int nPML;

	int gauss_hill;
	int useTerrain;
	int useMedium;
	int useMultiSource;
	int useSingleSource;
	float rickerfc;

	int ShenModel;
	int Crust_1Medel;

	
	int itSlice			;
	int itStep			;
	char waveOutput[64]	;
	char sliceName[64]	;
	int itStart			;
	int itEnd			;
	int igpu			;
    char OUT[256];
    
	
	char TerrainDir[256];
	
	
	int SRTM90;
	
	int lonStart;
	int latStart;
	int blockX;
	int blockY;
		
	float Depth;


	float MLonStart;
	float MLatStart;
	float MLonEnd;
	float MLatEnd;

	float MLonStep;
	float MLatStep;


	float CrustLonStep;
	float CrustLatStep;

	float MVeticalStep;
	
	char MediumDir[256];
	char crustDir[256];

	char sourceFile[256];
	char sourceDir[256];


}PARAMS;


typedef struct COORD
{	
	float * x;
	float * y;
	float * z;
}COORDINATE, COORD;


typedef struct LONLAT
{
	double * lon;
	double * lat;
	double * depth;
}LONLAT;

typedef struct GRID
{

	int PX;
	int PY;
	int PZ;

	int _NX_;
	int _NY_;
	int _NZ_;

	int _NX;
	int _NY;
	int _NZ;

	int NX;
	int NY;
	int NZ;

	int _nx_;
	int _ny_;
	int _nz_;

	int _nx;
	int _ny;
	int _nz;

	int nx;
	int ny;
	int nz;

	int frontNX;
	int frontNY;
	int frontNZ;

	int _frontNX;
	int _frontNY;
	int _frontNZ;

	int originalX;
	int originalY;

	int _originalX;
	int _originalY;
	//int originalZ;

	int halo;


	int nPML;
	
	float DH;
	float rDH;

}GRID;


typedef struct MPI_NEIGHBOR
{
	int X1; //left
	int X2; //right

	int Y1; //front
	int Y2; //back

	int Z1; //down
	int Z2; //up
	
}MPI_NEIGHBOR;


typedef struct NCFILE
{
	int ncID;

	int ntDimID;	
	int nzDimID;	
	int nyDimID;	
	int nxDimID;
	
	int VxVarID;
	int VyVarID;
	int VzVarID;
	
	int coordXVarID;
	int coordYVarID;
	int coordZVarID;
	
	
	int lonVarID;
	int latVarID;


}NCFILE;

typedef struct SOURCE_FILE_INPUT
{
	
	long long npts; //source point number
	int nt;   //number of source time sequences of every point
	float dt; //time sample interval

	float * lon;
	float * lat;
	float * coordZ;

	float * area;
	float * strike;
	float * dip;

	float * rake;
	float * rate;


}SOURCE_FILE_INPUT;


typedef struct POINT_INDEX
{
	int X;
	int Y;
	int Z;

}POINT_INDEX;

typedef struct STATION_INDEX
{
	int X;
	int Y;
	int Z;

}STATION_INDEX, STA_IDX, STATION_IDX;


typedef struct WAVE
{
	float * Vx; 
	float * Vy; 
	float * Vz; 
	float * Txx;
	float * Tyy;
	float * Tzz;
	float * Txy;
	float * Txz;
	float * Tyz;
}WAVE;

typedef struct WAVE4VAR
{
	WAVE h_W;
	WAVE W;
	WAVE t_W;
	WAVE m_W;

}WAVE4VAR, WAVE4;


typedef struct AUXILIARY
{
	float * Vx; 
	float * Vy; 
	float * Vz; 
	float * Txx;
	float * Tyy;
	float * Tzz;
	float * Txy;
	float * Txz;
	float * Tyz;
}AUXILIARY, AUX;

typedef struct AUXILIARY4VAR
{
	AUX h_Aux_x;
	AUX   Aux_x;
	AUX t_Aux_x;
	AUX m_Aux_x;

	AUX h_Aux_y;
	AUX   Aux_y;
	AUX t_Aux_y;
	AUX m_Aux_y;

	AUX h_Aux_z;
	AUX   Aux_z;
	AUX t_Aux_z;
	AUX m_Aux_z;

}AUXILARY4VAR, AUXILIARY4, AUX4;

typedef struct PML_ALPHA
{
	float * x; 
	float * y; 
	float * z; 
}PML_ALPHA;

typedef struct PML_BETA
{
	float * x; 
	float * y; 
	float * z; 
}PML_BETA;

typedef struct PML_D
{
	float * x;	
	float * y;	
	float * z;	
}PML_D;


typedef struct MPI_BORDER
{
	int isx1; int isx2;
	int isy1; int isy2;
	int isz1; int isz2;

}MPI_BORDER;

typedef struct CONTRAVARIANT
{
	float * xi_x;
	float * xi_y;
	float * xi_z;
	float * et_x;
	float * et_y;
	float * et_z;
	float * zt_x;
	float * zt_y;
	float * zt_z;
}CONTRAVARIANT;

typedef struct Mat3x3
{
	float * M11; float * M12; float * M13;
	float * M21; float * M22; float * M23;
	float * M31; float * M32; float * M33;
}Mat3x3;

typedef struct SLICE
{
	int X;
	int Y;
	int Z;
}SLICE;

typedef struct SLICE_DATA
{
	float * x;
	float * y;
	float * z;

}SLICE_DATA;

typedef struct SOURCE
{
	int X;
	int Y;
	int Z;
}SOURCE;

typedef struct SOURCE_INDEX
{
	int * X;
	int * Y;
	int * Z;
}SOURCE_INDEX;

typedef struct MEDIUM
{
	float * mu;
	float * lambda;
	float * buoyancy;
}MEDIUM;

typedef struct STRUCTURE
{
	float * Vs;
	float * Vp;
	float * rho;
}STRUCTURE;

 
typedef struct SEND_RECV_DATA
{
	float * thisXSend1;
	float * thisXRecv1;
	float * thisYSend1;
	float * thisYRecv1;
	float * thisZSend1;
	float * thisZRecv1;

	float * thisXSend2;
	float * thisXRecv2;
	float * thisYSend2;
	float * thisYRecv2;
	float * thisZSend2;
	float * thisZRecv2;

}SEND_RECV_DATA;

typedef struct  WSLICE{
	float * sliceX;
	float * sliceY;
	float * sliceZ;
}WSLICE;  


typedef struct MPI_COORDINATE
{
	int X;
	int Y;
	int Z;
}MPI_COORDINATE, MPI_COORD;


typedef struct DELTA_H_RANGE
{
	float * DT_min;
	float * DT_max;
}DELTA_H_RANGE;

typedef struct POINT_OR_VECTOR
{
	double x;
	double y;
	double z;
}POINT_OR_VECTOR;


typedef struct SOURCE_INFO
{	
	int npts;
	int nt;
	float dt;
}SOURCE_INFO;


typedef struct MOMENT_RATE
{	
	float * Mxx;
	float * Myy;
	float * Mzz;
	float * Mxy;
	float * Mxz;
	float * Myz;

}MOMENT_RATE;

typedef struct PGV
{
	float * pgvh;
	float * pgv;
}PGV;

typedef void (*WAVE_RK_FUNC )( float * h_W, float * W, float * t_W, float * m_W, long long num, float DT );

#endif //__STRUCT__
