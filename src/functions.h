/*================================================================
*   ESS, Southern University of Science and Technology
*   
*   File Name:functions.h
*   Author: Wenqiang Wang, 11849528@mail.sustech.edu.cn
*   Created Time:2021-09-05
*   Discription:
*
================================================================*/
#ifndef __FUNCTIONS__
#define __FUNCTIONS__

/*extern*/ void getParams( PARAMS * params);

void projTrans( double lon_0, double lat_0, GRID grid, COORD coord, LONLAT LonLat  );
double interp2d(double x[2], double y[2], double z[4], double x_, double y_ );
/*
double bilinear( double x, double y, double x1, double x2, double y1, double y2, double f11, double f12, double f21, double f22 ); 
double bilinearInterp( double x, double y, double x1, double y1, double x2, double y2, double A, double B, double C, double D );
void cart2LonLat( GRID grid, PJ * P,  PJ_DIRECTION PD, COORD coord, LONLAT LonLat );
void preprocessTerrain( PARAMS params, MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, COORD coord );
*/
void init_grid( PARAMS params, GRID * grid, MPI_COORD thisMPICoord );
void printInfo( GRID grid );
void createDir( PARAMS params );

void allocCoord( GRID grid, COORD * coord );
void allocCoord_cpu( GRID grid, COORD * coord );
void freeCoord( COORD coord );
void freeCoord_cpu( COORD coord );
#ifdef GPU_CUDA
void constructCoord(MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, PARAMS params, COORD coord, COORD cpu_coord );
#else                                                           
void constructCoord(MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, PARAMS params, COORD coord );
#endif 

void preprocessTerrain( PARAMS params, MPI_Comm comm_cart, MPI_COORD thisMPICoord, GRID grid, COORD coord );



void allocMedium( GRID grid, MEDIUM * medium );
void freeMedium( MEDIUM medium );
void allocMedium_cpu( GRID grid, MEDIUM * medium );
void freeMedium_cpu( MEDIUM medium );
#ifdef GPU_CUDA
void constructMedium( MPI_COORD thisMPICoord, PARAMS params, GRID grid, COORD coord, MEDIUM medium, MEDIUM cpu_medium );
#else
void constructMedium( MPI_COORD thisMPICoord, PARAMS params, GRID grid, COORD coord, MEDIUM medium );
#endif

void readWeisenShenModel( PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord, STRUCTURE structure );
void readCrustal_1( PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord, STRUCTURE structure );

void calc_CFL( GRID grid, COORD coord, MEDIUM medium, PARAMS params );

void vs_vp_rho2lam_mu_bou( MEDIUM medium, GRID grid );



void init_MultiSource( 
		PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord, 
		long long ** srcIndex, MOMENT_RATE * momentRate, MOMENT_RATE * momentRateSlice, SOURCE_FILE_INPUT * ret_src_in );
void addMomenteRate(  GRID grid,  SOURCE_FILE_INPUT src_in, WAVE hW, float * Jac, long long * srcIndex, MOMENT_RATE momentRate, MOMENT_RATE momentRateSlice, int it, int irk, float DT, float DH, float * gaussFactor, int nGauss, int flagSurf );

void finish_MultiSource( long long * srcIndex, MOMENT_RATE momentRate, MOMENT_RATE momentRateSlice,  long long npts );

void calculateMomentRate( SOURCE_FILE_INPUT src_in, MEDIUM medium, float * Jac, MOMENT_RATE momentRate, long long * srcIndex, float DH );

void allocContravariant( GRID grid, CONTRAVARIANT * contravariant );
void freeContravariant( CONTRAVARIANT contravariant );
void allocJac( GRID grid, float ** Jac );
void freeJac( float * Jac );
void allocMat3x3( GRID grid, Mat3x3 * _rDZ_DX, Mat3x3 * _rDZ_DY );
void freeMat3x3( Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY );
#ifdef FREE_SURFACE
void solveContravariantJac( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, SEND_RECV_DATA sr, CONTRAVARIANT contravariant, COORD coord, float * Jac, 
MEDIUM medium, Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY );
#else
void solveContravariantJac( MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, GRID grid, SEND_RECV_DATA sr, CONTRAVARIANT contravariant, COORD coord, float * Jac );
#endif



void allocWave( GRID grid, WAVE * h_W, WAVE * W, WAVE * t_W, WAVE * m_W );
void freeWave( WAVE h_W, WAVE W, WAVE t_W, WAVE m_W );

void waveDeriv( GRID grid, WAVE h_W, WAVE W, CONTRAVARIANT contravariant, MEDIUM medium, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3 );

void freeSurfaceDeriv( GRID grid, WAVE h_W, WAVE W, CONTRAVARIANT contravariant, 	
	MEDIUM medium, float * Jac, Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY, 
#ifdef PML
	PML_BETA pml_beta,
#endif
	int FB1, int FB2, int FB3 );

void pmlFreeSurfaceDeriv( GRID grid,
	WAVE h_W, WAVE W, CONTRAVARIANT contravariant, MEDIUM medium, AUX4 Aux4_1, AUX4 Aux4_2, 
	Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY, PML_D pml_d, MPI_BORDER border, int FB1, int FB2 );



void waveRk( GRID grid, int irk, float * h_W, float * W, float * t_W, float * m_W, float DT );

void allocPMLParameter( GRID grid, PML_ALPHA * pml_alpha, PML_BETA *pml_beta, PML_D * pml_d  );

void freePMLParamter( PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d );

void init_pml_parameter( PARAMS params, GRID grid, MPI_BORDER border, PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d );


void allocPML( GRID grid, AUX4 *Aux4_1, AUX4 *Aux4_2, MPI_BORDER border );



void freePML( MPI_BORDER border,  AUX4 Aux4_1, AUX4 Aux4_2 );




void pmlDeriv( 
	GRID grid,
	WAVE h_W, WAVE W, 
	CONTRAVARIANT contravariant,
	MEDIUM medium, AUX4 Aux4_1, AUX4 Aux4_2,	
	PML_ALPHA pml_alpha, PML_BETA pml_beta, PML_D pml_d,
	MPI_BORDER border, int FB1, int FB2, int FB3 );

void pmlRk( GRID grid, MPI_BORDER border, int irk, AUX4 Aux4_1, AUX4 Aux4_2, float DT );


void finalize_MPI( MPI_Comm * comm_cart );
void init_MPI( int *argc, char *** argv, PARAMS params, MPI_Comm * comm_cart, MPI_COORD * thisMPICoord, MPI_NEIGHBOR * mpiNeigbor );

void init_gpu( int PX, int PY, int PZ  );
void run( MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor, GRID grid, PARAMS params );
/*
void data2D_output_bin( PARAMS params, GRID grid, MPI_COORD thisMPICoord, float * data, const char * name );
void data2D_output_nc( PARAMS params, GRID grid, MPI_COORD thisMPICoord, COORD coord );
void dealSource( PARAMS params, GRID grid, COORD coord, MPI_COORD thisMPICoord );
void freeSourceIndex( SOURCE_INDEX srcIndex );
void readSource( PARAMS params, GRID grid, MPI_COORD thisMPICoord );
void dealWeisenShenMedium( PARAMS params, GRID grid, COORD coord, MPI_COORD thisMPICoord);
*/
void propagate( 
MPI_Comm comm_cart, MPI_COORD thisMPICoord, MPI_NEIGHBOR mpiNeighbor,
GRID grid, PARAMS params, SEND_RECV_DATA sr,
WAVE h_W, WAVE W, WAVE t_W, WAVE m_W,
Mat3x3 _rDZ_DX, Mat3x3 _rDZ_DY,
CONTRAVARIANT contravariant, float * Jac, MEDIUM medium, 
SOURCE_FILE_INPUT src_in, long long * srcIndex, MOMENT_RATE momentRate, MOMENT_RATE momentRateSlice,
SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu );

void loadPointSource( SOURCE S, WAVE h_W, int NX, int NY, int NZ, float * Jac, int it, int iRK, float DT, float DH, float rickerfc );

void GaussField( GRID grid, WAVE W );



void locateSlice( PARAMS params, GRID grid, SLICE * slice );

void locateSource( PARAMS params, GRID grid, SOURCE * source );
void locateFreeSurfSlice( GRID grid, SLICE * slice );

void allocSliceData( GRID grid, SLICE slice, SLICE_DATA * sliceData );
void freeSliceData( GRID grid, SLICE slice, SLICE_DATA sliceData );
void data2D_output_bin( GRID grid, SLICE slice, 
						MPI_COORD thisMPICoord, 
						float * datain, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu,	
						const char * name );
void data2D_XYZ_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, WAVE W, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu, char var, int it );
void data2D_Model_out( MPI_COORD thisMPICoord, PARAMS params, GRID grid, COORD coord, MEDIUM medium, SLICE slice, SLICE_DATA sliceData, SLICE_DATA sliceDataCpu );

#ifdef GPU_CUDA
void allocSliceData_cpu( GRID grid, SLICE slice, SLICE_DATA * sliceDataCpu );
void freeSliceData_cpu( GRID grid, SLICE slice, SLICE_DATA sliceDataCpu );
#endif






void allocSendRecv( GRID grid, MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA * sr );
void freeSendRecv( MPI_NEIGHBOR mpiNeighbor, SEND_RECV_DATA sr );

void mpiSendRecv( GRID grid, MPI_Comm comm_cart, MPI_NEIGHBOR mpiNeighbor, 
WAVE W, SEND_RECV_DATA sr );

void allocatePGV( GRID grid, PGV * pgv );
void freePGV( PGV pgv );
void allocatePGV_cpu( GRID grid, PGV * cpuPgv );
void freePGV_cpu( PGV cpuPgv );
void outputPGV( PARAMS params, GRID grid, MPI_COORDINATE thisMPICoord, PGV pgv, PGV cpuPgv );
void comparePGV( GRID grid, MPI_COORDINATE thisMPICoord, WAVE W, PGV pgv );

#endif //__FUNCTIONS__
