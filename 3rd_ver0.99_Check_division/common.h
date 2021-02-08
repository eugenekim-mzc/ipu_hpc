#pragma once

// #include <stdio.h>
// #include <vector>
// #include <sstream>

#define sim_dir			9
#define SIM_DIR_EXT		9
#define dim_s			2
#define cs				1.0/sqrt(3)
#define len				4

#define INPUT	(*pinput)
#define D2Q9	(*pmodel)

#define SAVE_PATH "Result/"

#define TECPLOTTYPE		0
#define CSVTYPE			1

#define tileNum 1470
#define div  	10000
#define div1 	10000

typedef struct _D2Q9model
{
	float			e[dim_s][sim_dir];
	int             eint[dim_s][sim_dir];
	int				oppi[sim_dir];
	float			tau;
	float			w[sim_dir];
} D2Q9model;

typedef struct _Inputval
{
	int h;
	int	nx, ny;
	int ngrid;
	int maxitr;
	int nwrite, irestart, irst;
	int m_itr;
	int icycle;
	float smagorinsky;
	float Re, length, nu;
	float pr, rho0, u0, v0, T0;
	int savetype;
} Inputval;

typedef struct _HostNode
{
	float *f;
	float *ftmp;
	float *rho;
	float *vel;
	float *prho;
	float *pvel;
} HostNode;

void initnode(HostNode *node, int length);
void freenode(HostNode *node);
