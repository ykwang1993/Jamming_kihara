
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define DIM 2
#define sqr(x) ((x)*(x))
#define D 1.00 //bead diameter
#define ETA 0.89e-3
#define PI 3.141592653590
#define NS 4 //Length of long chain
#define NR 4 //Length of short chain
#define KT 4.12e-3
#define DT 5.0e-5
#define CP 0.05


#if 1==DIM
#define index(ic,nc) ((ic)[0])
#elif 2==DIM
#define index(ic,nc) ((ic)[0] + (nc)[0]*(ic)[1])
#elif 3==DIM
#define index(ic,nc) ((ic)[0] + (nc)[0]*((ic)[1] + (nc)[1]*(ic)[2]))
#endif



typedef double real;


/*
define for kihara
*/ 
#define dot(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1])
#define cross(v1,v2) ((v1)[0]*(v2)[1]-(v1)[1]*(v2)[0])
#define length(v) sqrt((v)[0]*(v)[0]+(v)[1]*(v)[1])
#define swap(x, y) do{real tmp = x; x = y; y = tmp;} while(0)
	

typedef struct {

	real x[DIM]; // position
	real x_1[DIM]; //terminal of rod
	real x_2[DIM];
	real angle;
	real hl; //half length of rod 
	
	real v[DIM]; // velocity
	real w; //angular velocity
	real F[DIM]; // force
	real T; //torque
	
} Rod;

typedef struct RodList {
    Rod R;
    struct RodList *next;
} RodList;

typedef RodList* Cell;

void point_rod_dist(real *p, Rod *r, real *info);
void rod_dist(Rod *r1, Rod *r2, real* rod_dist_return);



extern real rzta; //rzta=119.21718
extern char dir[80];
extern real kl_long,kl_short,k_fene;


void outputResults_LC(int N, Rod **order, FILE *file);
real ran1(long *idum);
real gasdev(long *idum);

void insertList(RodList **root_list, RodList *i);
void deleteList(RodList **q);
void freeLists_LC(Cell* grid, int *nc);

void inputParameters_LC(real *delta_t, real *t_end, int *N, int *nc, real *l, real *r_cut);
void RanInitial(Rod **order,int N);
void initData_LC(int N, Cell *grid, int *nc, real *l, Rod **order);
void PBC_fail(Cell *grid, int *nc);

real makePBC(real x1, real x2, real xsize);
real makePBC_pst(real x1, real x2, real xsize);
void moveRods_LC(Cell *grid, int *nc, int *nc_moving, real *l);

//void force_WCA_wall(Rod *i, Rod *j, real r, int d, real *l);
void force_WCA(Rod *i, Rod *j, real *rod_dist_return, real *arm, real r, real *l);



void compF_Intr(Rod **order, int N, real *l, long *seed, real thrust);
void compF_LC(Cell *grid, int *nc, int *nc_moving, real *l, real r_cut, int PBC_state);
void updateX(Rod *R, real delta_t);
//void updateV(Rod *p, real delta_t);
void compX_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t);
//void compV_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t);
//void moving_boundary(Cell* grid, int *nc, int *nc_moving, real* l, real dx, real r_cut);
void timeIntegration_LC(real t, real delta_t, real t_end, Cell* grid, int *nc, real* l, real r_cut, int N, Rod **order);












































