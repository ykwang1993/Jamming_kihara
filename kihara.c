#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DIM 2
#define dot(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1])
#define cross(v1,v2) ((v1)[0]*(v2)[1]-(v1)[1]*(v2)[0])
#define length(v) sqrt((v)[0]*(v)[0]+(v)[1]*(v)[1])
#define swap(x, y) do{int tmp = x; x = y; y = tmp;} while(0)

typedef double real;

typedef struct {

	real x[DIM]; // position
	real x_1[DIM]; //terminal of rod
	real x_2[DIM];
	real angle[DIM];
	
	real v[DIM]; // velocity
	real F[DIM]; // force
	
} rod;

void point_rod_dist(real *p, rod *r, real *info){

	
	real v[DIM], v1[DIM], v2[DIM];
	for(int d=0;d<DIM;d++){
		v[d]=r->x_2[d] - r->x_1[d];
		v1[d]=p[d] - r->x_1[d];
		v2[d]=p[d] - r->x_2[d];
	}
	
	if(dot(v,v1) < 0.){
		
		info[0]=length(v1);
		info[1]=r->x_1[0];
		info[2]=r->x_1[1];
	
	}else if(dot(v,v2) > 0.){
		
		info[0]=length(v2);
		info[1]=r->x_2[0];
		info[2]=r->x_2[1];
	
	}else{
		
		info[0]=fabs(cross(v,v1))/length(v);
		
		//find the collision point
		real a=v[1];
		real b=-v[0];
		real c=-a*(r->x_1[0])-b*(r->x_1[1]); //c=-ax-by
		
		info[1]=p[0]-a*(a*p[0]+b*p[1]+c)/(a*a+b*b);
		info[2]=p[1]-b*(a*p[0]+b*p[1]+c)/(a*a+b*b); 
		
	}

}

void rod_dist(rod *r1, rod *r2){
	
	//intersect checking
	
	real dist[4][3]; //[min_dist,x,y]
	point_rod_dist(r1->x_1,r2,&dist[0][0]);
	point_rod_dist(r1->x_2,r2,&dist[1][0]);
	point_rod_dist(r2->x_1,r1,&dist[2][0]);
	point_rod_dist(r2->x_2,r1,&dist[3][0]);
	
	//for sorting the min dist 
	int idx[]={0,1,2,3};
	real min_dist[4];
	for(int i=0;i<4;i++)
		min_dist[i]=dist[i][0];
	
	
	for(int i = 0; i < 4; i++)
        for(int j = 1; j < 4 - i; j++)
            if(min_dist[j] < min_dist[j-1]){
				swap(idx[j], idx[j-1]);
				swap(min_dist[j], min_dist[j-1]);
			}
				

            

        

    
	
	for(int i=0;i<4;i++)
		printf("%f %f %f\n",dist[i][0],dist[i][1],dist[i][2]);
	for(int i=0;i<4;i++)
		printf("%d ",idx[i]);
}

int main(){
	
	rod r1,r2;

	r1.x_1[0] = 0.;
	r1.x_1[1] = 0.;
	r1.x_2[0] = 5.;
	r1.x_2[1] = 0.;
	
	r2.x_1[0] = 0.;
	r2.x_1[1] = 3.;
	r2.x_2[0] = 2.;
	r2.x_2[1] = 3.;
	
	rod_dist(&r1,&r2);
	
	
	//real x1[] = {2.,2.};
	
	//real info[3];
	//point_rod_dist(x1,&r1,info);
	
	//printf("min_dist:%f\n", info[0]);
	//printf("collision point:(%f,%f)\n",info[1],info[2]);
	
	return 0;
}