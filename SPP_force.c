
#include "SPP.h"

/* void force_WCA_wall(Rod *i, Rod *j, real r, int d, real *l){
	
	real sigma=1.;
    real epsilon=4.12e-3;

    real s = sqr(sigma) / r;
    s = sqr(s) * s;
    real f = 24 * epsilon * s / r * (1 - 2 * s);
	
	i->F[d] += f * rzta * makePBC(j->x[d],i->x[d],l[d]);
	
	
} */


void force_WCA(Rod *i, Rod *j, real *rod_dist_return, real *arm, real r, real *l) {
    real sigma=1.;
    real epsilon=4.12e-3;

    real s = sqr(sigma) / r;
    s = sqr(s) * s;
    real f = 24 * epsilon * s / r * (1 - 2 * s);
	real F_wca[DIM];
	
    for (int d=0; d<DIM; d++)		
		F_wca[d] = f * makePBC(rod_dist_return[3+d],rod_dist_return[1+d],l[d]); //rod_dist_return[3] and [4] is the pst of x2

	
	real F_project[DIM];
	real rzta_project[2] = {0.188/KT,0.154/KT} ;
	vector_project(F_wca, i->angle, F_project);
	
	for (int d=0; d<DIM; d++)
		i->F[d] += rzta_project[d]*F_project[d];
	
	i->T += cross(arm,F_wca) * 0.111/KT;
	
 	if(sqrt(r)<0.5){
		printf("crash\n");
		printf("wca=(%f,%f)\n",F_wca[0],F_wca[1]);
		printf("Torque=%f\n",cross(arm,F_wca));
		printf("r=%f\n",sqrt(r));
		
		printf("i->x_1=(%f,%f)\n",i->x_1[0],i->x_1[1]);
		printf("i->x_2=(%f,%f)\n",i->x_2[0],i->x_2[1]);
		
		printf("j->x_1=(%f,%f)\n",j->x_1[0],j->x_1[1]);
		printf("j->x_2=(%f,%f)\n",j->x_2[0],j->x_2[1]);
	} 
	
}




real vecAngle(real x, real y)
/* Return the angle of a vector with respect to the horizontal axis in
   a value within the close interval [-PI,PI]. The input values are the
   x and y component of the vector, respectively. */
{
    real ang;

    if (fabs(x) <= 1.0e-6) { /* avoid NaN when calculating y/x with x=0 */
        if (y > 0.) ang=0.5*PI;
        else if (y < 0.) ang=-0.5*PI;
    } else if (x > 1.0e-6) ang=atan(y/x);
    else if (x < -1.0e-6) {
        if (y >= 0.) ang=atan(y/x)+PI;
        else if (y < 0.) ang=atan(y/x)-PI;
    } return ang;
}

void point_rod_dist(real *p, Rod *r, real *info){

	
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

/*
rod_dist_return [min_dist, r1_x, r1_y, r2_x, r2_y]
*/
void rod_dist(Rod *r1, Rod *r2, real* rod_dist_return){
	
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
	
	rod_dist_return[0]=min_dist[idx[0]];
	
	switch(idx[0]){
		case 0:
			rod_dist_return[1]=r1->x_1[0];
			rod_dist_return[2]=r1->x_1[1];
			rod_dist_return[3]=dist[idx[0]][1];
			rod_dist_return[4]=dist[idx[0]][2];
			break;
		case 1:
			rod_dist_return[1]=r1->x_2[0];
			rod_dist_return[2]=r1->x_2[1];
			rod_dist_return[3]=dist[idx[0]][1];
			rod_dist_return[4]=dist[idx[0]][2];
			break;
		case 2:
			rod_dist_return[1]=dist[idx[0]][1];
			rod_dist_return[2]=dist[idx[0]][2];
			rod_dist_return[3]=r2->x_1[0];
			rod_dist_return[4]=r2->x_1[1];
			break;
		case 3:
			rod_dist_return[1]=dist[idx[0]][1];
			rod_dist_return[2]=dist[idx[0]][2];
			rod_dist_return[3]=r2->x_2[0];
			rod_dist_return[4]=r2->x_2[1];
			break;
	}	
}


/*
argument angle is the angle between rod and positive x axis.
if the input vector is (1,1), and angle is pi/4,
then the output vector is (1.414, 0)

that is, output = (project on long axis, project on short axis)
*/
void vector_project(real *input, real angle, real *output){
	
	output[0]=input[0]*cos(angle)+input[1]*sin(angle);
	output[1]=-input[0]*sin(angle)+input[1]*cos(angle);
	
}





