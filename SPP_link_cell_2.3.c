
/*
Author: Kuan
Date: 2016.09.03
Version: 2.3
Brief: moving boundary

os: windows 7/10
Compiler: MinGw
Compile command: gcc SPP_link_cell_2.3.c SPP_force.c SPP_sub.c -o jam -std=c99

os: Ubuntu
Compiler: gcc
Compile command: gcc SPP_link_cell_2.3.c SPP_force.c SPP_sub.c -o jam -std=gnu99 -lm




*/

#include "SPP.h"

real rzta=1./(3.*PI*ETA*D); //rzta=119.21718
char dir[80];
real kl_long,kl_short,k_fene;


int main(int argc, char *argv[])
{

	char ch; 
	while ((ch = getopt(argc, argv, "p:f")) != EOF)
		switch (ch) {
			case 'p': 
				printf("%c %s\n", ch, optarg);
				strcpy(dir,optarg);
				printf("%s\n", dir);
				break;
			case 'f':
				printf("%c\n", ch);
				break;
		}

    int nc[DIM];
    int N, pnc;

    real l[DIM], r_cut;
    real delta_t, t_end;
    inputParameters_LC(&delta_t, &t_end, &N, nc, l, &r_cut);
    pnc=1;
    for (int d=0; d<DIM; d++)
        pnc *= nc[d];
    Cell *grid = (Cell*)malloc(pnc*sizeof(Cell));

    Rod **order = malloc(N*sizeof(Rod*));

    initData_LC(N, grid, nc, l, order);
    timeIntegration_LC(0, delta_t, t_end, grid, nc, l, r_cut, N, order);
    freeLists_LC(grid, nc);

    free(grid);
    free(order);

    

    printf("Finish");
    return 0;

}

void inputParameters_LC(real *delta_t, real *t_end, int *N, int *nc, real *l, real *r_cut){

	char path[80];
	
	scanf("%d %lf %lf %lf",N,&kl_long,&kl_short,&k_fene);

    l[0] = 120.;
    l[1] = 120.;
    *delta_t = 5e-5; //5e-5
    *r_cut = 2.5+2*2;
    *t_end = 10.;

    nc[0] = (int)floor(l[0]/(*r_cut));
    nc[1] = (int)floor(l[1]/(*r_cut));
	
	int pnc=1;
    for (int d=0; d<DIM; d++)
        pnc *= nc[d];
	
	sprintf(path,"%s/para.txt",dir);
	
	FILE *file = fopen(path,"w");
	
	fprintf(file,"total N: %d\n",*N);
	fprintf(file,"t_end: %e\n",*t_end);
	fprintf(file,"delta_t: %e\n",*delta_t);
	fprintf(file,"kl_long: %e\n",kl_long);
	fprintf(file,"kl_short: %e\n",kl_short);
	fprintf(file,"k_fene: %e\n",k_fene);
	fprintf(file,"nc=[%d,%d] pnc=%d",nc[0],nc[1],pnc);
	
	fclose(file);

}

/* void RanInitial(Rod **order,int N)
{
	int num_chain = N/NR;
	int index_chain = 0;
	int index_bead = 0;
	Rod *temp_arr[NR+1];
	
	srand(time(NULL));
	
	for (index_chain=0;index_chain<=num_chain-1;index_chain++)
	{
		if(((rand()%2)==1) && (index_chain >= NS/NR))
		{			
			for (index_bead=1;index_bead<=NR;index_bead++)
			{
				temp_arr[index_bead] = order[NR*index_chain+index_bead-1];
				
			}
			for (index_bead=1;index_bead<=NR;index_bead++)
			{
				order[NR*index_chain+index_bead-1] = temp_arr[NR-index_bead+1];
				
			}
		}
	}
	
} */

void initData_LC(int N, Cell *grid, int *nc, real *l, Rod **order){

	char path[80];
    int pnc=1,kc[DIM],iniCol=21,yCount=1,ySign=1;
    real x0,y0,db=4*D+1.5;
    RodList *current;


    for (int d=0; d<DIM; d++)
        pnc *= nc[d];

    for (int i=0; i<pnc; i++)
        grid[i]=NULL;

    x0=0.5*db;
    y0=60.;

    for (int i=1; i<=N; i++){

        current = (RodList*)malloc(sizeof(RodList));
        order[i-1]=&(current->R);

        if (i==1){
            current->R.x[0]=x0;
            current->R.x[1]=y0;
            current->R.angle=0.;
			

        }
        else{
            current->R.x[0]=order[i-2]->x[0]+db;
            current->R.x[1]=0.;
            current->R.angle=0.;



            if (i%iniCol == 1) {
                current->R.x[0]=x0;
                current->R.x[1]=y0-(ySign*yCount*1.5);
                y0=current->R.x[1];
				ySign*=-1;
				yCount+=1;
            } else current->R.x[1]=y0;
        }
		
		current->R.hl=2.;
		current->R.x_1[0] = current->R.x[0]+current->R.hl*cos(current->R.angle);
		current->R.x_1[1] =	current->R.x[1]+current->R.hl*sin(current->R.angle);
		current->R.x_2[0] = current->R.x[0]-current->R.hl*cos(current->R.angle);
		current->R.x_2[1] = current->R.x[1]-current->R.hl*sin(current->R.angle);

        current->next=NULL;

        for (int d=0; d<DIM; d++)
            kc[d] = (int)floor(current->R.x[d] * nc[d] / l[d]);

        if(NULL==grid[index(kc,nc)])
            grid[index(kc,nc)]=current;
        else
            insertList(&grid[index(kc,nc)],current);

        
    }

    
	
	sprintf(path,"%s/ini.txt",dir);
	//RanInitial(order, N);
	
	FILE *file = fopen(path,"w");
	outputResults_LC(N,order,file);
	fclose(file);
	
    printf("initData_LC\n");

}

real makePBC(real x1, real x2, real xsize){
    real x=x1-x2;

    if (fabs(x) > 0.5*xsize) {
        if (x > 0.) x=x-xsize;
        else if (x < 0.) x=x+xsize;
    } return x;
}

real makePBC_pst(real x1, real x2, real xsize){
    real x=x2-x1;
	real x2_return;

    if (fabs(x) > 0.5*xsize) {
        if (x > 0.) x2_return=x2-xsize;
        else if (x < 0.) x2_return=x2+xsize;
    } return x2_return;
}



void compF_Intr(Rod **order, int N, real *l, long *seed, real thrust){

    for (int i=1;i<=N;i++){
     
        //Gaussian
        real ranForce_project;
		real rzta_project[2] = {0.188/KT,0.154/KT} ;
        real dc_project[3]={0.188,0.154,0.111};
        real invDT=1./DT;
		
		

        for (int d=0; d<DIM; d++){
            ranForce_project=sqrt(2.*dc_project[d]*invDT)*gasdev(seed);  //about 140
            order[i-1]->F[d] += ranForce_project;
        }
		
		order[i-1]->T += sqrt(2.*dc_project[3]*invDT)*gasdev(seed)*0.01;
		
        //if (i == N) printf("%e\n",ranForce);

        //Thrusting force
        
        /* real thrustAng;

        thrustAng=getAngle(i,order,l,seed);
        order[i-1]->F[0] += rzta*thrust*cos(thrustAng);
        order[i-1]->F[1] += rzta*thrust*sin(thrustAng);  */
		
		order[i-1]->F[0] += rzta_project[0]*thrust;

		
		


    }
}

void compF_LC(Cell *grid, int *nc, int *nc_moving, real *l, real r_cut, int PBC_state) {
    int ic[DIM], kc[DIM], kc_temp[DIM],flag;
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
            for (RodList *i=grid[index(ic,nc)]; NULL!=i; i=i->next) {
                for (int d=0; d<DIM; d++) //set force as 0
					i->R.F[d]=0.;
				i->R.T=0.;
				
                    
				
				//Reflecting Boundaries
				/* Rod wall_Rod;
				real r_wall=0.;
				int flag_wall=0;
				
				if(PBC_state)
					for (int d=0; d<DIM; d++){
						flag_wall=0;
						if(ic[d]==(nc_moving[d]-1)){
							wall_Rod.x[d] = 2*l[d]-(i->p.x[d]);
							r_wall = l[d]-(i->p.x[d]);
							flag_wall=1;
						}
						if(ic[d]==0){
							wall_Rod.x[d] = -(i->p.x[d]);
							r_wall = i->p.x[d];
							flag_wall=1;
						}
						if(flag_wall)	
							if(2*r_wall<=r_cut)
								force_WCA_wall(&i->p, &wall_Rod, sqr(2*r_wall), d, l);
					}	  */
				
				
				//Deal with the neighbors
                for (kc[0]=ic[0]-1; kc[0]<=ic[0]+1; kc[0]++)
                    for (kc[1]=ic[1]-1; kc[1]<=ic[1]+1; kc[1]++)
#if 3==DIM
                        for (kc[2]=ic[2]-1; kc[2]<=ic[2]+1; kc[2]++)
#endif
                        {
                            //treat kc[d]<0 and kc[d]>=nc[d] according to boundary conditions;

                            //PBC
                            flag=0;
							for (int d=0; d<DIM; d++){
					
								
								kc_temp[d]=kc[d];
								if (kc[d]<0){
                                    kc_temp[d]=nc_moving[d]-1;
                                    flag=1;
								}

								if (kc[d]>=nc_moving[d]){
                                    kc_temp[d]=0;
                                    flag=1;
								}
							}

							//if (distance of i->p to cell kc <= r_cut)
							for (RodList *j=grid[index(kc_temp,nc)];NULL!=j; j=j->next)
								if (i!=j) {
									real r = 0;
									real rod_dist_return[5];
									Rod ghost_rod;
									
									
									if (flag==0){																			
										for(int d=0; d<DIM; d++){											
											ghost_rod.x[d]=j->R.x[d];   
											ghost_rod.x_1[d]=j->R.x_1[d];
											ghost_rod.x_2[d]=j->R.x_2[d];
		
										}
																																												
									}else{ //PBC

										for(int d=0; d<DIM; d++){
											ghost_rod.x[d]=makePBC_pst(i->R.x[d],j->R.x[d], l[d]);   //here could be optimized
											ghost_rod.x_1[d]=makePBC_pst(i->R.x_1[d],j->R.x_1[d], l[d]);
											ghost_rod.x_2[d]=makePBC_pst(i->R.x_2[d],j->R.x_2[d], l[d]);
											
										}
																								
									} 
									rod_dist(&i->R,&ghost_rod,rod_dist_return);	
									r = sqr(rod_dist_return[0]);	
									
											
									if (r<=sqr(2.5)){
										real arm[DIM]={rod_dist_return[1]-i->R.x[0],rod_dist_return[2]-i->R.x[1]}; //arm of force
										
										/* printf("x_collision=(%f,%f)\n",rod_dist_return[1],rod_dist_return[2]);
										printf("x=(%f,%f)\n",i->R.x[0],i->R.x[1]);
										printf("arm=(%f,%f)\n",arm[0],arm[1]);		 */								
										
										if(r <= sqr(1.0e-6)){
											printf("sqr(r)=%f\n",r);
											
											printf("i->x=(%f,%f)\n",i->R.x[0],i->R.x[1]);
											printf("j->x=(%f,%f)\n",j->R.x[0],j->R.x[1]);
											
											printf("\n");	
											
											printf("i->x_1=(%f,%f)\n",i->R.x_1[0],i->R.x_1[1]);
											printf("i->x_2=(%f,%f)\n",i->R.x_2[0],i->R.x_2[1]);
											
											printf("\n");	
											
											printf("j->x_1=(%f,%f)\n",j->R.x_1[0],j->R.x_1[1]);
											printf("j->x_2=(%f,%f)\n",j->R.x_2[0],j->R.x_2[1]);
											
											printf("\n");	
											
											printf("i.collision=(%f,%f)\n",rod_dist_return[1],rod_dist_return[2]);
											printf("j.collision=(%f,%f)\n",rod_dist_return[3],rod_dist_return[4]);
											printf("dist=%f\n",sqr(rod_dist_return[3]-rod_dist_return[1])+sqr(rod_dist_return[4]-rod_dist_return[2]));
											
											printf("\n");	
											
											printf("r=%f\n",rod_dist_return[0]);
											
												
											exit(-1);
										}
										
										force_WCA(&i->R, &ghost_rod, rod_dist_return, arm, r, l);
									}
							}
                }
            }
}

/* void PBC_fail(Cell *grid, int *nc)
{
	int ic[DIM];
	char path[80];
	sprintf(path,"%s/PBC_fail.txt",dir);
	
	FILE *file = fopen(path,"w");
	
	for (ic[0]=0; ic[0]<nc[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc[2]; ic[2]++)
#endif
                for (RodList *i=grid[index(ic,nc)]; NULL!=i; i=i->next)
				{
                    fprintf(file,"%e    %e    %e    %e    %e    %e    %e    %e\n"\
					,i->p.x_old[0],i->p.x_old[1]\
					,i->p.F_wca[0],i->p.F_wca[1]\
					,i->p.F_bend[0],i->p.F_bend[1]\
					,i->p.F_fene[0],i->p.F_fene[1]);
				
				}
	
	fclose(file);
	
} */

void moveRods_LC(Cell *grid, int *nc, int *nc_moving, real *l) {
    int ic[DIM], kc[DIM];
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
        {   RodList **q = &grid[index(ic,nc)]; // pointer to predecessor
            RodList *i = *q;
            while (NULL != i) {

                //PBC
                for (int d=0; d<DIM; d++){
                    if (i->R.x[d] >= l[d]) i->R.x[d]-=l[d];
                    else if (i->R.x[d] < 0.) i->R.x[d]+=l[d];

                    if (i->R.x[d]<0. || i->R.x[d]>=l[d]){
                        printf("PBC fix fail\n");
                        printf("i->p.x=(%f,%f)\n",i->R.x[0],i->R.x[1]);
						printf("i->p.x_1=(%f,%f)\n",i->R.x_1[0],i->R.x_1[1]);
						printf("i->p.x_2=(%f,%f)\n",i->R.x_2[0],i->R.x_2[1]);
						printf("Force=(%f,%f)\n",i->R.F[0],i->R.F[1]);
						printf("Torque=%f\n",i->R.T);
						printf("l[0]=%e\n",l[0]);
						//PBC_fail(grid,nc);
                        exit(-2);
                    }
                }
				
				for (int d=0; d<DIM; d++){
                    kc[d] = (int)floor(i->R.x[d]*nc_moving[d]/l[d]);
                }


                
                if ((ic[0]!=kc[0])||(ic[1]!=kc[1])
                #if 3==DIM
                || (ic[2]!=kc[2])
                #endif
                ) {
                    deleteList(q);
                    insertList(&grid[index(kc,nc)], i);
                } else q = &i->next;
                i = *q;
            }
    }
}

void updateX(Rod *R, real delta_t) {

 
	real F_xy[DIM];
	
	vector_project(R->F, -(R->angle), F_xy);
		
	for (int d=0; d<DIM; d++) {		
		R->x[d] += DT * F_xy[d];	 
	}
	
	R->angle += DT * R->T;
	R->x_1[0] = R->x[0]+R->hl*cos(R->angle);
	R->x_1[1] =	R->x[1]+R->hl*sin(R->angle);
	R->x_2[0] = R->x[0]-R->hl*cos(R->angle);
	R->x_2[1] = R->x[1]-R->hl*sin(R->angle);
	
}

/* void updateV(Rod *R, real delta_t) {
	real a = delta_t * .5 / p->m;
	for (int d=0; d<DIM; d++)
		R->v[d] += R->F[d];
} */

void compX_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t) {
    int ic[DIM];
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
                for (RodList *i=grid[index(ic,nc)]; NULL!=i; i=i->next)
                    updateX(&i->R, delta_t);
    moveRods_LC(grid, nc, nc_moving, l);
}

/* void compV_LC(Cell *grid, int *nc, int *nc_moving, real *l, real delta_t) {
    int ic[DIM];
    for (ic[0]=0; ic[0]<nc_moving[0]; ic[0]++)
        for (ic[1]=0; ic[1]<nc_moving[1]; ic[1]++)
#if 3==DIM
            for (ic[2]=0; ic[2]<nc_moving[2]; ic[2]++)
#endif
                for (RodList *i=grid[index(ic,nc)]; NULL!=i; i=i->next)
                    updateV(&i->R, delta_t);
} */


void outputResults_LC(int N, Rod **order, FILE *file){

    for (int i=0;i<N;i++){
        fprintf(file,"%e    %e    %e    %e    %e    %e\n",order[i]->x[0],order[i]->x[1],order[i]->x_1[0],order[i]->x_1[1],order[i]->x_2[0],order[i]->x_2[1]);
    }

}

/*
void moving_boundary(Cell* grid, int *nc, int *nc_moving, real* l, real dx, real r_cut){
	
	int nc_old[2];
	for(int d=0;d<DIM;d++){
		nc_old[d]=nc_moving[d];
		l[d]-=dx;
		nc_moving[d] = (int)floor(l[d]/r_cut);	
	}
	
	int ic[DIM];
	int ic_target[DIM];
	if((nc_old[0]!=nc_moving[0])||(nc_old[1]!=nc_moving[1])){
		
				
		RodList *temp;
		ic[1]=nc_old[1]-1;
		ic_target[1]=nc_moving[1]-1;
		for(ic[0]=0;ic[0]<nc_old[0];ic[0]++){	
			ic_target[0]=ic[0];
			for (RodList *i=grid[index(ic,nc)]; NULL!=i;){
				temp=i->next;
				insertList(&grid[index(ic_target,nc)], i);
				i=temp;
			}
			grid[index(ic,nc)]=NULL;
		}
		
		ic[0]=nc_old[0]-1;
		ic_target[0]=nc_old[0]-2;
		for(ic[1]=0;ic[1]<nc_old[1]-1;ic[1]++){	
			ic_target[1]=ic[1];
			for (RodList *i=grid[index(ic,nc)]; NULL!=i;){
				temp=i->next;
				insertList(&grid[index(ic_target,nc)], i);
				i=temp;
			}
			grid[index(ic,nc)]=NULL;	
		}
			
		
			
				
	}
	
	

	
	
	
}*/

void timeIntegration_LC(real t, real delta_t, real t_end, Cell* grid, int *nc, real* l, real r_cut, int N, Rod **order) {

    int count=0,file_name=0;
    char path[12];
	int PBC_state=0; 
	real dx=0.0001;
	real thrust=0.8;
	int nc_moving[DIM];
	for(int d=0;d<DIM;d++)
		nc_moving[d]=nc[d];
	

    long storeTime,*seed=NULL;
    storeTime=-time(NULL);
    seed=&storeTime;

	//compF_LC(grid,nc,r_cut);
	while (t < t_end) {
		t += delta_t;
		
		
		//if(((float)count)/200==250.) PBC_state=0;
		/* if(l[0]<90. && PBC_state){
			PBC_state=0;
			thrust=0.05;
		} */
			
		
		compF_LC(grid,nc,nc_moving,l,r_cut,PBC_state);
		compF_Intr(order,N,l,seed,thrust);
		//compV_LC(grid,nc,nc_moving,l,delta_t);
		compX_LC(grid,nc,nc_moving,l,delta_t);
		
/* 		if(l[0]>90.){
			if(l[0]<100.)
				dx=0.00002;
			moving_boundary(grid,nc,nc_moving,l,dx,r_cut);
		} */

		count++;
		if(count%200 == 0){
            printf("%d\n",file_name);
            sprintf(path,"%s/%03d.txt",dir,file_name);
            FILE *file = fopen(path, "w");
            outputResults_LC(N, order, file);
            fclose(file);
            file_name++;
		}
	}
	

}


