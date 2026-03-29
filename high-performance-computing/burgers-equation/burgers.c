/*==========================================================//
//  ANSI-C code (un optimised) for PX457 assignment 2 2025  //
//  Evolves a function u in 2D via finite differences       //
//  of the 2D Burgers Equation (simplified)                 //
//                                                          //
//  d u          du   du           d^2       d^2            //
//  ----- = -u ( -- + -- ) + nu ( ---- u  +  ---- u )       //
//  d t          dx   dy          dx^2       dy^2           //
//                                                          //
//  Based originally on code created by D. Quigley          //
//  Adapted by N. Hine and T. Latham                        //
//==========================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "makePNG.h"     /* For visualisation */
#include "mt19937ar.h"   /* Random number generator */

/* Function prototypes for memory management routines */
void allocate2d(double ***a,int num_rows,int num_cols);
void free2d(double ***a,int num_rows);

int main () {

  /* Function u on new and current grid */
  double **u_new, **u;

  /* Approximate Laplacian */
  double laplacian, grad;

  /* Number of grid points */
  int Nx = 256;
  int Ny = 256;

  /* Loop counters */
  int ix,iy,istep;

  /* Filename to which the grid is drawn */
  int  isnap=0;
  char filename[25];

  /*--------------*/
  /* Initial time */
  /*--------------*/
  clock_t t1 = clock();

  /*------------------------------------*/
  /* Initialise random number generator */
  /*------------------------------------*/
  unsigned long seed = 120549784972;
  init_genrand(seed);

  /*--------------------------*/
  /* Set grid spacing         */
  /*--------------------------*/
  double dx = 1.0;
  double dy = 1.0;

  /*-----------------------------*/
  /* Set timestep and K          */
  /*-----------------------------*/
  double dt = 0.001;
  double nu = 5.0;

  /*--------------------------*/
  /* Number of steps to run   */
  /*--------------------------*/
  int nstep = 10000;

  /*--------------------------------------*/
  /* Allocate memory for a bunch of stuff */
  /*--------------------------------------*/
  allocate2d(&u,Nx,Ny);
  allocate2d(&u_new,Nx,Ny);

  /*--------------------------------*/
  /* Initialise with random numbers */
  /*--------------------------------*/
  for(ix=0;ix<Nx;ix++) {
    for(iy=0;iy<Ny;iy++) {
      u[ix][iy] = 2.0*genrand() - 1.0;
    }
  }
      
  /*------------------------------------*/
  /* Write an image of the initial grid */
  /*------------------------------------*/
  int stepincr = 10; 
  sprintf(filename,"snapshot%08d.png",isnap); 
  writePNG(filename,u,Nx,Ny); 
  isnap++;

  /* setup time */
  clock_t t2 = clock();
  printf("Setup time                    : %15.6f seconds\n",(double)(t2-t1)/(double)CLOCKS_PER_SEC);
  t1 = t2;

  /*===============================*/
  /* BEGIN SECTION TO BE OPTIMISED */
  /*===============================*/
  // Pre-calculated constants to reduce the number of calculations
  const double inverse_dx = 1.0/dx;
  const double inverse_dy = 1.0/dy;
  const double half_inverse_dx = 0.5*inverse_dx;
  const double half_inverse_dy = 0.5*inverse_dy;
  const double inverse_dx_squared = inverse_dx*inverse_dx;
  const double inverse_dy_squared = inverse_dy*inverse_dy;


  /*------------------------------------------*/
  /* Loop over the number of output timesteps */
  /*------------------------------------------*/
  for (istep=1;istep<nstep;istep++) {
    
    /*-----------------------*/
    /* Loop over grid points */
    /*-----------------------*/
    for(iy=0;iy<Ny;iy++) {
	// use modulo arithmetic to find next y index
	int next_y = (iy+1)%Ny;
	int previous_y = (iy-1+Ny) % Ny;
	    
        for(ix=0;ix<Nx;ix++) {
	    // use modulo arithmetic to find next x index
	    int next_x = (ix+1)%Nx;
	    int previous_x = (ix-1+Nx)%Nx;
	    // get current and surrounding values
	    double u_current = u[ix][iy];
	    double u_right = u[next_x][iy];
	    double u_left = u[previous_x][iy];
	    double u_up = u[ix][next_y];
	    double u_down = u[ix][previous_y];

	    // compute derivatives
	    // grad term (du/dx+du/dy)
	    grad = (u_right - u_left)*half_inverse_dx + (u_up - u_down) * half_inverse_dy;
	    laplacian = (u_left + u_right - 2.0*u_current)*inverse_dx_squared + (u_down + u_up - 2.0*u_current)*inverse_dy_squared;
	    // update u_new values
	    u_new[ix][iy] = u_current - dt*u_current*grad + dt*nu*laplacian;
      }
    }
    // pointer swapping for O(1) allocation   
    double **temporary_grid = u;
    u = u_new;
    u_new = temporary_grid;

    /*-----------------------------*/
    /* Snapshots of grid to file   */
    /*-----------------------------*/
    if ( istep==isnap)  { 
        sprintf(filename,"snapshot%08d.png",isnap); 
        writePNG(filename,u,Nx,Ny); 
        isnap *= stepincr; 
    }  

  }

  /*=============================*/
  /* END SECTION TO BE OPTIMISED */
  /*=============================*/
  
  /* calculation time */
  t2 = clock();
  printf("Time taken for %8d steps : %15.6f seconds\n",nstep,(double)(t2-t1)/(double)CLOCKS_PER_SEC);
    
  /*----------------------------------*/
  /* Write an image of the final grid */
  /*----------------------------------*/
  sprintf(filename,"snapshot%08d.png",istep); 
  writePNG(filename,u,Nx,Ny); 

  /*--------------------------------------------*/
  /* Write final time-evolved solution to file. */
  /*--------------------------------------------*/
  FILE *fp = fopen("final_grid.dat","w");
  if (fp==NULL) printf("Error opening final_grid.dat for output\n");

  for(ix=0;ix<Nx-1;ix++) {
    for(iy=0;iy<Ny-1;iy++) {
      /* x and y at the current grid points */
      double x = dx*(double)ix;
      double y = dy*(double)iy;
      fprintf(fp,"%8.4f %8.4f %8.4e\n",x,y,u[ix][iy]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  /* Release memory */
  free2d(&u,Nx);
  free2d(&u_new,Nx);
  
  return 0;
  
}
  


/*===========================================*/
/* Auxilliary routines for memory management */ 
/*===========================================*/
void allocate2d(double ***a,int Nx,int Ny) {

  double **b_loc; 

  b_loc = (double **)calloc(Nx,sizeof(double *));
  if (b_loc==NULL) printf("malloc error in allocate2d\n"); 

  int iy;
  for (iy=0;iy<Nx;iy++) {
    
    b_loc[iy] = (double *)calloc(Ny,sizeof(double));
    if (b_loc[iy]==NULL) printf("malloc error for row %d of %d in allocate2d\n",iy,Nx);

  }

  *a = b_loc;

}

void free2d(double ***a,int Nx) {

  int iy;

  double **b_loc = *a;

  /* Release memory */
  for (iy=0;iy<Nx;iy++) { 
    free(b_loc[iy]);
  }
  free(b_loc);
  *a = NULL;

}








