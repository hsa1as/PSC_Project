int DEBUG = 0;
/*
Serial
Program to study unsteady heat conduction in 3 dimensions.

Consider a box of size 1m * 1m * 1m
The box extends from 0,0,0 to 1,1,1

The box is at an initial temperature of 300K
The box is placed in a medium of 800K

Heat conduction to be studied in the box over a period of time.
Numerical calculation done using GS Scheme


*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#ifdef _OPENMP
#include<omp.h>
#endif



// Defining Constants
#define PI 3.14159265359 
#define rho 2700
#define kappa 237
#define cp 0.9
#define h 242.3
#define Ta 800
#define dx 0.04
#define dx_sq dx*dx
#define dy 0.04
#define dy_sq dy*dy
#define dz 0.04
#define dz_sq dz*dz
#define dt 0.001
#define dt_save 10     // dt_save - time step for saving output
#define c 1.0
// Useful helper macros
#define pd2x(T,i,j,k,in) ((in==0) ? 0 : (T[i+1][j][k] - 2*T[i][j][k] + T[i-1][j][k])/dx_sq)
#define pd2y(T,i,j,k,in) ((in==0) ? 0 : (T[i][j+1][k] - 2*T[i][j][k] + T[i][j-1][k])/dy_sq)
#define pd2z(T,i,j,k,in) ((in==0) ? 0 : (T[i][j][k+1] - 2*T[i][j][k] + T[i][j][k-1])/dz_sq)
#define pdex(T,i,j,k) ((i == 0) ? ( (T[i+1][j][k] - T[i][j][k])/dx ) : ( (i == nx-1) ? ( (T[i-1][j][k] - T[i][j][k])/dx ) : 0 ))
#define pdey(T,i,j,k) ((j == 0) ? ( (T[i][j+1][k] - T[i][j][k])/dy ) : ( (j == ny-1) ? ( (T[i][j-1][k] - T[i][j][k])/dy ) : 0 ))
#define pdez(T,i,j,k) ((k == 0) ? ( (T[i][j][k+1] - T[i][j][k])/dz ) : ( (k == nz-1) ? ( (T[i][j][k-1] - T[i][j][k])/dz ) : 0 ))

#define isin(idx, nidx) ( idx !=0 && idx != nidx-1)

#define needs_boundary(i,j,k) (i == 0 | j== 0 | k == 0 | i == nx-1 | j == ny-1 | k == nz-1)

enum {nx = (int) (1/dx + 1),
      ny = (int) (1/dy + 1),
      nz = (int) (1/dz + 1)};

double T[nx][ny][nz] = {0};
double T_new[nx][ny][nz] = {0};

int main(int argc, char* argv[])
{
    if(argc < 2){
        printf("Usage: %s stop_time\n", argv[0]);
        exit(1);
    }
    for(int i = 0; i < argc; i++){
        if(!strcmp(argv[i],"-d")) DEBUG = 1;
    }
    int i, j, k;
    double x, y, z;
    // Defining initial temperatures
    for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {
                T[i][j][k] = 300;
            }

    double t=0.0; double t_save = 0.0;
    double t_max = strtod(argv[1], NULL);
    double start = omp_get_wtime();
    for(t=0; t < t_max; t += dt)
    {

        /*if(fabs(t-t_save)<dt/2)
        {
            int l = (int) (t_save/dt_save);

            FILE *fpt;

            char filename[15];
            sprintf(filename, "plot.csv.%d", l);

            fpt = fopen(filename, "w+");

            x=0; y=0; z=0;

            for(i=0;i<nx;x=i*dx,i++)
                for(j=0;j<ny;y=j*dy,j++)
                    for(k=0;k<nz;z=k*dz,k++)
                    {
                        fprintf(fpt, "%+.2f, %+.2f, %+.2f, %+.4f \n", x, y, z, T[i][j][k]);
                        // fprintf(fpt, "T[%d][%d][%d]: \t %+.4lf\n",i,j,k,T[i][j][k]);
                    }
            
            fclose(fpt);
            printf("%f\n", t_save);

            t_save += dt_save;
        }*/
        // For interior points:
        for(i=0;i<nx;i++)
            for(j=0;j<ny;j++)
                for(k=0;k<nz;k++)
                {
                    // Check if we are on the boundary
                    if(needs_boundary(i,j,k) == 0){
                        T_new[i][j][k] = T[i][j][k] +
                            kappa*dt*(pd2x(T,i,j,k,1) +  pd2y(T,i,j,k,1) + pd2z(T,i,j,k,1))/(rho * cp);
                    }else{
                        // If we are in the boundary, the heat balance equation changes. We can write a generalised equation for the heat flux
                        // using selectors that would use the index to determine whether each term in the sequence is active or not. This would
                        // simplify coding edge cases greatly
                        T_new[i][j][k] = kappa*dx*dy*dz*( pd2x(T,i,j,k,isin(i,nx)) + pd2y(T,i,j,k,isin(j,ny)) + pd2z(T,i,j,k,isin(k,nz)) )
                            + h*(Ta - T[i][j][k])*( dy*dz*(!isin(i,nx)) + dx*dz*(!isin(j,ny)) + dx*dy*(!isin(k,nz)) )
                            + kappa*(dy*dz*pdex(T,i,j,k) + dx*dz*pdey(T,i,j,k) + dx*dy*pdez(T,i,j,k));
                        T_new[i][j][k] /= rho*cp*dx*dy*dz/dt;
                        T_new[i][j][k] += T[i][j][k];
                        if(DEBUG){
                            printf("DEBUG: T[%d][%d][%d] @ t = %lf: \t %+.3lf\n", i ,j, k, t, T[i][j][k]);
                            printf("DEBUG: pdex(T,%d,%d,%d) @ t = %lf: \t %+.3lf\n", i, j, k, t, pdex(T,i,j,k));
                            printf("DEBUG: pdey(T,%d,%d,%d) @ t = %lf: \t %+.3lf\n", i, j, k, t, pdey(T,i,j,k));
                            printf("DEBUG: pdez(T,%d,%d,%d) @ t = %lf: \t %+.3lf\n", i, j, k, t, pdez(T,i,j,k));

                        }
                    }
               }

        memcpy(T, T_new, nx*ny*nz*sizeof(double));


    }
    double end = omp_get_wtime();
    for(i=0; i<nx; i++)
        for(j=0; j<ny; j++)
            for(k=0; k<nz; k++)
                printf("T[%d][%d][%d]: \t %+.4lf\n",i,j,k,T[i][j][k]);
    printf("[*] %s took %lf seconds to solve for t_max = %lf", argv[0], end-start, t_max);
    

    return 0;
}
