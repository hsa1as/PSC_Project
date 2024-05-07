#define DEBUG 0
/*
Parallel MPI
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
#include<mpi.h>


// Defining Constants
#define PI 3.14159265359 
#define rho 2700
#define kappa 237
#define cp 0.9
#define h 242.3
#define Ta 800
#define dx 0.1
#define dx_sq dx*dx
#define dy 0.1
#define dy_sq dy*dy
#define dz 0.1
#define dz_sq dz*dz
#define dt 0.001
#define dt_save 0.2     // dt_save - time step for saving output
#define c 1.0
// Useful helper macros

#define pd2x(T,i,j,k,in) ((in==0) ? 0 : (T[i+1][j][k] - 2*T[i][j][k] + T[i-1][j][k])/dx_sq)
#define pd2y(T,i,j,k,in) ((in==0) ? 0 : (T[i][j+1][k] - 2*T[i][j][k] + T[i][j-1][k])/dy_sq)
#define pd2z(T,i,j,k,in) ((in==0) ? 0 : (T[i][j][k+1] - 2*T[i][j][k] + T[i][j][k-1])/dz_sq)
#define pdex(T,i,j,k,local_nx,rank,size) ((i==0&&rank==0) ? ( (T[i+1][j][k] - T[i][j][k])/dx ) : ( (i==nx-1&&rank==size-1) ? ( (T[i-1][j][k] - T[i][j][k])/dx ) : 0 ))
#define pdey(T,i,j,k) ((j == 0) ? ( (T[i][j+1][k] - T[i][j][k])/dy ) : ( (j == ny-1) ? ( (T[i][j-1][k] - T[i][j][k])/dy ) : 0 ))
#define pdez(T,i,j,k) ((k == 0) ? ( (T[i][j][k+1] - T[i][j][k])/dz ) : ( (k == nz-1) ? ( (T[i][j][k-1] - T[i][j][k])/dz ) : 0 ))

#define isin(idx, nidx) ( idx !=0 && idx != nidx-1)
#define isxin(idx, nidx, rank, size) ( (isin(idx,nidx)==0)? !(idx==0&&rank==0 || (idx==nidx-1)&&(rank==size-1)) : 1)

#define needs_boundary(i,j,k,local_nx) (i == 0 | j== 0 | k == 0 | i == local_nx-1  | j == ny-1 | k == nz-1)

// #define temp(a,i,j,k) *(a + i*ny*nz + j*nz + k)

enum {nx = (int) (1/dx + 1),
      ny = (int) (1/dy + 1),
      nz = (int) (1/dz + 1)};


double pd2x_spec(double ***T, int i,int j,int k,int local_nx,double T_u[ny][nz],double T_d[ny][nz],int rank,int size)
{
    if(i==0)
        if(rank==0)
            return 0.0;
        // return (T[i+1][j][k] - 2*T[i][j][k] + T_d[j][k])/dx_sq;
        return 0.0;
    
    if(i==local_nx-1)
        if(rank==size-1)
            return 0.0;
        return (T_u[j][k] - 2*T[i][j][k] + T[i-1][j][k])/dx_sq;
    
    return (T[i+1][j][k] - 2*T[i][j][k] + T[i-1][j][k])/dx_sq;
}






int main(int argc, char* argv[])
{
    int i, j, k;
    double x, y, z;

    //Initialize MPI
    MPI_Status status;
    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0)
    {
        if(argc < 2){
        printf("Usage: %s stop_time not given\n", argv[0]);
        exit(1);
        }
    }

    double t1 = MPI_Wtime();   // Starting time

    // Define and Initialize variables per proc
    int *displs = (int *)malloc(size*sizeof(int));
    int *scounts = (int *)malloc(size*sizeof(int));

    int rem = nx%size, qu = nx/size;
    int offset = 0;
    for(i=0; i<size; i++)
    {
        scounts[i] = qu*ny*nz;
        if(rem>0)
        {
            scounts[i] += ny*nz;
            rem--;
        }
        displs[i] = offset;
        offset += scounts[i];
    }

    int local_nx = scounts[rank]/(ny*nz);

    // Defining T_loc
    // double *T = (double *)malloc(local_nx*ny*nz*sizeof(double));
    // double *T_new = (double *)malloc(local_nx*ny*nz*sizeof(double));
    

    double*** T = (double ***)malloc(local_nx*sizeof(double **));
    double*** T_new = (double ***)malloc(local_nx*sizeof(double **));

    for(i=0;i<local_nx;i++)
    {
        T[i] = (double **)malloc(ny*sizeof(double *));
        T_new[i] = (double **)malloc(ny*sizeof(double *));
        
        for(j=0;j<ny;j++)
        {
            T[i][j] = (double *)malloc(nz*sizeof(double));
            T_new[i][j] = (double *)malloc(nz*sizeof(double));
        }
    }

    double T_up[ny][nz] = {0};
    double T_down[ny][nz] = {0};

    // double T[nx][ny][nz] = {0};
    // double T_new[nx][ny][nz] = {0};

    // Defining initial temperatures
    for(i=0;i<local_nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {
                if(i==0)
                {
                    T_up[j][k] = 300;
                    T_down[j][k] = 300;
                }
                T[i][j][k] = 300;
            }
    

    double t=0.0; double t_save=0.0;
    double t_max = strtod(argv[1], NULL);
    for(t=0; t < t_max; t += dt)
    {

        /*if(fabs(t-t_save)<dt/2)
        {
            int l = (int) (t_save/dt_save);

            // memcpy(T_save[k], T, nx*ny*nz*sizeof(double));
            FILE *fpt;

            char filename[15];
            sprintf(filename, "test.csv.%d", l);

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

        if(rank!=0)
        {
            MPI_Recv(&T_down, ny*nz, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&T[0], ny*nz, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if(rank!=size-1)
        {
            MPI_Send(&T[local_nx-1], ny*nz, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            MPI_Recv(&T_up, ny*nz, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
        }


        // For interior points:
        for(i=0;i<local_nx;i++)
            for(j=0;j<ny;j++)
                for(k=0;k<nz;k++)
                {
                    // Check if we are on the boundary
                    if(needs_boundary(i,j,k,local_nx) == 0){
                        T_new[i][j][k] = T[i][j][k] +
                            kappa*dt*(pd2x(T,i,j,k,1) +  pd2y(T,i,j,k,1) + pd2z(T,i,j,k,1))/(rho * cp);
                    }
                    else{
                        // If we are in the boundary, the heat balance equation changes. We can write a generalised equation for the heat flux
                        // using selectors that would use the index to determine whether each term in the sequence is active or not. This would
                        // simplify coding edge cases greatly
                        T_new[i][j][k] = kappa*dx*dy*dz*( pd2x_spec(T,i,j,k,local_nx,T_up,T_down,rank,size) + pd2y(T,i,j,k,isin(j,ny)) + pd2z(T,i,j,k,isin(k,nz)) )
                            + h*(Ta - T[i][j][k])*( dy*dz*(!isxin(i,nx,rank,size)) + dx*dz*(!isin(j,ny)) + dx*dy*(!isin(k,nz)) )
                            + kappa*(dy*dz*pdex(T,i,j,k,local_nx,rank,size) + dx*dz*pdey(T,i,j,k) + dx*dy*pdez(T,i,j,k));
                        T_new[i][j][k] /= rho*cp*dx*dy*dz/dt;
                        T_new[i][j][k] += T[i][j][k];
                        /*if(DEBUG){
                            printf("DEBUG: T[%d][%d][%d] @ t = %lf: \t %+.3lf\n", i ,j, k, t, T[i][j][k]);
                            printf("DEBUG: pdex(T,%d,%d,%d) @ t = %lf: \t %+.3lf\n", i, j, k, t, pdex(T,i,j,k,local_nx,rank,size));
                            printf("DEBUG: pdey(T,%d,%d,%d) @ t = %lf: \t %+.3lf\n", i, j, k, t, pdey(T,i,j,k));
                            printf("DEBUG: pdez(T,%d,%d,%d) @ t = %lf: \t %+.3lf\n", i, j, k, t, pdez(T,i,j,k));

                        }*/
                    }
               }
        
        

        // memcpy(T, T_new, nx*ny*nz*sizeof(double));
        for(i=0;i<local_nx;i++)
            for(j=0;j<ny;j++)
                for(k=0;k<nz;k++)
                    T[i][j][k] = T_new[i][j][k];


    }

    // for(i=0; i<nx; i++)
    //     for(j=0; j<ny; j++)
    //         for(k=0; k<nz; k++)
    //             printf("T[%d][%d][%d]: \t %+.4lf\n",i,j,k,T[i][j][k]);

    free(displs);
    free(scounts);
    free(T);
    free(T_new);

    MPI_Finalize();

    return 0;
}
