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
#define pdex(T,i,j,k) ((i == 0) ? ( (T[i+1][j][k] - T[i][j][k])/dx ) : ( (i == nx-1) ? ( (T[i-1][j][k] - T[i][j][k])/dx ) : 0 ))
#define pdey(T,i,j,k) ((j == 0) ? ( (T[i][j+1][k] - T[i][j][k])/dy ) : ( (j == ny-1) ? ( (T[i][j-1][k] - T[i][j][k])/dy ) : 0 ))
#define pdez(T,i,j,k) ((k == 0) ? ( (T[i][j][k+1] - T[i][j][k])/dz ) : ( (k == nz-1) ? ( (T[i][j][k-1] - T[i][j][k])/dz ) : 0 ))

#define isin(idx, nidx) ( idx !=0 && idx != nidx-1)

#define needs_boundary(i,j,k,i_left,i_right) (i == i_left | j== 0 | k == 0 | i == i_right-1  | j == ny-1 | k == nz-1)

// #define temp(a,i,j,k) *(a + i*ny*nz + j*nz + k)

enum {nx = (int) (1/dx + 1),
      ny = (int) (1/dy + 1),
      nz = (int) (1/dz + 1)};






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
    // int local_nx = scounts[rank]/(ny*nz);

    int i_left = displs[rank]/(ny*nz), i_right = (displs[rank] + scounts[rank])/(ny*nz);
    // if(rank==0)
    //     printf("nx,ny,nz, %d,%d,%d", nx,ny,nz);
    // printf("displs, scounts, rank, i_left, i_right = %d,%d,%d,%d,%d\n", displs[rank],scounts[rank],rank,i_left, i_right);


    // Defining T_loc
    double T[nx][ny][nz]={0};
    double T_new[nx][ny][nz] = {0};

    // Defining initial temperatures
    for(i=i_left-1;i<i_right+1;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {
                if(i<0)
                    continue;
                if(i>=nx)
                    continue;
                T[i][j][k] = 300;
                
            }
    
    memcpy(T_new[i_left], T[i_left], (i_right-i_left)*ny*nz*sizeof(double));
    
    // Note: T[i_left-1] = T_down, T[i_right+1] = T_up
    

    double t=0.0; double t_save=0.0;
    double t_max = strtod(argv[1], NULL);
    for(t=0; t < t_max; t += dt)
    {

        if(fabs(t-t_save)<dt_save/2)
        {
            int l = (int) (t_save/dt_save);

            for(j=0;j<ny;j++)
                for(k=0;k<nz;k++)
                {
                    i=i_left-1;
                    if(i>=0)
                        T_new[i][j][k] = 0.0;
                    i=i_right+1;
                    if(i<nx)
                        T_new[i][j][k] = 0.0;
                }
            if(rank==0)
                for(i=i_right+1;i<nx;i++)
                    for(j=0;j<ny;j++)
                        for(k=0;k<nz;k++)
                            T_new[i][j][k]=0.0;
            
            if(rank==0)
                MPI_Reduce(MPI_IN_PLACE, T_new, nx*ny*nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            else
                MPI_Reduce(T_new, T_new, nx*ny*nz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            

            if(rank==0)
            {
                FILE *fpt;

                char filename[15];
                sprintf(filename, "mpi_test.csv.%d", l);

                fpt = fopen(filename, "w+");

                x=0; y=0; z=0;

                for(i=0;i<nx;x=i*dx,i++)
                    for(j=0;j<ny;y=j*dy,j++)
                        for(k=0;k<nz;z=k*dz,k++)
                        {
                            fprintf(fpt, "%+.2f, %+.2f, %+.2f, %+.4f \n", x, y, z, T_new[i][j][k]);
                            // fprintf(fpt, "T[%d][%d][%d]: \t %+.4lf\n",i,j,k,T[i][j][k]);
                        }
                
                fclose(fpt);
                printf("%f\n", t_save);

            }
            t_save += dt_save;
        }

        if(rank!=0)
        {
            MPI_Recv(&T[i_left-1], ny*nz, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
            MPI_Send(&T[i_left], ny*nz, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
        }
        if(rank!=size-1)
        {
            MPI_Send(&T[i_right], ny*nz, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
            MPI_Recv(&T[i_right+1], ny*nz, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
        }


        // For interior points:
        for(i=i_left;i<i_right;i++)
            for(j=0;j<ny;j++)
                for(k=0;k<nz;k++)
                {
                    // Check if we are on the boundary
                    if(needs_boundary(i,j,k,i_left,i_right) == 0){
                        T_new[i][j][k] = T[i][j][k] +
                            kappa*dt*(pd2x(T,i,j,k,1) +  pd2y(T,i,j,k,1) + pd2z(T,i,j,k,1))/(rho * cp);
                    }
                    else{
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

        memcpy(T[i_left], T_new[i_left], (i_right-i_left)*ny*nz*sizeof(double));

    }

    // for(i=0; i<nx; i++)
    //     for(j=0; j<ny; j++)
    //         for(k=0; k<nz; k++)
    //             printf("T[%d][%d][%d]: \t %+.4lf\n",i,j,k,T[i][j][k]);

    free(displs);
    free(scounts);

    MPI_Finalize();

    return 0;
}
