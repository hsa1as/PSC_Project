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


// Defining Constants
#define PI 3.14159265359 
#define rho 8000
#define kappa 50
#define cp 420
#define h 242.3
#define Ta 800
#define dx 0.01
#define dx_sq dx*dx
#define dy 0.01
#define dy_sq dy*dy
#define dz 0.01
#define dz_sq dz*dz
#define dt 0.01
#define dt_save 0.5     // dt_save - time step for saving output
#define t_max 2.0       // t_max - time upto which computation is done
#define c 1.0
enum {nx = (int) (1/dx + 1),
      ny = (int) (1/dx + 1),
      nz = (int) (1/dx + 1),
      nt= (int) (t_max/dt + 1), 
      nt_save= (int) (t_max/dt_save + 1)};



int main(int argc, char* argv[])
{
    int i, j, k;
    double x, y, z;
    double T[nx][ny][nz] = {0};
    double T_new[nx][ny][nz] = {0};
    double T_save[nt_save][nx][ny][nz] = {0};


    // Defining initial temperatures 
    for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
            for(k=0;k<nz;k++)
            {
                T[i][j][k] = 300;
                // if(i>1 && i<nx-1)
                //     if(j>1 && j<ny-1)
                //         if(k>1 && k<nz-1)
                //             T[i][j][k] = 300;
                // else
                //     T[i][j][k] = 800;
            }
    

    
    // double var = 0.0;
    // int iter=0;
    // double norm_err = 0.0;
    double t=0.0, t_save=0.0;

    for(t=0; t<t_max; t+=dt)
    {

        if(fabs(t-t_save)<dt/2)
        {
            int l = (int) (t_save/dt_save);

            memcpy(T_save[k], T, nx*ny*nz*sizeof(double));

            t_save += dt_save;
        }

        for(i=1;i<nx-1;i++)
            for(j=1;j<ny-1;j++)
                for(k=1;k<nz-1;k++)
                {
                    if(i>1 && i<nx-1)
                        if(j>1 && j<ny-1)
                            if(k>1 && k<nz-1)
                                T_new[i][j][k] = T[i][j][k] + kappa*dt*( T[i+1][j][k] + T[i-1][j][k] + T[i][j+1][k] + T[i][j-1][k]
                                                    T[i][j][k+1] + T[i][j][k-1] - 6*T[i][j][k])/(dx_sq * rho * cp);
                    else
                        // T[i][j][k] = 800;
                        // T_new[i][j][k] = T[i][j][k] + kappa*dt*( T[i+1][j][k] + T[i-1][j][k] + T[i][j+1][k] + T[i][j-1][k]
                        //         T[i][j][k+1] + T[i][j][k-1] - 6*T[i][j][k])/(dx_sq * rho * cp);
                        T_new[i][j][k] = T[i][j][k] - (h*dx_sq*(T[i][j][k]-Ta)*dt)/(rho*cp);    
                }
        
        memcpy(T, T_new, nx*ny*nz*sizeof(double));


    }
    

    

    return 0;
}