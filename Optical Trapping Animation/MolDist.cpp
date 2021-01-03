// compile with: g++ -o MolDist.out MolDist.cpp Function.c -O2 -lm -lgsl -lgslcblas -Wall
// The output is the set of coordinates of N molecules: (x,y,z,vx,vy,vz)
// Space coordinates: you can choose molecules uniformly distributed within a cylinder or within an ellipsoid
// Speed coordinates: you can choose dvx, dvy, dvz (the FWHM per each Gaussian distribution)

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <sys/time.h>
#include <string.h>

#include "Parameters.hpp"                                                         // All defined parameters
#include "Functions.hpp"                                                          // All defined functions

using namespace std;

int main (){
        int i, microseconds, Excitation=1;                                      // All molecules are initialized at Excitation=1 (a3pi1 state)
        double y[6], phi, theta, rho, sigmav, Microtrap_X;
        srand(time(NULL));

        Microtrap_X=-3.*IPG_WY(Chip_DY);                                        // Mimimum chip-trap distance [µm]
        sigmav=SIGMA_V(Microtrap_T);                                            // σ [µm/µs] for the Maxwell-Boltzmann distribution of speed: σv=σvx=σvy=σvz

        struct timeval st;
        gettimeofday(&st,NULL);
        microseconds=st.tv_usec;

        //set the random number generators
        gsl_rng * r1=gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r1,microseconds);
        gsl_rng * r2=gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r2,microseconds+42);
        gsl_rng * r3=gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r3,microseconds+97);

        cout<<"#name"<<"\t\t"<<"t"<<"\t\t\t\t"<<"x"<<"\t\t\t"<<"y"<<"\t\t\t"<<"z"<<"\t\t\t"<<"vx"<<"\t\t\t"<<"vy"<<"\t\t\t"<<"vz"<<"\t\t"<<"Exc"<<"\t\t"<<"color"<<endl;
        if (strcmp( Microtrap_Shape, "Cyl" )==0)
        {
                for( i=0; i< Number_Molecules; i++ )
                {
                        //Uniform cylindrical distribution
                        phi=2*M_PI*gsl_rng_uniform(r1);
                        rho=sqrt(gsl_rng_uniform(r2));                          //the sqrt function is to obtain the uniform distribution inside a cyrcle of unitary radius
                        y[3] = Microtrap_X + Microtrap_DX*rho*sin(phi);         //x  a*rho*sin(phi) is to scale to the x dimension of the ellipse
                        y[4] = Microtrap_Y + 2.*Microtrap_DY*(gsl_rng_uniform(r3)-0.5); //y
                        y[5] = Microtrap_Z + Microtrap_DZ*rho*cos(phi);         //z  c*rho*cos(phi) is to scale to the z dimension of the ellipse
                        //Gaussian distribution on vx, vy, vz
                        y[0] = Microtrap_V+gsl_ran_gaussian(r1, sigmav);        //vx
                        y[1] = gsl_ran_gaussian(r2, sigmav);                    //vy
                        y[2] = gsl_ran_gaussian(r3, sigmav);                    //vz

                        printf("%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %d\t %d\n", i, 0., y[3], y[4], y[5], y[0], y[1], y[2], Excitation, RGB_QS(Excitation) );
                }
        }
        else if (strcmp( Microtrap_Shape, "Ell" )==0)
        {
                for( i=0; i< Number_Molecules; i++ )
                {
                        //Uniform ellipsoidal distribution
                        theta=2*M_PI*gsl_rng_uniform(r1);
                        phi=acos(2*gsl_rng_uniform(r2)-1);
                        rho=cbrt(gsl_rng_uniform(r3));                          //the cbrt function is to obtain the uniform distribution inside a sphere of unitary radius
                        y[3] = Microtrap_X + Microtrap_DX*rho*cos(theta)*sin(phi); //x     a*rho*cos(theta)*sin(phi) is to scale to the x dimension of the ellipse
                        y[4] = Microtrap_Y + Microtrap_DY*rho*sin(theta)*sin(phi); //y     b*rho*sin(theta)*sin(phi) is to scale to the y dimension of the ellipse
                        y[5] = Microtrap_Z + Microtrap_DZ*rho*cos(phi);         //z     c*rho*cos(phi) is to scale to the z dimension of the ellipse
                        //Gaussian distribution on vx, vy, vz
                        y[0] = Microtrap_V+gsl_ran_gaussian(r1, sigmav);        //vx
                        y[1] = gsl_ran_gaussian(r2, sigmav);                    //vy
                        y[2] = gsl_ran_gaussian(r3, sigmav);                    //vz

                        printf("%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %d\t %d\n", i, 0., y[3], y[4], y[5], y[0], y[1], y[2], Excitation, RGB_QS(Excitation) );
                }
        }

        return 0;
}
