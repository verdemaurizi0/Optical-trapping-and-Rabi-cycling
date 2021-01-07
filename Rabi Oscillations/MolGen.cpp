// compile with: g++ -o MolGen MolGen.cpp-L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lconfig

// The output is the set of coordinates of N molecules: (x,y,z,vx,vy,vz)
// Space coordinates: you can choose molecules uniformly distributed within a cylinder or within an ellipsoid
// Speed coordinates: you can choose dvx, dvy, dvz (the FWHM per each Gaussian distribution)

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <libconfig.h>
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
#include "microtraps.h"

using namespace std;

int main (){
        int i, microseconds;
        double y[6], phi, theta, rho, sigmav;
        struct microtrap m;
        srand(time(NULL));
        struct timeval st;
        gettimeofday(&st,NULL);
        microseconds=st.tv_usec;

        config_t cfg;

        config_init(&cfg);
        // if config.cfg is not present an error msg occurs
        if(!config_read_file(&cfg, "config.cfg"))
        {
                fprintf(stderr, "%s:%d - %s\n", config_error_file(&cfg),
                        config_error_line(&cfg), config_error_text(&cfg));
                config_destroy(&cfg);
                return(EXIT_FAILURE);
        }

        m.num=config_setting_get_float(config_lookup(&cfg, "microtraps.num"));
        m.t=config_setting_get_float(config_lookup(&cfg, "microtraps.t"));
        m.x0=config_setting_get_float(config_lookup(&cfg, "microtraps.x0"));
        m.y0=config_setting_get_float(config_lookup(&cfg, "microtraps.y0"));
        m.z0=config_setting_get_float(config_lookup(&cfg, "microtraps.z0"));
        m.dx=config_setting_get_float(config_lookup(&cfg, "microtraps.dx"));
        m.dy=config_setting_get_float(config_lookup(&cfg, "microtraps.dy"));
        m.dz=config_setting_get_float(config_lookup(&cfg, "microtraps.dz"));
        m.vel=config_setting_get_float(config_lookup(&cfg, "microtraps.vel"));
        m.temp=config_setting_get_float(config_lookup(&cfg, "microtraps.temp"));

        //set the random number generators
        gsl_rng * r1=gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r1,microseconds);
        gsl_rng * r2=gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r2,microseconds+42);
        gsl_rng * r3=gsl_rng_alloc (gsl_rng_mt19937);
        gsl_rng_set(r3,microseconds+97);

        sigmav = sqrt(0.2969453*m.temp); // σ [µm/µs] for the Maxwell-Boltzmann distribution of speed: σv=σvx=σvy=σvz

        cout<<"#name"<<"\t\t"<<"t"<<"\t\t\t\t"<<"x"<<"\t\t\t"<<"y"<<"\t\t\t"<<"z"<<"\t\t\t"<<"vx"<<"\t\t\t"<<"vy"<<"\t\t\t"<<"vz"<<endl;
        for( i=0; i< m.num; i++ )
        {
                //Uniform cylindrical distribution
                phi=2*M_PI*gsl_rng_uniform(r1);
                rho=sqrt(gsl_rng_uniform(r2));                              //the sqrt function is to obtain the uniform distribution inside a cyrcle of unitary radius
                y[3] = m.x0 + m.dx*rho*sin(phi);                            //x  a*rho*sin(phi) is to scale to the x dimension of the ellipse
                y[4] = m.y0 + 2.*m.dy*(gsl_rng_uniform(r3)-0.5);            //y
                y[5] = m.z0 + m.dz*rho*cos(phi);                            //z  c*rho*cos(phi) is to scale to the z dimension of the ellipse
                //Gaussian distribution on vx, vy, vz
                y[0] = gsl_ran_gaussian(r1, sigmav);                //vx
                y[1] = gsl_ran_gaussian(r2, sigmav);                        //vy
                y[2] = m.vel + gsl_ran_gaussian(r3, sigmav);                        //vz

                printf("%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\n", i, m.t, y[3], y[4], y[5], y[0], y[1], y[2]);
        }

        return 0;
}
