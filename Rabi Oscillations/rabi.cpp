// g++ -o rabi rabi.cpp feynman.c -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm -lpthread -lconfig

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <pthread.h>
#include <math.h>
#include <libconfig.h>
#include "feynman.h"
#include "molecules.h"
#include "lasers.h"
#include "microtraps.h"

using namespace std;

#define NUM_THREADS 1

pthread_mutex_t readmutex, writemutex;

struct laser _6um;
struct microtrap m;
double omega0, mu, GammaDecay, BranchingRatio;
double El0, w01, lambda, lambda_cb, M_PIcb, w0_sq, w0_4th, lpiw, lpiwsq, omega, wy6um;


int main (int argc, char *argv[]){
        if (argc != 1) {
                fprintf(stderr, "Usage: just define the parameters within the config.cfg file and execute the rabi_ensemble.sh from the shell.\n");
                return 1;
        }

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

        // take parameters from config.cfg
        _6um.power=config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.power"));
        _6um.w0=config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.w0"));
        _6um.nu=config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.nu"));
        _6um.lambda=2.99792458e8 / config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.nu"));
        _6um.omega=2.0*M_PI*config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.nu"));
        _6um.x0=config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.x0"));
        _6um.y0=config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.y0"));
        _6um.z0=config_setting_get_float(config_lookup(&cfg, "lasers.QCL_6um.z0"));

        m.t=config_setting_get_float(config_lookup(&cfg, "microtraps.t"));
        m.y0=config_setting_get_float(config_lookup(&cfg, "microtraps.y0"));
        m.z0=config_setting_get_float(config_lookup(&cfg, "microtraps.z0"));
        m.vel=config_setting_get_float(config_lookup(&cfg, "microtraps.vel"));

        // w(y) QCL calculated at the y-coordinate of the microchip
        wy6um = _6um.w0*sqrt(1.+((((m.y0-_6um.y0)*_6um.lambda)/(M_PI*_6um.w0 *_6um.w0))*(((m.y0-_6um.y0)*_6um.lambda)/(M_PI*_6um.w0*_6um.w0))));

        omega0=2.0*M_PI*config_setting_get_float(config_lookup(&cfg, "mol_params.nu0"));
        mu=config_setting_get_float(config_lookup(&cfg, "mol_params.mu"));
        GammaDecay=1.0 / config_setting_get_float(config_lookup(&cfg, "mol_params.a3PItau"));             // decay rate expresses in MHz ( inverse of timelife expressed in us),
        BranchingRatio=config_setting_get_float(config_lookup(&cfg, "mol_params.a3PIBranchingRatio"));    // percentage of decays from the upper state towards the lower state

        El0 = 0.69258063 * sqrt( _6um.power/(_6um.w0*_6um.w0) );  // V/um  0.69258063 V/um = sqrt (mW 4.0/(c epsilon0 pi um^2))
        w01 = _6um.w0;                                            // um
        lambda = _6um.lambda;                                     //um
        lambda_cb=lambda*lambda*lambda;
        M_PIcb=M_PI*M_PI*M_PI;
        w0_sq = w01*w01;
        w0_4th = w0_sq*w0_sq;
        lpiw = lambda/(M_PI*w0_sq);                               //lambda0 = 5.9 um, 100 um waist
        lpiwsq = lpiw*lpiw;                                       //lambda0 = 5.9 um, 100 um waist
        omega = _6um.omega;

        // rabi dynamics calculation
        struct mol cc;
        char buf[0x1000];

        gsl_odeiv2_system sys = {func, jac, 6, &cc};
        //  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
        gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4imp, 1e-6, 1e-6, 0.0);

        int i, molnr=0, nsteps=40000;
        double t, t1, ti, tinit, tfin, tlaser, dtlaser;
        double temp[nsteps];
        double y[6];

        tlaser = (_6um.z0 - m.z0)/m.vel;
        dtlaser = (3.0*wy6um)/m.vel;
        tinit = m.t + tlaser - dtlaser;
        tfin = m.t + tlaser + dtlaser;

        for(i = 1; i <= nsteps; i++) temp[i]=0.0;

        for (;;) {
                pthread_mutex_lock (&readmutex);
                if (fgets(buf, sizeof(buf), stdin) == NULL) {
                        pthread_mutex_unlock (&readmutex);
                        break;
                }
                pthread_mutex_unlock (&readmutex);
                if (buf[0] == '#') {
                        cout<<"#time [us]"<<"\t"<<"Population v1"<<endl;
                        continue;
                }
                sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf",
                       &cc.molname, &cc.tinit, &cc.pos[0], &cc.pos[1], &cc.pos[2], &cc.vel[0], &cc.vel[1], &cc.vel[2]);

                //System described by the density matrix rho=(rho_uu,rho_ul;rho_lu,rho_ll) : u=upper level v1, l=lower level v0
                // y[0]=Real(upper state population)=upper state population=rho_uu
                // y[1]=Real(lower state population)=lower state population=rho_ll
                // y[2]=Real(rho_lu);
                // y[3]=Imaginary(rho_lu);
                // y[4]=Real(rho_ul);
                // y[5]=Imaginary(rho_ul);

                molnr++;
                y[0]=0.0;
                y[1]=1.0;         // initialized in the v0 state
                y[2]=0.0;
                y[3]=0.0;
                y[4]=0.0;
                y[5]=0.0;

                t = tinit;

                for (i = 1; i <= nsteps; i++) {
                        //every time it starts from the t of last iteration, i.e. last ti
                        ti = tinit + i * (tfin-tinit) / nsteps;
                        //    double ti = t + i * t1 / 1000.0;
                        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
                        if (status != GSL_SUCCESS) {
                                printf ("error, return value=%d\n", status);
                                break;
                        }
                        temp[i] += y[0];
                        //printf ("%.5e\t %.5e\t %.5e\n", t, y[0], y[1]);         // y[0] and y[1] are the population of the v1 and v0 state
                }

        }
        for (i = 1; i <= nsteps; i++) {
                ti = tinit + i * (tfin-tinit) / nsteps;
                printf ("%.5e %.5e\n", ti, temp[i]/molnr);
        }
        pthread_exit(NULL);
        gsl_odeiv2_driver_free (d);

        return 0;
}
