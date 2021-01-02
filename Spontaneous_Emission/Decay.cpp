// compile with: g++ -o Decay.out Decay.cpp -O2 -lm -lgsl -lgslcblas -Wall

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <string.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

#include "Parameters.hpp"
#include "Functions.hpp"

int func (double t, const double y[], double f[], void *params)
{
        // put here the force field if present
        f[0] = 0.; //a_x
        f[1] = 0.; //a_y
        f[2] = 0.; //a_z

        f[3] = y[0]; // v_x
        f[4] = y[1]; // v_y
        f[5] = y[2]; // v_z
        return GSL_SUCCESS;
}

int main (int argc, char *argv[])
{
        double t, tfin, h, MaxTimeStep, y[6];
        int a3pi1_N=0, countMol=0;
        int microseconds, name, Excitation;
        char buf[0x200];                                                // To contain 512 byte

        struct timeval st;
        gettimeofday(&st,NULL);
        microseconds=st.tv_usec;                                        // to seed the random number generators by taking the machine time expressed in µs

        //set the random number generators
        gsl_rng * r1=gsl_rng_alloc (gsl_rng_mt19937);                   // to account for the spontaneously decay of a3pi1 state
        gsl_rng_set(r1, microseconds );

        if(argc != 2)
        {
                fprintf(stderr, "Usage: <Final_Time_us> \nIt evolves molecules until a final time <Final_Time_us>.\nAll units in um and us.\nParameters:\n\t<Final_Time_us>: final time.\n");
                return 1;
        }

        tfin = atof(argv[1]);                                           // Final time input

        //from SDE
        const gsl_odeiv_step_type * T = gsl_odeiv_step_rkck;
        gsl_odeiv_step * s = gsl_odeiv_step_alloc (T,6);                // 6 is the vector's size defined by gsl
        gsl_odeiv_control * c = gsl_odeiv_control_y_new (1.e-8,0.0);    // control eps err 10^-8
        gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (6);              // 6 is the vector's size defined by gsl

        // to read and copy the comment in the new file
        if(fgets(buf, sizeof(buf),stdin)==NULL) {perror("Error Input Data\n"); return -1;}
        while (buf[0] == '#')
        {
                fputs(buf,stdout);
                if(fgets(buf, sizeof(buf),stdin)==NULL) {perror("Error Input Data\n"); return -1;}
        }

        do
        {
                sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf %d", &name, &t, &y[3], &y[4], &y[5], &y[0], &y[1], &y[2], &Excitation);
                //Excitation is a external variable not defined in the main, it is then better to use an ausiliary variable

                gsl_odeiv_system sys = {func,NULL,6};

                countMol++;
                if ( Excitation==1 ) a3pi1_N++;

                h = 1.e-6;
                MaxTimeStep=0.1;

                while ( t < tfin )
                {
                        h=fminf(h,MaxTimeStep);

                        if ( gsl_rng_uniform(r1) >= exp(-h/a3pi1_T) && Excitation == 1 )
                        {
                                a3pi1_N--;
                                Excitation=2;
                        }

                        gsl_odeiv_system sys = {func,NULL,6};
                        int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tfin, &h, y);
                        if(status != GSL_SUCCESS) break;
                }

                printf("%d\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %d\t %lf\n", name, t, y[3], y[4], y[5], y[0], y[1], y[2], Excitation, (float)a3pi1_N/countMol);
        }

        while (fgets(buf, sizeof(buf), stdin) != NULL);

        gsl_odeiv_evolve_free (e);
        gsl_odeiv_control_free (c);
        gsl_odeiv_step_free (s);

        return 0;
}
