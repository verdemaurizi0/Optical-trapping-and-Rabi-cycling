// compile with: g++ -o Trajectories.out Trajectories.cpp Functions.cpp -O2 -lm -lgsl -lgslcblas -Wall

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <string.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

#include "Parameters.hpp"                                                 // All defined parameters
#include "Functions.hpp"                                                  // All defined functions


// Quantum state dependent mehanical potential:
// - Excitation = 0 describes the ground state X1Sigma+
// - Excitation = 1 describes the excited state a3pi1
// - Excitation = 2 describes lost molecules (no force)

int Excitation = -1;

int func (double t, const double y[], double f[], void *params)
    {
    if (Excitation==0)                                                  // X1Sigma+: DC Stark U_S + polarizability potential U_P (in unit of COmass)
        { 
        double DU_P = gs_p * EL_IPG(y);
        double DU_S = -2. * k_gs * EL_STARK(y);
                                                              
        f[0] = DU_P * DX_EL_IPG(y) - DU_S * D_EL_STARK(y);              //a_x   dU_P/dEl_IPG * dEl_IPG/dx  +  dU_S/dEl_S * dEl_S/dx
        f[1] = DU_P * DY_EL_IPG(y);                                     //a_y   dU_P/dEl_IPG * dEl_IPG/dy
        f[2] = DU_P * DZ_EL_IPG(y);                                     //a_z   dU_P/dEl_IPG * dEl_IPG/dz
  
        f[3] = y[0]; // v_x
        f[4] = y[1]; // v_y
        f[5] = y[2]; // v_z
        return GSL_SUCCESS;
        }
    else if (Excitation==1)                                             // a3pi1: DC Stark U_S + polarizability potential U_P (in unit of COmass)
        {
        double DU_P = a3pi1_p * EL_IPG(y);                                                               
        double DU_S = a3pi1_Mu * a3pi1_Mu * EL_STARK(y) / hypot( Half_L , ( a3pi1_Mu * EL_STARK(y) ));

        f[0] = DU_P * DX_EL_IPG(y) - DU_S * D_EL_STARK(y);              //a_x   dU_P/dEl_IPG * dEl_IPG/dx  +  dU_S/dEl_S * dEl_S/dx
        f[1] = DU_P * DY_EL_IPG(y);                                     //a_y   dU_P/dEl_IPG * dEl_IPG/dy
        f[2] = DU_P * DZ_EL_IPG(y);                                     //a_z   dU_P/dEl_IPG * dEl_IPG/dz
  
        f[3] = y[0]; // v_x
        f[4] = y[1]; // v_y
        f[5] = y[2]; // v_z
        return GSL_SUCCESS;
        }
    else if (Excitation==2)                                             // other quantum states: lost molecules, no forces
        {                               
        f[0] = 0.; //a_x
        f[1] = 0.; //a_y
        f[2] = 0.; //a_z
  
        f[3] = y[0]; // v_x
        f[4] = y[1]; // v_y
        f[5] = y[2]; // v_z
        return GSL_SUCCESS;
        }
        else return 0; 
    }

int main ()
    {
    double t, tfin, h, MaxTimeStep, y[6];                               // y = {v_x, v_y, v_z, x, y, z}  f = {a_x, a_y, a_z, v_x, v_y, v_z}
    int microseconds, name, ExcitationAus;
    int a3pi1_N=Number_Molecules, Lost_N=0, GS_N=0, Trapped_N=0;
    char buf[0x200];                                                    // To contain 512 byte

    tfin=10.*(3.*IPG_WY(Chip_DY))/Microtrap_V;                           // Final time calculated as approximatively ten times the period needed to reach the optical trap 

    struct timeval st;
    gettimeofday(&st,NULL);
    microseconds=st.tv_usec;                                            // to seed the random number generators by taking the machine time expressed in µs

    //set the random number generators
    gsl_rng * r1=gsl_rng_alloc (gsl_rng_mt19937);                       // to account for the spontaneous decay of a3pi1 state 
    gsl_rng_set(r1, microseconds );
    gsl_rng * r2=gsl_rng_alloc (gsl_rng_mt19937);                       // to account for the RDL-driven transition towards the ground state
    gsl_rng_set(r2, microseconds + 44 );

    //from SDE
    const gsl_odeiv_step_type * T = gsl_odeiv_step_rkck;
    gsl_odeiv_step * s = gsl_odeiv_step_alloc (T,6);                    // 6 is the vector's size defined by gsl
    gsl_odeiv_control * c = gsl_odeiv_control_y_new (1.e-8,0.0);        // control eps err 10^-8
    gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (6);                  // 6 is the vector's size defined by gsl

    do 
        {
        sscanf(buf, "%d %lf %lf %lf %lf %lf %lf %lf %d", &name, &t, &y[3], &y[4], &y[5], &y[0], &y[1], &y[2], &ExcitationAus); 
        //Excitation is a external variable not defined in the main, it is then better to use an ausiliary variable 
        Excitation=ExcitationAus; 

        gsl_odeiv_system sys = {func,NULL,6};
                
        h = 1.e-6;
        
        while ( t < tfin ) 
            {
            if ( NEAR_RDL(y,Excitation) )     h=fminf(h,1.e-2);             // to control the time step when a3pi1 molecules transit to x1Sigma+ state
            
            if ( Excitation ==1   &&  y[0] < 0 )         break;             // stop the evolution for a3pi1 molecules which inverted their longitudinal motion                       

            if ( Excitation == 1  &&  gsl_rng_uniform(r1) >= exp(-h/a3pi1_T) ) 
                {
                Excitation=2;
                a3pi1_N--;
                Lost_N++;
                break;
                }

            if ( Excitation == 1  &&  INSIDE_RDL(y) )
                {
                if ( gsl_rng_uniform(r2) <= Trans_P ) 
                    {
                    Excitation=0;
                    a3pi1_N--;
                    GS_N++;
                    }
                else 
                    {
                    Excitation=2;
                    a3pi1_N--;
                    Lost_N++;
                    break;
                    }
                }

            gsl_odeiv_system sys = {func,NULL,6};
            int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &t, tfin, &h, y);
            if(status != GSL_SUCCESS) break;
            }
        
        if ( Excitation == 0 && GS_TRAPPED(y) ) Trapped_N++;

        } 
    while (fgets(buf, sizeof(buf), stdin) != NULL);

    printf("%lf\t %lf\n", IPG_W, (float)Trapped_N/(a3pi1_N + Lost_N + GS_N) );

gsl_odeiv_evolve_free (e);
gsl_odeiv_control_free (c);
gsl_odeiv_step_free (s);
  
return 0;
}

