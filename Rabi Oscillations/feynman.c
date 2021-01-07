#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>
#include "feynman.h"
#include "molecules.h"
#include "lasers.h"



// Two-levels system interacting with a laser field
// omega0 is the transition frequency
// omega is the QCL frequency


double Efield, dtEfield, Delta, dtDelta;  // omega0 is related to the transition v0->v1, omega and omega2 are the QCL frequencies

void *fields (double *t, struct mol *cc){
        //the following variables depend on the molecules position and need to be recalculated at every step.
        double time = *t - cc->tinit;
        double x=cc->pos[0]+cc->vel[0]*time - _6um.x0; //di quale laser stiamo parlando qui? abbiamo duplicato le coordinate, riferendole a due diversi sistemi di riferimento (i centri dei due laser..)
        double y=cc->pos[1]+cc->vel[1]*time - _6um.y0;


        if(fabs(y)<1e-10) y=1e-10;
        double z=cc->pos[2]+cc->vel[2]*time - _6um.z0;
        double ylpiwsq = (1.0+y*y*lpiwsq);
        double expo =exp( (-x*x-z*z)/( w0_sq*ylpiwsq ));
        double frac = ( 1.0+M_PI*M_PI*w0_4th/(y*y*lambda*lambda) );
        double frac_sq = frac*frac;

        Efield = mu* El0*expo / hypot(1.0,y*lpiw);

        dtEfield =-El0 * cc->vel[1] * y * lambda*lambda * expo / ( M_PI*M_PI*w0_4th * sqrt(ylpiwsq)*ylpiwsq );
        dtEfield += El0 * expo * ( ( (-2.0 * cc->vel[0] * x - 2.0 * cc->vel[2] * z)/(w0_sq * ylpiwsq) ) - ( (2.0 * cc->vel[1] *y *(-x*x- z*z)*lpiwsq)/(w0_sq*ylpiwsq*ylpiwsq) ) ) / sqrt(ylpiwsq);
        dtEfield *= mu;


        Delta = omega - omega0 - (2.0 *M_PIcb* cc->vel[1] *(w0_4th)*( x*x + z*z ) ) / (
                (y*y*y*y)*(frac_sq)*lambda_cb) - ( 2.0 * M_PI* cc->vel[1] )/lambda + (M_PI * cc->vel[1] * (x*x + z*z) ) / (
                (y*y)*frac*lambda) - (M_PI*(2.0 * cc->vel[0] * x + 2.0 * cc->vel[2] * z))/( y*frac*lambda) + (cc->vel[1] * lambda)/(M_PI*w0_sq*ylpiwsq);

        dtDelta=  -(8.0 *(M_PIcb*M_PI*M_PI)* cc->vel[1] * cc->vel[1] * (w0_4th*w0_4th)*( x*x + z*z ) ) / (
                (y*y*y*y*y*y*y) * (frac_sq*frac)*(lambda_cb*lambda*lambda)) + (  10.0 *M_PIcb* cc->vel[1] * cc->vel[1] *(w0_sq*w0_sq)*( x*x + z*z ) ) / (
                ( y*y*y*y*y )*( frac_sq )*lambda_cb) - (  4.*M_PIcb * cc->vel[1] * (w0_4th)*(2.0 * cc->vel[0] * x + 2.0 * cc->vel[2] * z) ) / (
                ( y*y*y*y )*(frac_sq)*lambda_cb) - (2.0 *M_PI * cc->vel[1] * cc->vel[1] * ( x*x + z*z ) ) / (
                (y*y*y)*frac*lambda) + (2.0 * M_PI * cc->vel[1] * (2.0 * cc->vel[0] * x + 2.0 * cc->vel[2] * z) ) / (
                (y*y)*frac*lambda) - (M_PI*(2.0 * cc->vel[0] * cc->vel[0] + 2.0 * cc->vel[2] * cc->vel[2]) ) / (
                y*frac*lambda) - (  2.0 * cc->vel[1] * cc->vel[1] * y *lambda_cb)/(M_PIcb*(w0_4th*w0_sq)*(ylpiwsq*ylpiwsq) );

        return NULL;
}

int jac (double t, const double y[], double *dfdy,double dfdt[], void *params){
        struct mol cc = *(struct mol *)params;
        fields(&t, &cc);
        //temp:
        (void)(t); /* avoid unused parameter warning */
        //  double mu = *(double *)params;
        double sindelta=sin(Delta*t);
        double cosdelta=cos(Delta*t);
        double GammaDecay=0.0, BranchingRatio=0.0;


        gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 6, 6);
        gsl_matrix * m = &dfdy_mat.matrix;

        gsl_matrix_set (m, 0, 0, -GammaDecay );
        gsl_matrix_set (m, 0, 1, 0.0 );
        gsl_matrix_set (m, 0, 2, 0.0 );
        gsl_matrix_set (m, 0, 3, 0.0 );
        gsl_matrix_set (m, 0, 4, 0.0 );
        gsl_matrix_set (m, 0, 5, -Efield );

        gsl_matrix_set (m, 1, 0, BranchingRatio * GammaDecay );
        gsl_matrix_set (m, 1, 1, 0.0 );
        gsl_matrix_set (m, 1, 2, 0.0 );
        gsl_matrix_set (m, 1, 3, 0.0 );
        gsl_matrix_set (m, 1, 4, 0.0 );
        gsl_matrix_set (m, 1, 5, Efield );

        gsl_matrix_set (m, 2, 0, 0.0 );
        gsl_matrix_set (m, 2, 1, 0.0 );
        gsl_matrix_set (m, 2, 2, -0.5 * GammaDecay );
        gsl_matrix_set (m, 2, 3, Delta );
        gsl_matrix_set (m, 2, 4, 0.0 );
        gsl_matrix_set (m, 2, 5, 0.0 );

        gsl_matrix_set (m, 3, 0, -0.5 * Efield );
        gsl_matrix_set (m, 3, 1,  0.5 * Efield );
        gsl_matrix_set (m, 3, 2, -Delta );
        gsl_matrix_set (m, 3, 3, -0.5 * GammaDecay );
        gsl_matrix_set (m, 3, 4, 0.0 );
        gsl_matrix_set (m, 3, 5, 0.0 );

        gsl_matrix_set (m, 4, 0, 0.0 );
        gsl_matrix_set (m, 4, 1, 0.0 );
        gsl_matrix_set (m, 4, 2, 0.0 );
        gsl_matrix_set (m, 4, 3, 0.0 );
        gsl_matrix_set (m, 4, 4, -0.5 * GammaDecay );
        gsl_matrix_set (m, 4, 5, -Delta );

        gsl_matrix_set (m, 5, 0,  0.5 * Efield );
        gsl_matrix_set (m, 5, 1, -0.5 * Efield );
        gsl_matrix_set (m, 5, 2, 0.0 );
        gsl_matrix_set (m, 5, 3, 0.0 );
        gsl_matrix_set (m, 5, 4, Delta );
        gsl_matrix_set (m, 5, 5, -0.5 * GammaDecay );

        /* Debugging the matrix
           int j,yy;
           for (j = 0; j < 4; j++){
           for(yy=0; yy<4; yy++) printf ("%f\t",gsl_matrix_get (m,j,yy));
           printf("\n");
           }
         */

         dfdt[0] = -dtEfield * y[5];
         dfdt[1] = dtEfield * y[5];
         dfdt[2] = dtDelta * y[3];
         dfdt[3] = -dtDelta * y[2] - 0.5 * dtEfield * (y[0]-y[1]);
         dfdt[4] = -dtDelta * y[5];
         dfdt[5] = dtDelta * y[4] + 0.5 * dtEfield * (y[0]-y[1]);


        return GSL_SUCCESS;
}

int func (double t, const double y[], double f[], void *params){
        (void)(t); /* avoid unused parameter warning */
        struct mol cc = *(struct mol *)params;
        fields(&t, &cc);
        double sindelta=sin(Delta*t);
        double cosdelta=cos(Delta*t);
        double GammaDecay=0.0, BranchingRatio=0.0;

        f[0] = -Efield * y[5] - GammaDecay * y[0];
        f[1] = Efield * y[5] + BranchingRatio * GammaDecay * y[0];
        f[2] = Delta * y[3] - 0.5 * GammaDecay * y[2];
        f[3] = -Delta * y[2] - 0.5 * Efield * (y[0]-y[1]) - 0.5 * GammaDecay * y[3];
        f[4] = -Delta * y[5] - 0.5 * GammaDecay * y[4];
        f[5] = Delta * y[4] + 0.5 * Efield * (y[0]-y[1]) - 0.5 * GammaDecay * y[5];


        return GSL_SUCCESS;
}
