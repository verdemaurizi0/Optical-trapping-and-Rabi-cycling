#ifndef FEYNMAN_H
#define FEYNMAN_H

//when solving the von-neumann diffeq the positions and velocities become parameters
//thus here the structure to pass infos about the whereabouts of the molecules

int jac (double t, const double y[], double *dfdy,double dfdt[], void *cc);
int func (double t, const double y[], double f[], void *params);

extern struct laser _6um;
extern double omega0, mu;
extern double El0, w01, lambda, lambda_cb, M_PIcb, w0_sq, w0_4th, lpiw, lpiwsq, omega;

#endif /* FEYNMAN_H */
