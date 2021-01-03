#ifndef FUNCTIONS_H
#define FUNCTIONS_H

extern int Excitation;

double SIGMA_V(const double T);  
double IPG_WY(const double y);
double RDL_WY(const double y);

double EL_0(const double v); 
double K0(const double v);                        
double EL_STARK(const double y[]);
double D_EL_STARK(const double y[]);

int INSIDE_RDL(const double y[]);

double EL_IPG(const double y[]);
double DX_EL_IPG(const double y[]);
double DY_EL_IPG(const double y[]);
double DZ_EL_IPG(const double y[]);

int RGB_QS(const int exc);
int GS_TRAPPED(const double y[]);

#endif /* FUNCTIONS_H */



