#include <math.h>
#include <string.h>

double SIGMA_V(const double T)
{
        return sqrt(0.2969453*T);                   // 0.2969453 = ( kB * 10^(-3) ) / ( mCO )
}
