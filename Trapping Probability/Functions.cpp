#include <math.h>
#include <string.h>
#include "Parameters.hpp"
#include "Functions.hpp"

double SIGMA_V(const double T)
    {
    return sqrt(0.2969453*T);                       // 0.2969453 = ( kB * 10^(-3) ) / ( mCO )    
    }

double IPG_WY(const double y)
    {
    return IPG_W*sqrt(1.+( (( ( y - IPG_Y ) * IPG_L )/( M_PI * IPG_W * IPG_W ))*(( ( y - IPG_Y ) * IPG_L )/( M_PI * IPG_W * IPG_W )) ));
    }

double RDL_WY(const double y)
    {
    return RDL_W*sqrt(1.+( (( ( y - RDL_Y ) * RDL_L )/( M_PI * RDL_W * RDL_W ))*(( ( y - RDL_Y ) * RDL_L )/( M_PI * RDL_W * RDL_W )) ));
    }

double EL_0(const double v)
    {
    double y[6]={0.,0.,0.,IPG_X,IPG_Y,IPG_Z};
    return sqrt( ( ( v*v/2. + Half_L + 0.5*a3pi1_p*EL_IPG(y)*EL_IPG(y) )*( v*v/2. + Half_L + 0.5*a3pi1_p*EL_IPG(y)*EL_IPG(y) ) - Half_L * Half_L )/( a3pi1_Mu * a3pi1_Mu ) );
    }

double K0(const double v)
    {
    double IPG_W_min=10., Microtrap_T_max=10.;
    return ( EL_0( v + SIGMA_V( Microtrap_T_max ) ) - EL_0( v - SIGMA_V( Microtrap_T_max) ) )/( IPG_W_min * EL_0( v ) );
    }

double EL_STARK(const double y[])
    {
    if ( y[3] <  IPG_X )
        return EL_0( Microtrap_V ) * exp( K0( Microtrap_V ) * (y[3] - IPG_X) );
    else 
        return EL_0( Microtrap_V ) * ( 1. + K0( Microtrap_V ) * (y[3] - IPG_X) );
    }

double D_EL_STARK(const double y[])
    {
    if ( y[3] <  IPG_X )
        return EL_0( Microtrap_V ) * K0( Microtrap_V ) * exp( K0( Microtrap_V ) * (y[3] - IPG_X) );
    else
        return EL_0( Microtrap_V ) * K0( Microtrap_V );    
    }

int INSIDE_RDL(const double y[])
    {
    if ( strcmp( RDL_Shape, "Cutted" )==0)
        {
        if ( hypot( y[3]-RDL_X, y[5]-RDL_Z ) <  RDL_WY(y[4])  &&  y[3] > RDL_X )
        return 1;
        else return 0;
        }
    else if ( strcmp( RDL_Shape, "TEM00" )==0)
        {
        if ( hypot( y[3]-RDL_X, y[5]-RDL_Z ) <  RDL_WY(y[4]) )
        return 1;
        else return 0;
        }
    else return 0;
    }

int NEAR_RDL(const double y[], const int exc)
    {
    if ( exc == 1 )
        {
        if ( strcmp( RDL_Shape, "Cutted" )==0)
            {
            if ( y[3] > ( RDL_X - 1. )  &&   y[3] < ( RDL_X + 1. )   &&   y[5] < ( RDL_Z + RDL_WY(y[4]) )   &&   y[5] > ( RDL_Z - RDL_WY(y[4]) )  ) return 1;
            else return 0;
            }
        else if ( strcmp( RDL_Shape, "TEM00" )==0)
            {
            if ( y[3] > ( RDL_X - RDL_WY(y[4]) )  &&   y[3] < RDL_X   &&   y[5] < ( RDL_Z + RDL_WY(y[4]) )   &&   y[5] > ( RDL_Z - RDL_WY(y[4]) )  ) return 1;
            else return 0;
            }
        else return 0;
        }
    else return 0;
    }

double EL_IPG(const double y[])
    {
    double wsq, El0, lpiw, lpiwsq, ylpiwsq, expo;
    wsq = IPG_W * IPG_W;
    El0 = sqrt( 4.0 * IPG_P / ( M_PI * cEps0 * wsq ));
    lpiw = IPG_L / ( M_PI * wsq );
    lpiwsq = lpiw * lpiw;
    ylpiwsq = ( 1.0 + ( y[4] - IPG_Y ) * ( y[4] - IPG_Y ) * lpiwsq );
    expo = exp( (-( y[3] - IPG_X ) * ( y[3] - IPG_X ) - ( y[5] - IPG_Z ) * ( y[5] - IPG_Z ) ) / ( wsq * ylpiwsq ) );
    return El0 * expo / hypot( 1.0 , ( y[4] - IPG_Y ) * lpiw ); 
    }

double DX_EL_IPG(const double y[])
    {
    double El0, expo, wsq, lpiw, lpiwsq, ylpiwsq;
    wsq = IPG_W * IPG_W;
    El0 = sqrt( 4.0 * IPG_P / ( M_PI * cEps0 * wsq ));
    lpiw = IPG_L / ( M_PI * wsq );
    lpiwsq = lpiw * lpiw;
    ylpiwsq = ( 1.0 + ( y[4] - IPG_Y ) * ( y[4] - IPG_Y ) * lpiwsq );
    expo = exp( (-( y[3] - IPG_X ) * ( y[3] - IPG_X ) - ( y[5] - IPG_Z ) * ( y[5] - IPG_Z ) ) / ( wsq * ylpiwsq ) );
    return -( 2.0 * ( y[3] - IPG_X ) / ( wsq * ylpiwsq ) ) * El0 * expo / hypot( 1.0 , ( y[4] - IPG_Y ) * lpiw ); 

}

double DY_EL_IPG(const double y[])
    {
    double El0, expo, wsq, lpiw, lpiwsq, ylpiwsq;
    wsq = IPG_W * IPG_W;
    El0 = sqrt( 4.0 * IPG_P / ( M_PI * cEps0 * wsq ));
    lpiw = IPG_L / ( M_PI * wsq );
    lpiwsq = lpiw * lpiw;
    ylpiwsq = ( 1.0 + ( y[4] - IPG_Y ) * ( y[4] - IPG_Y ) * lpiwsq );
    expo = exp( (-( y[3] - IPG_X ) * ( y[3] - IPG_X ) - ( y[5] - IPG_Z ) * ( y[5] - IPG_Z ) ) / ( wsq * ylpiwsq ) );
    return (2.0 * (y[4] - IPG_Y) * lpiwsq * ( (y[3] - IPG_X) * (y[3] - IPG_X) + (y[5] - IPG_Z) * (y[5] - IPG_Z) ) / ( wsq * ylpiwsq * ylpiwsq ) - ( (y[4] - IPG_Y) * lpiwsq / ylpiwsq ) ) * El0 * expo / hypot( 1.0 , (y[4] - IPG_Y) * lpiw );
    }

double DZ_EL_IPG(const double y[])
    {
    double El0, expo, wsq, lpiw, lpiwsq, ylpiwsq;
    wsq = IPG_W * IPG_W;
    El0 = sqrt( 4.0 * IPG_P / ( M_PI * cEps0 * wsq ));
    lpiw = IPG_L / ( M_PI * wsq );
    lpiwsq = lpiw * lpiw;
    ylpiwsq = ( 1.0 + ( y[4] - IPG_Y ) * ( y[4] - IPG_Y ) * lpiwsq );
    expo = exp( (-( y[3] - IPG_X ) * ( y[3] - IPG_X ) - ( y[5] - IPG_Z ) * ( y[5] - IPG_Z ) ) / ( wsq * ylpiwsq ) );
    return -( 2.0 * (y[5] - IPG_Z) / ( wsq * ylpiwsq ) ) * El0 * expo / hypot( 1.0 , (y[4] - IPG_Y) * lpiw ); 
    }

int RGB_QS(const int exc)
    {
    if ( exc == 0) return ( 65536 * 0 + 256 * 0 + 255);
    if ( exc == 1) return ( 65536 * 255 + 256 * 140 + 0);
    else if ( exc == 2 ) return ( 65536 * 255 + 256 * 0 + 0);
    else return 0;
    }

int GS_TRAPPED(const double y[])
    {
    if ( hypot( y[3]-IPG_X, y[5]-IPG_Z ) <  2.*IPG_WY(y[4]) )
       return 1;
    else 
        return 0; 
    }

