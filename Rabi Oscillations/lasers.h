#ifndef LASERS_H
#define LASERS_H

struct laser
{
        double El0;         // [V/µm]
        double power;       // [mW]
        double w0;          // [µm]
        double omega;       // [MHz]
        double nu;          // [MHz]
        double lambda;      // [µm]
        double x0;          // [µm]
        double y0;          // [µm]
        double z0;          // [µm]
        double tinit;       // [µs]
        double PulseLength; // [µs]
};

#endif
