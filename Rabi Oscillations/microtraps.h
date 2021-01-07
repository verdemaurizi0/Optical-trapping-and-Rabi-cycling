#ifndef MICROTRAPS_H
#define MICROTRAPS_H

struct microtrap
{
        double num;       // Total number of CO molecules
        double t;         // Time creation
        double x0;        // Microtrap's x-coordinate when it leaves the Chip     [µm]
        double y0;        // Microtrap's y-coordinate when it leaves the Chip     [µm]
        double z0;        // Microtrap's z-coordinate when it leaves the Chip     [µm]
        double dx;        // Microtrap's dx = semiaxis of an ellipsoid/cylinder   [µm]
        double dy;        // Microtrap's dy = semiaxis of an ellipsoid/cylinder   [µm]  
        double dz;        // Microtrap's dz = semiaxis of an ellipsoid/cylinder   [µm]  
        double halfwidth; // Microchip's half width                               [µm]
        double vel;       // Microtrap's velocity                                 [µm/µs]
        double temp;      // Microtrap's temperature                              [mK]
}; 

#endif
