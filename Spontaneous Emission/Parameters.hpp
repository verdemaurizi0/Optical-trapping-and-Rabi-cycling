#ifndef PARAMETERS_H
#define PARAMETERS_H

// Microtraps

#define Number_Molecules 500    // Total number of CO molecules
#define Microtrap_Shape "Cyl"   // Microtrap's shape: "Cyl" or "Ell"
#define Microtrap_Y 0.          // Microtrap's y-coordinate when it leves the Chip     [µm]
#define Microtrap_Z 0.          // Microtrap's z-coordinate when it leaves the Chip     [µm]
#define Microtrap_DX 10.        // Microtrap's dx = semiaxis of an ellipsoid/cylinder   [µm]
#define Microtrap_DY 2000.      // Microtrap's dy = semiaxis of an ellipsoid/cylinder   [µm]
#define Microtrap_DZ 10.        // Microtrap's dz = semiaxis of an ellipsoid/cylinder   [µm]
#define Chip_DY    20000.       // Microchip's half width                               [µm]
#define Microtrap_V 2.          // Microtrap's velocity                                 [µm/µs]
#define Microtrap_T 0.1         // Microtrap's temperature                              [mK]

// Physical constants

#define a3pi1_T 2630.           // lifetime of molecules in the a3pi1 state [µs]

#endif /* PARAMS_H */
