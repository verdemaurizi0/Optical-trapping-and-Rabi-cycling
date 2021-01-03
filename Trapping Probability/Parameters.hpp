#ifndef PARAMETERS_H
#define PARAMETERS_H

// Microtraps

#define Number_Molecules 5000    // Total number of CO molecules
#define Microtrap_Shape "Cyl"   // Microtrap's shape: "Cyl" or "Ell"
#define Microtrap_Y 0.          // Microtrap's y-coordinate when it leaves the Chip     [µm]
#define Microtrap_Z 0.          // Microtrap's z-coordinate when it leaves the Chip     [µm]
#define Microtrap_DX 10.        // Microtrap's dx = semiaxis of an ellipsoid/cylinder   [µm]
#define Microtrap_DY 2000.      // Microtrap's dy = semiaxis of an ellipsoid/cylinder   [µm]
#define Microtrap_DZ 10.        // Microtrap's dz = semiaxis of an ellipsoid/cylinder   [µm]
#define Chip_DY    20000.       // Microchip's half width                               [µm]
#define Microtrap_V 20.0         // Microtrap's velocity                                 [µm/µs]
#define Microtrap_T 1.0         // Microtrap's temperature                              [mK]

// IPG optical trapping laser (assumed to be oriented along the y-axis)

#define IPG_X 0.                // Focal spot x-coordinate [µm]
#define IPG_Y 0.                // Focal spot y-coordinate [µm]
#define IPG_Z 0.                // Focal spot z-coordinate [µm]
#define IPG_P 240.              // Maximum available power [W]
#define IPG_L 1.07              // Wavelength              [µm]
#define IPG_W 30.               // Waist                   [µm]          

// RDL driving quantum transitions (assumed to be oriented along the y-axis)   

#define RDL_X 0.               // Focal spot x-coordinate [µm]
#define RDL_Y 0.                // Focal spot y-coordinate [µm]
#define RDL_Z 0.                // Focal spot z-coordinate [µm]
#define RDL_L 0.563             // Wavelength              [µm]
#define RDL_W 30.               // Waist                   [µm]
#define RDL_Shape "Cutted"      // Two shapes: "TEM00" or "Cutted" 

// Physical constants

#define cEps0 0.0026544187      // c*eps0  [J / (V^2*s)]
#define a3pi1_T 2630.           // lifetime of molecules in the a3pi1 state [µs]
#define Trans_P 1.            // transit probability to the ground state gs
#define Half_L 2.80747          // Λ/(2*mCO) = 2.80747  [µm^2/µs^2]                                 (a3pi1 DC Stark)
#define a3pi1_Mu 49.3278        // µ/(2*mCO) = 49.1431  [µm^3/(µs^2*V)] (1.37 Debye -> 49.3278)     (a3pi1 DC Stark)
#define a3pi1_p 0.0047          // α/mCO = 0.0047       [um^4/(us^2*V^2)]                           (a3pi1 polarizability)
#define k_gs 0.014              // kdc/mCO = 0.014      [µm^4/(µs^2*V^2)]                           (gs DC Stark)
#define gs_p 0.0047             // α/mCO = 0.0047       [um^4/(us^2*V^2)]                           (gs polarizability)

#endif /* PARAMS_H */
