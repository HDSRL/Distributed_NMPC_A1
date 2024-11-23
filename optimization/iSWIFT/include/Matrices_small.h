// Place to put all the Input Matrices in Sparse Format

#ifndef __MATRICES_SMALL_H__
#define __MATRICES_SMALL_H__
#include "GlobalOptions.h"

realqp sigma_d = 0.0;

realqp Ppr[120] = {
0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,20,20,20,20,20,20,5.0833,0.14132,0.029636,-0.18235,3.6087,0.73325,-0.16659,3.2968,0.66988,1.3819,0.14132,0.31627,-0.95267,0.024363,-0.48215,-0.097969,-0.044437,0.87943,0.17869,-1.4578,0.029636,-0.95267,4.8063,-0.16525,3.2703,0.6645,0.17727,-3.5082,-0.71284,-0.28285,-0.18235,0.024363,-0.16525,2.3768,0.26252,-2.4833,-2.4837,0.33185,-2.2508,1.3819,3.6087,-0.48215,3.2703,0.26252,5.0577,0.3498,-0.66254,0.088521,-0.60042,-1.4578,0.73325,-0.097969,0.6645,-2.4833,0.3498,2.7713,2.643,-0.35313,2.3952,-0.28285,-0.16659,-0.044437,0.17727,-2.4837,-0.66254,2.643,2.8264,-0.706,2.3951,1.3819,3.2968,0.87943,-3.5082,0.33185,0.088521,-0.35313,-0.706,4.9448,0.66472,-1.4578,0.66988,0.17869,-0.71284,-2.2508,-0.60042,2.3952,2.3951,0.66472,2.4345,-0.28285,1.3819,-1.4578,-0.28285,20,1.3819,-1.4578,-0.28285,20,1.3819,-1.4578,-0.28285,20
};
idxint Pir[120] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,18,19,20,21,22,23,24,25,26,27,18,19,20,21,22,23,24,25,26,27,18,19,20,21,22,23,24,25,26,28,18,19,20,21,22,23,24,25,26,28,18,19,20,21,22,23,24,25,26,28,18,19,20,21,22,23,24,25,26,29,18,19,20,21,22,23,24,25,26,29,18,19,20,21,22,23,24,25,26,29,18,19,20,27,21,22,23,28,24,25,26,29
};
idxint Pjc[31] = {
0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,28,38,48,58,68,78,88,98,108,112,116,120
};
realqp Gpr[40] = {
1,-1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,-1
};
idxint Gir[40] = {
0,1,2,3,0,1,2,3,16,17,4,5,6,7,4,5,6,7,18,19,8,9,10,11,8,9,10,11,20,21,12,13,14,15,12,13,14,15,22,23
};
idxint Gjc[31] = {
0,2,4,10,12,14,20,22,24,30,32,34,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40
};
realqp Apr[66] = {
-0.0022222,0.15748,0.092595,-0.019377,-0.0022222,0.061067,-0.039223,0.071813,-0.0022222,-0.13968,0.0017738,-0.071402,-0.0022222,0.15896,0.078814,-0.020292,-0.0022222,0.061067,-0.039223,0.071813,-0.0022222,-0.12959,0.0011203,-0.058337,-0.0022222,0.15748,0.092595,-0.019377,-0.0022222,0.064245,-0.068754,0.069851,-0.0022222,-0.20266,-0.00019648,-0.061898,-0.0022222,0.15896,0.078814,-0.020292,-0.0022222,0.064245,-0.068754,0.069851,-0.0022222,-0.19258,-0.00084995,-0.048833,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
};
idxint Air[66] = {
3,15,16,17,4,15,16,17,5,15,16,17,3,15,16,17,4,15,16,17,5,15,16,17,3,15,16,17,4,15,16,17,5,15,16,17,3,15,16,17,4,15,16,17,5,15,16,17,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17
};
idxint Ajc[31] = {
0,4,8,12,16,20,24,28,32,36,40,44,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66
};
realqp c[30] = {
0.17519,0.11003,0.1245,0.11741,0.041548,0.060249,0.094185,0.046098,0.16886,0.038953,0.045184,0.034142,-14.894,-3.7791,-13.736,-3.6702,-7.3697,-12.512,0,0,0,0,0,0,0,0,0,0,0,0
};
realqp h[24] = {
10.783,12.535,11.109,12.209,10.75,11.925,11.13,11.545,11.41,12.351,11.65,12.111,11.012,11.402,10.981,11.433,0.37752,1.6225,0.69875,1.3012,0.15569,1.8443,0.82929,1.1707
};
realqp b[18] = {
0.51721,0.82141,0.80295,0.64905,0.38131,0.81589,-0.32372,0.945,0.046756,0.047061,0.065438,-0.99675,-0.94498,-0.32047,-0.065656,7.305,-0.11164,2.8261
};
idxint P[72] = {
12,30,13,31,14,32,15,33,49,48,61,60,57,56,53,52,17,35,67,66,55,54,69,68,59,58,71,70,63,62,16,34,65,64,51,50,42,43,44,29,39,40,41,28,36,37,38,27,18,19,20,21,22,23,24,25,26,1,10,7,4,5,8,0,2,3,6,9,11,45,46,47
};
idxint n = 30;
idxint m = 24;
idxint p = 18;











#endif