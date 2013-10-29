#ifndef DIAGNOSTIC_H
#define DIAGNOSTIC_H

#include "calibration.h"
#include "parameters.h"
#include <stdio.h>

#define DIAG_PRINT_F(fname) diag_print_ ## fname


// print to screen functions for structs in calibration.h

void DIAG_PRINT_F(Exterior) (Exterior * Ex);
//void diag_print_Exterior (Exterior * Ex);
void DIAG_PRINT_F(Interior) (Interior * I);
void DIAG_PRINT_F(Glass) (Glass * G);

//print to screen functions for structs in parameters.h

void DIAG_PRINT_F(mm_np)(mm_np * mm);


#endif
