/*
 * CalcCorrRad.h
 *
 * Code generation for function 'CalcCorrRad'
 *
 */

#ifndef CALCCORRRAD_H
#define CALCCORRRAD_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "CalcCorrRad_types.h"

/* Function Declarations */
extern void CalcCorrRad(const emlrtStack *sp, const emxArray_real_T *rad_bins,
  const emxArray_real_T *data, emxArray_real_T *G2);

#endif

/* End of code generation (CalcCorrRad.h) */
