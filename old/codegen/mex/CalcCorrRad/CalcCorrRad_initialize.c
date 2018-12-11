/*
 * CalcCorrRad_initialize.c
 *
 * Code generation for function 'CalcCorrRad_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CalcCorrRad.h"
#include "CalcCorrRad_initialize.h"
#include "_coder_CalcCorrRad_mex.h"
#include "CalcCorrRad_data.h"

/* Function Definitions */
void CalcCorrRad_initialize(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (CalcCorrRad_initialize.c) */
