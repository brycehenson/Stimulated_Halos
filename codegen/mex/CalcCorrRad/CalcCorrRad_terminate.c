/*
 * CalcCorrRad_terminate.c
 *
 * Code generation for function 'CalcCorrRad_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CalcCorrRad.h"
#include "CalcCorrRad_terminate.h"
#include "_coder_CalcCorrRad_mex.h"
#include "CalcCorrRad_data.h"

/* Function Definitions */
void CalcCorrRad_atexit(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void CalcCorrRad_terminate(void)
{
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (CalcCorrRad_terminate.c) */
