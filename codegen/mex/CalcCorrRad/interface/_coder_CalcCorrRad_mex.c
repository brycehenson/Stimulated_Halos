/*
 * _coder_CalcCorrRad_mex.c
 *
 * Code generation for function '_coder_CalcCorrRad_mex'
 *
 */

/* Include files */
#include "CalcCorrRad.h"
#include "_coder_CalcCorrRad_mex.h"
#include "CalcCorrRad_terminate.h"
#include "_coder_CalcCorrRad_api.h"
#include "CalcCorrRad_initialize.h"
#include "CalcCorrRad_data.h"

/* Function Declarations */
static void CalcCorrRad_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[2]);

/* Function Definitions */
static void CalcCorrRad_mexFunction(int32_T nlhs, mxArray *plhs[1], int32_T nrhs,
  const mxArray *prhs[2])
{
  int32_T n;
  const mxArray *inputs[2];
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL, NULL, NULL };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        11, "CalcCorrRad");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 11,
                        "CalcCorrRad");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
  }

  /* Call the function. */
  CalcCorrRad_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  CalcCorrRad_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(CalcCorrRad_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  CalcCorrRad_initialize();

  /* Dispatch the entry-point. */
  CalcCorrRad_mexFunction(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_CalcCorrRad_mex.c) */
