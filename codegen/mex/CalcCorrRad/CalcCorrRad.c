/*
 * CalcCorrRad.c
 *
 * Code generation for function 'CalcCorrRad'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "CalcCorrRad.h"
#include "CalcCorrRad_emxutil.h"
#include "error.h"
#include "eml_int_forloop_overflow_check.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 19, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtRSInfo b_emlrtRSI = { 23, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtRSInfo c_emlrtRSI = { 27, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtRSInfo d_emlrtRSI = { 49, "power",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\ops\\power.m" };

static emlrtRSInfo e_emlrtRSI = { 12, "sqrt",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m"
};

static emlrtRSInfo f_emlrtRSI = { 39, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRSInfo g_emlrtRSI = { 61, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRSInfo h_emlrtRSI = { 108, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtRSInfo i_emlrtRSI = { 20, "eml_int_forloop_overflow_check",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\eml\\eml_int_forloop_overflow_check.m"
};

static emlrtRSInfo j_emlrtRSI = { 9, "sum",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\datafun\\sum.m"
};

static emlrtRSInfo k_emlrtRSI = { 58, "sumprod",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\datafun\\private\\sumprod.m"
};

static emlrtRSInfo l_emlrtRSI = { 69, "combine_vector_elements",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\datafun\\private\\combine_vector_elements.m"
};

static emlrtRTEInfo emlrtRTEI = { 1, 13, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtRTEInfo b_emlrtRTEI = { 15, 1, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtRTEInfo d_emlrtRTEI = { 48, 9, "sumprod",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\datafun\\private\\sumprod.m"
};

static emlrtRTEInfo e_emlrtRTEI = { 20, 15, "sumprod",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\datafun\\private\\sumprod.m"
};

static emlrtRTEInfo f_emlrtRTEI = { 71, 15, "reshape",
  "C:\\Program Files\\MATLAB\\R2016a\\toolbox\\eml\\lib\\matlab\\elmat\\reshape.m"
};

static emlrtECInfo emlrtECI = { -1, 27, 15, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtRTEInfo g_emlrtRTEI = { 18, 5, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m"
};

static emlrtBCInfo emlrtBCI = { -1, -1, 17, 12, "data", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

static emlrtDCInfo emlrtDCI = { 13, 11, "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  4 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 27, 19, "rad_bins", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 27, 39, "rad_bins", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 27, 5, "G2", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 19, 35, "data", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 19, 13, "rad", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 19, 15, "rad", "CalcCorrRad",
  "C:\\Users\\Bryce\\Dropbox\\UNI\\project\\programs\\Halo Correlations\\CalcCorrRad.m",
  0 };

/* Function Definitions */
void CalcCorrRad(const emlrtStack *sp, const emxArray_real_T *rad_bins, const
                 emxArray_real_T *data, emxArray_real_T *G2)
{
  int32_T i0;
  int32_T maxdimlen;
  emxArray_real_T *rad;
  int32_T n;
  int32_T nx;
  uint32_T varargin_1[2];
  int32_T m;
  uint32_T b_m;
  emxArray_real_T *b_rad;
  real_T b_rad_bins;
  int32_T i1;
  emxArray_boolean_T *r0;
  emxArray_boolean_T *x;
  boolean_T overflow;
  boolean_T p;
  int32_T exitg1;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* ======================================90char============================================= */
  /* +++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++ */
  /* there a number of ways to approach counting the number of points within a */
  /* spcified radius */
  /* simplest approach */
  /* calc the distance from each point to every other */
  /* prevent duplicates */
  /* then for each point calc the number within that bin using the rad */
  /* data=vertcat(halo_centered_cells{1:50}); */
  i0 = G2->size[0];
  maxdimlen = rad_bins->size[1] - 1;
  if (!(maxdimlen >= 0)) {
    emlrtNonNegativeCheckR2012b(maxdimlen, &emlrtDCI, sp);
  }

  G2->size[0] = maxdimlen;
  emxEnsureCapacity(sp, (emxArray__common *)G2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  maxdimlen = rad_bins->size[1] - 1;
  if (!(maxdimlen >= 0)) {
    emlrtNonNegativeCheckR2012b(maxdimlen, &emlrtDCI, sp);
  }

  for (i0 = 0; i0 < maxdimlen; i0++) {
    G2->data[i0] = 0.0;
  }

  emxInit_real_T(sp, &rad, 2, &b_emlrtRTEI, true);
  i0 = rad->size[0] * rad->size[1];
  rad->size[0] = data->size[0];
  rad->size[1] = data->size[0];
  emxEnsureCapacity(sp, (emxArray__common *)rad, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  maxdimlen = data->size[0] * data->size[0];
  for (i0 = 0; i0 < maxdimlen; i0++) {
    rad->data[i0] = 0.0;
  }

  for (n = 0; n < data->size[0]; n++) {
    i0 = data->size[0] * 3;
    maxdimlen = n + 1;
    if (!((maxdimlen >= 1) && (maxdimlen <= i0))) {
      emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, i0, &emlrtBCI, sp);
    }

    i0 = data->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)n) + 1.0));
    emlrtForLoopVectorCheckR2012b((1.0 + (real_T)n) + 1.0, 1.0, data->size[0],
      mxDOUBLE_CLASS, i0, &g_emlrtRTEI, sp);
    for (m = 0; m < i0; m++) {
      b_m = ((uint32_T)n + m) + 2U;
      st.site = &emlrtRSI;
      maxdimlen = data->size[0] * 3;
      nx = (int32_T)b_m;
      if (!((nx >= 1) && (nx <= maxdimlen))) {
        emlrtDynamicBoundsCheckR2012b(nx, 1, maxdimlen, &e_emlrtBCI, &st);
      }

      b_rad_bins = data->data[n] - data->data[nx - 1];
      b_st.site = &d_emlrtRSI;
      b_rad_bins *= b_rad_bins;
      st.site = &emlrtRSI;
      if (b_rad_bins < 0.0) {
        b_st.site = &e_emlrtRSI;
        error(&b_st);
      }

      maxdimlen = rad->size[0];
      nx = 1 + n;
      if (!((nx >= 1) && (nx <= maxdimlen))) {
        emlrtDynamicBoundsCheckR2012b(nx, 1, maxdimlen, &f_emlrtBCI, sp);
      }

      maxdimlen = rad->size[1];
      i1 = (int32_T)b_m;
      if (!((i1 >= 1) && (i1 <= maxdimlen))) {
        emlrtDynamicBoundsCheckR2012b(i1, 1, maxdimlen, &g_emlrtBCI, sp);
      }

      rad->data[(nx + rad->size[0] * (i1 - 1)) - 1] = muDoubleScalarSqrt
        (b_rad_bins);
    }
  }

  st.site = &b_emlrtRSI;
  nx = rad->size[0] * rad->size[1];
  b_st.site = &f_emlrtRSI;
  for (i0 = 0; i0 < 2; i0++) {
    varargin_1[i0] = (uint32_T)rad->size[i0];
  }

  maxdimlen = (int32_T)varargin_1[0];
  if ((int32_T)varargin_1[1] > (int32_T)varargin_1[0]) {
    maxdimlen = (int32_T)varargin_1[1];
  }

  maxdimlen = muIntScalarMax_sint32(nx, maxdimlen);
  if (nx > maxdimlen) {
    b_st.site = &g_emlrtRSI;
    b_error(&b_st);
  }

  if (1 > maxdimlen) {
    b_st.site = &g_emlrtRSI;
    b_error(&b_st);
  }

  emxInit_real_T1(&st, &b_rad, 1, &b_emlrtRTEI, true);
  i0 = b_rad->size[0];
  b_rad->size[0] = nx;
  emxEnsureCapacity(&st, (emxArray__common *)b_rad, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  if (nx == b_rad->size[0]) {
  } else {
    emlrtErrorWithMessageIdR2012b(&st, &f_emlrtRTEI,
      "Coder:MATLAB:getReshapeDims_notSameNumel", 0);
  }

  b_st.site = &h_emlrtRSI;
  if ((!(1 > nx)) && (nx > 2147483646)) {
    c_st.site = &i_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (maxdimlen = 0; maxdimlen + 1 <= nx; maxdimlen++) {
    b_rad->data[maxdimlen] = rad->data[maxdimlen];
  }

  emxFree_real_T(&rad);

  /* rad=rad(rad~=0); */
  n = 1;
  emxInit_boolean_T(sp, &r0, 1, &emlrtRTEI, true);
  emxInit_boolean_T(sp, &x, 1, &emlrtRTEI, true);
  while (n - 1 <= rad_bins->size[1] - 2) {
    i0 = rad_bins->size[1];
    maxdimlen = 1 + n;
    if (!((maxdimlen >= 1) && (maxdimlen <= i0))) {
      emlrtDynamicBoundsCheckR2012b(maxdimlen, 1, i0, &b_emlrtBCI, sp);
    }

    b_rad_bins = rad_bins->data[maxdimlen - 1];
    i0 = x->size[0];
    x->size[0] = b_rad->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)x, i0, (int32_T)sizeof(boolean_T),
                      &emlrtRTEI);
    maxdimlen = b_rad->size[0];
    for (i0 = 0; i0 < maxdimlen; i0++) {
      x->data[i0] = (b_rad->data[i0] < b_rad_bins);
    }

    i0 = rad_bins->size[1];
    if (!((n >= 1) && (n <= i0))) {
      emlrtDynamicBoundsCheckR2012b(n, 1, i0, &c_emlrtBCI, sp);
    }

    b_rad_bins = rad_bins->data[n - 1];
    i0 = r0->size[0];
    r0->size[0] = b_rad->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)r0, i0, (int32_T)sizeof(boolean_T),
                      &emlrtRTEI);
    maxdimlen = b_rad->size[0];
    for (i0 = 0; i0 < maxdimlen; i0++) {
      r0->data[i0] = (b_rad->data[i0] > b_rad_bins);
    }

    i0 = x->size[0];
    maxdimlen = r0->size[0];
    if (i0 != maxdimlen) {
      emlrtSizeEqCheck1DR2012b(i0, maxdimlen, &emlrtECI, sp);
    }

    st.site = &c_emlrtRSI;
    i0 = x->size[0];
    emxEnsureCapacity(&st, (emxArray__common *)x, i0, (int32_T)sizeof(boolean_T),
                      &emlrtRTEI);
    maxdimlen = x->size[0];
    for (i0 = 0; i0 < maxdimlen; i0++) {
      x->data[i0] = (x->data[i0] && r0->data[i0]);
    }

    b_st.site = &j_emlrtRSI;
    if ((x->size[0] == 1) || (x->size[0] != 1)) {
      overflow = true;
    } else {
      overflow = false;
    }

    if (overflow) {
    } else {
      emlrtErrorWithMessageIdR2012b(&b_st, &e_emlrtRTEI,
        "Coder:toolbox:autoDimIncompatibility", 0);
    }

    overflow = false;
    p = false;
    maxdimlen = 0;
    do {
      exitg1 = 0;
      if (maxdimlen < 2) {
        if (maxdimlen + 1 <= 1) {
          i0 = x->size[0];
        } else {
          i0 = 1;
        }

        if (i0 != 0) {
          exitg1 = 1;
        } else {
          maxdimlen++;
        }
      } else {
        p = true;
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    if (!p) {
    } else {
      overflow = true;
    }

    if (!overflow) {
    } else {
      emlrtErrorWithMessageIdR2012b(&b_st, &d_emlrtRTEI,
        "Coder:toolbox:UnsupportedSpecialEmpty", 0);
    }

    c_st.site = &k_emlrtRSI;
    if (x->size[0] == 0) {
      b_rad_bins = 0.0;
    } else {
      b_rad_bins = x->data[0];
      d_st.site = &l_emlrtRSI;
      overflow = ((!(2 > x->size[0])) && (x->size[0] > 2147483646));
      if (overflow) {
        e_st.site = &i_emlrtRSI;
        check_forloop_overflow_error(&e_st);
      }

      for (maxdimlen = 2; maxdimlen <= x->size[0]; maxdimlen++) {
        b_rad_bins += (real_T)x->data[maxdimlen - 1];
      }
    }

    i0 = G2->size[0];
    if (!((n >= 1) && (n <= i0))) {
      emlrtDynamicBoundsCheckR2012b(n, 1, i0, &d_emlrtBCI, sp);
    }

    G2->data[n - 1] = b_rad_bins;
    n++;
  }

  emxFree_boolean_T(&x);
  emxFree_boolean_T(&r0);
  emxFree_real_T(&b_rad);

  /*   */
  /*  for n=1:(length(rad_bins)-1) */
  /*      for m=1:size(rad,1) */
  /*          if rad(m)<rad_bins(n+1) && rad(m)>rad_bins(n) */
  /*              G2(n)= G2(n)+1; */
  /*          end */
  /*      end */
  /*  end */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (CalcCorrRad.c) */
