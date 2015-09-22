/****************************************************************
 *                                                              *
 * File          : nvector.c                                    *
 * Programmers   : Scott D. Cohen,  Alan C. Hindmarsh, and      *
 *               : Allan G. Taylor, LLNL                        *
 * Version of    : 17 December 1999                             *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the implementation file for a generic serial NVECTOR *
 * package. It contains the implementation of the N_Vector      *
 * kernels listed in nvector.h.                                 *
 *                                                              *
 ****************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include "nvector.h"
#include "llnltyps.h"
#include "llnlmath.h"

namespace CVODE {

#define ZERO RCONST(0.0)
#define HALF RCONST(0.5)
#define ONE  RCONST(1.0)
#define ONEPT5 RCONST(1.5)


/* Private Helper Prototypes */

static void VCopy(N_Vector x, N_Vector z); /* z=x */
static void VSum(N_Vector x, N_Vector y, N_Vector z); /* z=x+y */
static void VDiff(N_Vector x, N_Vector y, N_Vector z); /* z=x-y */
static void VNeg(N_Vector x, N_Vector z); /* z=-x */
/* z=c(x+y) */
static void VScaleSum(float c, N_Vector x, N_Vector y, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff(float c, N_Vector x, N_Vector y, N_Vector z);
static void VLin1(float a, N_Vector x, N_Vector y, N_Vector z); /* z=ax+y */
static void VLin2(float a, N_Vector x, N_Vector y, N_Vector z); /* z=ax-y */
static void Vaxpy(float a, N_Vector x, N_Vector y); /* y <- ax+y */
static void VScaleBy(float a, N_Vector x); /* x <- ax */

/********************* Exported Functions ************************/


N_Vector N_VNew(int N, void *machEnv)
{
  N_Vector v;

  if (N <= 0) return(NULL);

  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  v->data = (float *) malloc(N * sizeof(float));
  if (v->data == NULL) {
    free(v);
    return(NULL);
  }

  v->length = N;

  return(v);
}


void N_VFree(N_Vector x)
{
  free(x->data);
  free(x);
}


void N_VLinearSum(float a, N_Vector x, float b, N_Vector y, N_Vector z)
{
  int i, N;
  float c, *xd, *yd, *zd;
  N_Vector v1, v2;
  bool test;

  if ((b == ONE) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
    Vaxpy(a,x,y);
    return;
  }

  if ((a == ONE) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
    Vaxpy(b,y,x);
    return;
  }

  /* Case: a == b == 1.0 */

  if ((a == ONE) && (b == ONE)) {
    VSum(x, y, z);
    return;
  }

  /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */

  if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
    v1 = test ? y : x;
    v2 = test ? x : y;
    VDiff(v2, v1, z);
    return;
  }

  /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
  /* if a or b is 0.0, then user should have called N_VScale */

  if ((test = (a == ONE)) || (b == ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin1(c, v1, v2, z);
    return;
  }

  /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */

  if ((test = (a == -ONE)) || (b == -ONE)) {
    c = test ? b : a;
    v1 = test ? y : x;
    v2 = test ? x : y;
    VLin2(c, v1, v2, z);
    return;
  }

  /* Case: a == b */
  /* catches case both a and b are 0.0 - user should have called N_VConst */

  if (a == b) {
    VScaleSum(a, x, y, z);
    return;
  }

  /* Case: a == -b */

  if (a == -b) {
    VScaleDiff(a, x, y, z);
    return;
  }

  /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = a * (*xd++) + b * (*yd++);
}


void N_VConst(float c, N_Vector z)
{
  int i, N;
  float *zd;

  N = z->length;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = c;
}


void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = (*xd++) * (*yd++);
}


void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = (*xd++) / (*yd++);
}


void N_VScale(float c, N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  if (z == x) {       /* BLAS usage: scale x <- cx */
    VScaleBy(c, x);
    return;
  }

  if (c == ONE) {
    VCopy(x, z);
  } else if (c == -ONE) {
    VNeg(x, z);
  } else {
    N = x->length;
    xd = x->data;
    zd = z->data;
    for (i=0; i < N; ++i ) *zd++ = c * (*xd++);
  }
}


void N_VAbs(N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; i++, xd++, zd++)
    *zd = ABS(*xd);
}


void N_VInv(N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = ONE / (*xd++);
}


void N_VAddConst(N_Vector x, float b, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; ++i ) *zd++ = (*xd++) + b;
}


float N_VDotProd(N_Vector x, N_Vector y)
{
  int i, N;
  float sum = ZERO, *xd, *yd;

  N = x->length;
  xd = x->data;
  yd = y->data;

  for (i=0; i < N; ++i )
    sum += (*xd++) * (*yd++);

  return(sum);
}


float N_VMaxNorm(N_Vector x)
{
  int i, N;
  float max = ZERO, *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i < N; i++, xd++) {
    if (ABS(*xd) > max) max = ABS(*xd);
  }

  return(max);
}


float N_VWrmsNorm(N_Vector x, N_Vector w)
{
  int i, N;
  float sum = ZERO, prodi, *xd, *wd;

  N = x->length;
  xd = x->data;
  wd = w->data;

  for (i=0; i < N; ++i ) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  return(RSqrt(sum / N));
}



float N_VMin(N_Vector x)
{
  int i, N;
  float min, *xd;

  N = x->length;
  xd = x->data;
  min = xd[0];

  for (i=1; i < N; i++, xd++) {
    if ((*xd) < min) min = *xd;
  }

  return(min);
}


float N_VWL2Norm(N_Vector x, N_Vector w)
{
  int i, N;
  float sum = ZERO, prodi, *xd, *wd;

  N = x->length;
  xd = x->data;
  wd = w->data;

  for (i=0; i < N; ++i ) {
    prodi = (*xd++) * (*wd++);
    sum += prodi * prodi;
  }

  return(RSqrt(sum));
}


float N_VL1Norm(N_Vector x)
{
 int i, N;
 float sum = ZERO, *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i<N; ++i )
   sum += ABS(xd[i]);

 return(sum);
}


void N_VOneMask(N_Vector x)
{
  int i, N;
  float *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i<N; i++,xd++) {
    if (*xd != ZERO) *xd = ONE;
  }
}


void N_VCompare(float c, N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; i++, xd++, zd++) {
    *zd = (ABS(*xd) >= c) ? ONE : ZERO;
  }
}


bool N_VInvTest(N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; ++i ) {
    if (*xd == ZERO) return(FALSE);
    *zd++ = ONE / (*xd++);
  }

  return(TRUE);
}


bool N_VConstrProdPos(N_Vector c, N_Vector x)
{
  int i, N;
  float  *xd, *cd;
  bool test;

  N = x->length;
  xd = x->data;
  cd = c->data;
  test = TRUE;

  for (i=0; i < N; i++, xd++,cd++) {
    if (*cd != ZERO) {
      if ((*xd)*(*cd) <= ZERO) {
	test = FALSE;
        break;
      }
    }
  }

  return(test);
}

bool N_VConstrMask(N_Vector c, N_Vector x, N_Vector m)
{
  int i, N;
  bool test;
  float *cd, *xd, *md;

  N = x->length;
  cd = c->data;
  xd = x->data;
  md = m->data;

  test = TRUE;

  for (i=0; i<N; i++, cd++, xd++, md++) {
    if ( *cd == ZERO) *md = ZERO;
    else {
      if ( *cd > ONEPT5 || (*cd) < -ONEPT5)
        if ( (*xd)*(*cd) <= ZERO) {
	   test = FALSE;
	    *md = ONE;
        }
        else
	  *md = ZERO;
      else if ( (*cd) > HALF || (*cd) < -HALF)
        if ( (*xd)*(*cd) < ZERO ) {
	  test = FALSE;
	  *md = ONE;
	}
       else
         *md = ZERO;
    }
  }
return(test);
}


float N_VMinQuotient(N_Vector num, N_Vector denom)
{
  bool notEvenOnce;
  int i, N;
  float *nd, *dd, min;

  N = num->length;
  nd = num->data;
  dd = denom->data;
  notEvenOnce = TRUE;

  for (i=0; i<N; i++, nd++, dd++) {
    if (*dd == ZERO) continue;
    else {
      if (notEvenOnce) {
	min = *nd / *dd ;
	notEvenOnce = FALSE;
      }
      else
	min = MIN(min, (*nd)/(*dd));
    }
  }
  if (notEvenOnce) min = 1.e99;

  return(min);
}


void N_VPrint(N_Vector x)
{
  int i, N;
  float *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i < N; ++i ) printf("%11.8g\n", *xd++);

  printf("\n");
}


/***************** Private Helper Functions **********************/


static void VCopy(N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = *xd++;
}


static void VSum(N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = (*xd++) + (*yd++);
}


static void VDiff(N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = (*xd++) - (*yd++);
}


static void VNeg(N_Vector x, N_Vector z)
{
  int i, N;
  float *xd, *zd;

  N = x->length;
  xd = x->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = -(*xd++);
}


static void VScaleSum(float c, N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = c * ((*xd++) + (*yd++));
}


static void VScaleDiff(float c, N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = c * ((*xd++) - (*yd++));
}


static void VLin1(float a, N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = a * (*xd++) + (*yd++);
}


static void VLin2(float a, N_Vector x, N_Vector y, N_Vector z)
{
  int i, N;
  float *xd, *yd, *zd;

  N = x->length;
  xd = x->data;
  yd = y->data;
  zd = z->data;

  for (i=0; i < N; ++i )
    *zd++ = a * (*xd++) - (*yd++);
}

static void Vaxpy(float a, N_Vector x, N_Vector y)
{
  int i, N;
  float *xd, *yd;

  N = x->length;
  xd = x->data;
  yd = y->data;

  if (a == ONE) {
    for (i=0; i < N; ++i )
      *yd++ += (*xd++);
    return;
  }

  if (a == -ONE) {
    for (i=0; i < N; ++i )
      *yd++ -= (*xd++);
    return;
  }

  for (i=0; i < N; ++i )
    *yd++ += a * (*xd++);
}

static void VScaleBy(float a, N_Vector x)
{
  int i, N;
  float *xd;

  N = x->length;
  xd = x->data;

  for (i=0; i < N; ++i )
    *xd++ *= a;
}

} // namespace CVODE

