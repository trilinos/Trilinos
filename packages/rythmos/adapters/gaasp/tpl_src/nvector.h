 #ifndef INC_nvector_h
#define INC_nvector_h

namespace CVODE {

/****************************************************************
 *                                                              *
 * File          : nvector.h                                    *
 * Programmers   : Scott D. Cohen, Alan C. Hindmarsh, and       *
 *               : Allan G. Taylor, LLNL                        *
 * Version of    : 17 December 1999                             *
 *--------------------------------------------------------------*
 *                                                              *
 * This is the header file for a generic serial NVECTOR package.*
 * It exports the type N_Vector.                                *
 *                                                              *
 * Part I of this file contains declarations which are specific *
 * to the particular machine environment in which this version  *
 * of the vector package is to be used. This includes the       *
 * typedef for the type N_Vector, as well as accessor macros    *
 * that allow the user to use efficiently the type N_Vector     *
 * without making explicit references to its underlying         *
 * representation. The underlying type of N_Vector will always  *
 * be some pointer type.                                        *
 *                                                              *
 * Part II of this file contains the prototypes for the vector  *
 * kernels which operate on the type N_Vector. These prototypes *
 * are fixed for all implementations of the vector package. The *
 * definitions of the types float and int are in the header  *
 * file llnltyps.h and these may be changed according to the    *
 * user's needs. The llnltyps.h file also contains the          *
 * definition for the type bool (short for boolean) that is the*
 * return type for the routine N_VInvTest.                      *
 *                                                              *
 * Important Note: N_Vector arguments to arithmetic kernels     *
 * need not be distinct. Thus, for example, the call            *
 *         N_VLinearSum(a,x,b,y,y);    y <- ax+by               *
 * is legal.                                                    *
 *                                                              *
 * This version of nvector.h is for the ordinary sequential     *
 * machine environment. In the documentation given below, N is  *
 * the length of all N_Vector parameters and x[i] denotes the   *
 * ith component of the N_Vector x, where 0 <= i <= N-1.        *
 *                                                              *
 ****************************************************************/


#ifdef __cplusplus     /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "llnltyps.h"


/* Part I: Machine Environment-Dependent Declarations */


/* Environment: Sequential */

 typedef struct {
  int dummy;          /* dummy element */
} *machEnvType;       /* dummy machEnvType definition */


/***************************************************************
 *                                                             *
 * Type: N_Vector                                              *
 *-------------------------------------------------------------*
 * The type N_Vector is an abstract vector type. The fields of *
 * its concrete representation should not be accessed          *
 * directly, but rather through the macros given below.        *
 *                                                             *
 * A user may assume that the N components of an N_Vector      *
 * are stored contiguously. A pointer to the first component   *
 * can be obtained via the macro N_VDATA.                      *
 *                                                             *
 ***************************************************************/

struct N_Vector_struct
{
  int length;
  float   *data;

  N_Vector_struct ()
    : length (0),
      data (0)
    {
    }
};

typedef N_Vector_struct* N_Vector;

/***************************************************************
 *                                                             *
 * Macros: N_VMAKE, N_VDISPOSE, N_VDATA, N_VLENGTH, N_VIth     *
 *-------------------------------------------------------------*
 * In the descriptions below, the following user               *
 * declarations are assumed:                                   *
 *                                                             *
 * N_Vector v; float *v_data, r; int v_len, i;              *
 *                                                             *
 * (1) N_VMAKE, N_VDISPOSE                                     *
 *                                                             *
 *     These companion routines are used to create and         *
 *     destroy an N_Vector with a component array v_data       *
 *     allocated by the user.                                  *
 *                                                             *
 *     The call N_VMAKE(v, v_data, v_len) makes v an           *
 *     N_Vector with component array v_data and length v_len.  *
 *     N_VMAKE stores the pointer v_data so that changes       *
 *     made by the user to the elements of v_data are          *
 *     simultaneously reflected in v. There is no copying of   *
 *     elements.                                               *
 *                                                             *
 *     The call N_VDISPOSE(v) frees all memory associated      *
 *     with v except for its component array. This memory was  *
 *     allocated by the user and, therefore, should be         *
 *     deallocated by the user.                                *
 *                                                             *
 * (2) N_VDATA, N_VLENGTH                                      *
 *                                                             *
 *     These routines give individual access to the parts of   *
 *     an N_Vector.                                            *
 *                                                             *
 *     The assignment v_data=N_VDATA(v) sets v_data to be      *
 *     a pointer to the first component of v. The assignment   *
 *     N_VDATA(v)=v_data sets the component array of v to      *
 *     be v_data by storing the pointer v_data.                *
 *                                                             *
 *     The assignment v_len=N_VLENGTH(v) sets v_len to be      *
 *     the length of v. The call N_VLENGTH(v)=len_v sets       *
 *     the length of v to be len_v.                            *
 *                                                             *
 * (3) N_VIth                                                  *
 *                                                             *
 *     In the following description, the components of an      *
 *     N_Vector are numbered 0..N-1, where N is the length of  *
 *     v.                                                      *
 *                                                             *
 *     The assignment r=N_VIth(v,i) sets r to be the value of  *
 *     the ith component of v. The assignment N_VIth(v,i)=r    *
 *     sets the value of the ith component of v to be r.       *
 *                                                             *
 * Notes..                                                     *
 *                                                             *
 * Users who use the macros (1) must #include<stdlib.h>        *
 * since these macros expand to calls to malloc and free.      *
 *                                                             *
 * When looping over the components of an N_Vector v, it is    *
 * more efficient to first obtain the component array via      *
 * v_data=N_VDATA(v) and then access v_data[i] within the      *
 * loop than it is to use N_VIth(v,i) within the loop.         *
 *                                                             *
 * N_VMAKE and N_VDISPOSE are similar to N_VNew and N_VFree.   *
 * The difference is one of responsibility for component       *
 * memory allocation and deallocation. N_VNew allocates memory *
 * for the N_Vector components and N_VFree frees the component *
 * memory allocated by N_VNew. For N_VMAKE and N_VDISPOSE, the *
 * component memory is allocated and freed by the user of      *
 * this package.                                               *
 *                                                             *
 ***************************************************************/

#define N_VMAKE(v, v_data, v_len) \
	v = (N_Vector) malloc(sizeof(*v)); \
	v->data   = v_data; \
	v->length = v_len

#define N_VDISPOSE(v) free(v)

#define N_VDATA(v) (v->data)

#define N_VLENGTH(v) (v->length)

#define N_VIth(v,i) ((v->data)[i])


/* Part II: N_Vector Kernel Prototypes (Machine Environment-Independent) */


/***************************************************************
 *                                                             *
 * Memory Allocation and Deallocation: N_VNew, N_VFree         *
 *                                                             *
 ***************************************************************/


/***************************************************************
 *                                                             *
 * Function : N_VNew                                           *
 * Usage    : x = N_VNew(N, machEnv);                          *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns a new N_Vector of length N. The parameter machEnv   *
 * is a pointer to machine environment-specific information.   *
 * It is ignored in the sequential machine environment and the *
 * user in this environment should simply pass NULL for this   *
 * argument. If there is not enough memory for a new N_Vector, *
 * then N_VNew returns NULL.                                   *
 *                                                             *
 ***************************************************************/

N_Vector N_VNew(int n, void *machEnv);


/***************************************************************
 *                                                             *
 * Function : N_VFree                                          *
 * Usage    : N_VFree(x);                                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Frees the N_Vector x. It is illegal to use x after the call *
 * N_VFree(x).                                                 *
 *                                                             *
 ***************************************************************/

void N_VFree(N_Vector x);


/***************************************************************
 *                                                             *
 * N_Vector Arithmetic: N_VLinearSum, N_VConst, N_VProd,       *
 *                      N_VDiv, N_VScale, N_VAbs, N_VInv,      *
 *                      N_VAddConst                            *
 *                                                             *
 ***************************************************************/


/***************************************************************
 *                                                             *
 * Function  : N_VLinearSum                                    *
 * Operation : z = a x + b y                                   *
 *                                                             *
 ***************************************************************/

void N_VLinearSum(float a, N_Vector x, float b, N_Vector y, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VConst                                        *
 * Operation : z[i] = c for i=0, 1, ..., N-1                   *
 *                                                             *
 ***************************************************************/

void N_VConst(float c, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VProd                                         *
 * Operation : z[i] = x[i] * y[i] for i=0, 1, ..., N-1         *
 *                                                             *
 ***************************************************************/

void N_VProd(N_Vector x, N_Vector y, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VDiv                                          *
 * Operation : z[i] = x[i] / y[i] for i=0, 1, ..., N-1         *
 *                                                             *
 ***************************************************************/

void N_VDiv(N_Vector x, N_Vector y, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VScale                                        *
 * Operation : z = c x                                         *
 *                                                             *
 ***************************************************************/

void N_VScale(float c, N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VAbs                                          *
 * Operation : z[i] = |x[i]|,   for i=0, 1, ..., N-1           *
 *                                                             *
 ***************************************************************/

void N_VAbs(N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VInv                                          *
 * Operation : z[i] = 1.0 / x[i] for i = 0, 1, ..., N-1        *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine does not check for division by 0. It should be *
 * called only with an N_Vector x which is guaranteed to have  *
 * all non-zero components.                                    *
 *                                                             *
 ***************************************************************/

void N_VInv(N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VAddConst                                     *
 * Operation : z[i] = x[i] + b   for i = 0, 1, ..., N-1        *
 *                                                             *
 ***************************************************************/

void N_VAddConst(N_Vector x, float b, N_Vector z);


/***************************************************************
 *                                                             *
 * N_Vector Measures: N_VDotProd, N_VMaxNorm, VWrmsNorm,       *
 *                    N_VMin, N_VWL2Norm, N_VL1Norm            *
 *                                                             *
 *                                                             *
 ***************************************************************/


/***************************************************************
 *                                                             *
 * Function : N_VDotProd                                       *
 * Usage    : dotprod = N_VDotProd(x, y);                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the value of the ordinary dot product of x and y:   *
 *                                                             *
 * -> sum (i=0 to N-1) {x[i] * y[i]}                           *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

float N_VDotProd(N_Vector x, N_Vector y);


/***************************************************************
 *                                                             *
 * Function : N_VMaxNorm                                       *
 * Usage    : maxnorm = N_VMaxNorm(x);                         *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the maximum norm of x:                              *
 *                                                             *
 * -> max (i=0 to N-1) |x[i]|                                  *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

float N_VMaxNorm(N_Vector x);


/***************************************************************
 *                                                             *
 * Function : N_VWrmsNorm                                      *
 * Usage    : wrmsnorm = N_VWrmsNorm(x, w);                    *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the weighted root mean square norm of x with        *
 * weight vector w:                                            *
 *                                                             *
 * -> sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) / N]          *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

float N_VWrmsNorm(N_Vector x, N_Vector w);


/***************************************************************
 *                                                             *
 * Function : N_VMin                                           *
 * Usage    : min = N_VMin(x);                                 *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns min x[i] if N > 0 and returns 0.0 if N <= 0.        *
 *          i                                                  *
 *                                                             *
 ***************************************************************/

float N_VMin(N_Vector x);


/***************************************************************
 *                                                             *
 * Function : N_VWL2Norm                                       *
 * Usage    : wl2norm = N_VWL2Norm(x, w);                      *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns the weighted Euclidean L2 norm of x with            *
 * weight vector w:                                            *
 *                                                             *
 * -> sqrt [(sum (i=0 to N-1) {(x[i] * w[i])^2}) ]             *
 *                                                             *
 * Returns 0.0 if N <= 0.                                      *
 *                                                             *
 ***************************************************************/

float N_VWL2Norm(N_Vector x, N_Vector w);


/***************************************************************
 *                                                             *
 * Function : N_VL1Norm                                        *
 * Usage    : l1norm = N_VL1Norm(x);                           *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns sum of ABS(x[i]) if N > 0 and returns 0.0 if N <= 0.*
 *          i                                                  *
 *
 *     i.e., calculates and returns the L1 norm of x           *
 *                                                             *
 ***************************************************************/

float N_VL1Norm(N_Vector x);


/***************************************************************
 *                                                             *
 * Miscellaneous : N_VOneMask, N_VCompare, N_VInvTest,         *
 * N_VConstrProdPos, N_VConstrMask, and N_VMinQuotient         *
 *                                                             *
 ***************************************************************/

/***************************************************************
 *                                                             *
 * Function  : N_VOneMask                                      *
 * Operation : x[i] = 1.0 if |x[i]| != 0.  i = 0, 1, ..., N-1  *
 *                    0.0 otherwise                            *
 *                                                             *
 ***************************************************************/

void N_VOneMask(N_Vector x);


/***************************************************************
 *                                                             *
 * Function  : N_VCompare                                      *
 * Operation : z[i] = 1.0 if |x[i]| >= c   i = 0, 1, ..., N-1  *
 *                    0.0 otherwise                            *
 *                                                             *
 ***************************************************************/


void N_VCompare(float c, N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function  : N_VInvTest                                      *
 * Operation : z[i] = 1.0 / x[i] with a test for x[i]==0.0     *
 *             before inverting x[i].                          *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine returns TRUE if all components of x are        *
 * non-zero (successful inversion) and returns FALSE           *
 * otherwise.                                                  *
 *                                                             *
 ***************************************************************/

bool N_VInvTest(N_Vector x, N_Vector z);


/***************************************************************
 *                                                             *
 * Function : N_VConstrProdPos                                 *
 * Usage    : booltest = N_VConstrProdPos(c,x);                *
 *-------------------------------------------------------------*
 *                                                             *
 * Returns a boolean FALSE if some c[i]!=0.0 and x[i]*c[i]<=0.0*
 *         and       TRUE otherwise                            *
 *                                                             *
 * This routine is used for constraint checking.               *
 *                                                             *
 ***************************************************************/

bool N_VConstrProdPos(N_Vector c, N_Vector x);


/***************************************************************
 *                                                             *
 * Function  : N_VConstrMask                                   *
 * Operation : m[i] = 1.0 , if  constraint test fails, for i   *
 *             m[i] = 0.0 , if  constraint test passes, for i  *
 *                where the constraint tests parallel those    *
 *                of routine N_VConstrProdPos                  *
 *-------------------------------------------------------------*
 * This routine returns a bool FALSE if any element failed    *
 * the constraint test, TRUE if all passed. It also creates a  *
 * mask vector, m, which has all elements whose corresponding  *
 * constraint test failed, marked with 1.0, passed with 0.0    *
 * This routine is specialized in that it is used only for     *
 * constraint checking.                                        *
 ***************************************************************/

bool N_VConstrMask(N_Vector c, N_Vector x, N_Vector m);


/***************************************************************
 *                                                             *
 * Function  : N_VMinQuotient                                  *
 * Operation : minq  = min ( num[i]/denom[i]) over all i such  *
 *             that   denom[i] != 0.                           *
 *-------------------------------------------------------------*
 *                                                             *
 * This routine returns the minimum of the quotients obtained  *
 * by term-wise dividing num[i] by denom[i]. A zero element    *
 * in denom will be skipped. If no such quotients are found,   *
 * then the large value 1.e99 is returned.                     *
 *                                                             *
 ***************************************************************/

float N_VMinQuotient(N_Vector num, N_Vector denom);


/***************************************************************
 *                                                             *
 * Debugging Tools : N_VPrint                                  *
 *                                                             *
 ***************************************************************/

/***************************************************************
 *                                                             *
 * Function : N_VPrint                                         *
 * Usage    : N_VPrint(x);                                     *
 *-------------------------------------------------------------*
 *                                                             *
 * Prints the N_Vector x to stdout. Each component of x is     *
 * printed on a separate line using the %g specification. This *
 * routine is provided as an aid in debugging code which uses  *
 * this vector package.                                        *
 *                                                             *
 ***************************************************************/

void N_VPrint(N_Vector x);

#ifdef __cplusplus
}
#endif

} // namespace CVODE

#endif // INC_nvector_h
