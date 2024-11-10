/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef __MLDEFS__
#define __MLDEFS__

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

#include "ml_common.h"

/*
#ifdef ML_MATCHED
#define MLFORTRAN(aaa) aaa
#else
#define MLFORTRAN(aaa) aaa ## _
#endif
*/

#define ML_VERSION        ml4_0

#define ML_LOCATION  __FILE__,"(", __LINE__,"): "

#ifndef MB_MODIF
#define MB_MODIF
#endif
#ifndef MB_MODIF_QR
#define MB_MODIF_QR
#endif
#ifndef RST_MODIF
#define RST_MODIF
#endif
#define ML_NONE           10
#define ML_MGV            11
#define ML_MG2CGC         12
#define ML_MGFULLV        13
#define ML_RSAMG          14
#define ML_SAAMG          15
#define ML_MGW            16
#define ML_PAMGV          17
/*MS*/
/* those are to mimic one-level and two-level DD methods.
   They are a bit experimental now */
#define ML_ONE_LEVEL_DD            -30
#define ML_TWO_LEVEL_DD_ADD        -31
#define ML_TWO_LEVEL_DD_HYBRID     -32
#define ML_TWO_LEVEL_DD_HYBRID_2   -33
/*ms*/

#define ML_GRID_DIMENSION   21
#define ML_GRID_NVERTICES   22
#define ML_GRID_NELEMENTS   23
#define ML_GRID_ELEM_GLOBAL 24
#define ML_GRID_ELEM_NVERT  25
#define ML_GRID_ELEM_VLIST  26
#define ML_GRID_VERT_GLOBAL 27
#define ML_GRID_VERT_COORD  28
#define ML_GRID_BASISFCNS   29
#define ML_GRID_ELEM_VOLUME 30
#define ML_GRID_ELEM_MATRIX 31

#define ML_ID_ML            101
#define ML_ID_SL            102
#define ML_ID_SLSOL         103
#define ML_ID_KLUDGE        104
#define ML_ID_VEC           105
#define ML_ID_OPAGX         106
#define ML_ID_ILIST         107
#define ML_ID_COMM          108
#define ML_ID_COMMINFOAGX   109
#define ML_ID_COMMINFOOP    110
#define ML_ID_GRIDAGX       111
#define ML_ID_GRIDFCN       112
#define ML_ID_GGRAPH        113
#define ML_ID_MATRIX        114
#define ML_ID_DESTROYED     115
#define ML_ID_SOLVER        116
#define ML_ID_MAPPER        117
#define ML_ID_MATCSR        118
#define ML_ID_BC            119
#define ML_ID_GRID          120
#define ML_ID_SMOOTHER      121
#define ML_ID_CSOLVE        122
#define ML_ID_OP            123
#define ML_ID_MATRIXDCSR    124
#define ML_ID_AGGRE         125
#define ML_ID_KRYLOVDATA    126
#define ML_ID_AMG           127

#define ML_COMPUTE_RES_NORM 129
#define ML_NO_RES_NORM      179

#define ML_INCREASING       717
#define ML_DECREASING       718

#define ML_PRESMOOTHER      201
#define ML_POSTSMOOTHER     202
#define ML_BOTH             203

#define ML_BDRY_DIRICHLET   1
#define ML_BDRY_NEUMANN     2
#define ML_BDRY_INSIDE      'I'
#define ML_BDRY_RIDGE       'R'
#define ML_BDRY_FACE        'F'
#define ML_BDRY_CORNER      'C'

#define ML_FALSE              0
#define ML_TRUE               1

#ifndef ML_MAX_NEIGHBORS
#define ML_MAX_NEIGHBORS    250
#endif
#ifndef ML_MAX_MSG_BUFF_SIZE
#define ML_MAX_MSG_BUFF_SIZE 100000  /* max usable message buffer size */
#endif

#define ML_OVERWRITE          0
#define ML_ADD                1
#define ML_NO                 0
#define ML_YES                1
#define ML_Set                111

#define ML_EMPTY             -1
#define ML_DEFAULT           -2
#define ML_DDEFAULT          -2.0
#define ML_ONE_STEP_CG       -100
#define ML_ZERO               3
#define ML_NONZERO            4
#define ML_CONVERGE          -2
#define ML_NOTSET            -1
#define ML_NONEMPTY         111

#define ML_CHAR               1
#define ML_INT                2
#define ML_DOUBLE             3

#define ML_ALL_LEVELS     -1237
#define ML_MSR_MATRIX      -201
#define ML_CSR_MATRIX      -203
#define ML_EpetraCRS_MATRIX  -205

/* JJH -- these are message tags */
#define ML_TAG_BASE        27000
#define ML_TAG_PRESM       1
#define ML_TAG_POSTSM      101

/*MS*/
/* those are for the aggregate ordering within METIS or ParMETIS */
/*IK Also used for columns_loc_glob in the get row structure*/
#define ML_LOCAL_INDICES     0
#define ML_GLOBAL_INDICES    1
/*ms*/

/* RST: These are for doing block diagonal scaling of a matrix */
#define ROW_SCALE_WITH_D     1
#define ROW_SCALE_WITH_DT    2
#define COL_SCALE_WITH_D     3
#define COL_SCALE_WITH_DT    4
/* RST */


/* maximum dimension of the subspace associated with an ML_Operator type */
#define ML_MAX_SUBSPACE_DIM 3

/* block partitioning option to use Parmetis */
enum ml_partitioner_enum {ML_USEPARMETIS,ML_USEMETIS,ML_USEZOLTAN,ML_USEJOSTLE};
typedef enum ml_partitioner_enum ML_Partitioner;



/* Allow 64-bit integer support.  This is needed when the ML library is used
 * with MPSalsa.  MPSalsa requires 64-bit global element numbers. */
#ifndef ML_BIG_INT
#define ML_BIG_INT int
#endif
typedef ML_BIG_INT ml_big_int;

/* Eigenvalue options */
#define ML_DIAGSCALE           1
#define ML_NO_SCALE            0
#define ML_SYMMETRIC           1
#define ML_NONSYMM             3
#define ML_USE_MATRIX_NORM     9
#define ML_USE_CG              10
#define ML_USE_ANASAZI         20
#define ML_USE_POWER           30
#define ML_NO_SYMMETRIZE       0

#define ML_Get_MySmootherData(smoother_obj) ((smoother_obj)->smoother->data)
#define ML_Get_MyGetrowData(matrix_obj) ((matrix_obj)->data)
#define ML_Get_MyMatvecData(matrix_obj) ((matrix_obj)->data)

#ifdef HAVE_ML_MATLAB
#include "mex.h"
#endif /*ifdef HAVE_ML_MATLAB*/


#endif
