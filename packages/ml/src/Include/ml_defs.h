#ifndef __MLDEFS__
#define __MLDEFS__

#ifdef matched
#define MLFORTRAN(aaa) aaa
#else
#define MLFORTRAN(aaa) aaa ## _
#endif

#define ML_VERSION        ml2_0_0_5

#ifndef MB_MODIF
#define MB_MODIF
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
#define ML_INTERNAL         111
#define ML_EXTERNAL         112
#define ML_CONVERGE          -2
#define ML_NOTSET            -1

#define ML_CHAR               1
#define ML_INT                2
#define ML_DOUBLE             3

#define ML_ALL_LEVELS     -1237
#define ML_MSR_MATRIX      -201
#define ML_CSR_MATRIX      -203

/* JJH -- these are message tags */
#define ML_TAG_BASE        27000
#define ML_TAG_PRESM       1
#define ML_TAG_POSTSM      101

/* maximum dimension of the subspace associated with an ML_Operator type */
#define ML_MAX_SUBSPACE_DIM 3

#endif
