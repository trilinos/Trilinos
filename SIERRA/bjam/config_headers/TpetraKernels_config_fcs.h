#ifndef TPETRAKERNELS_CONFIG_H
#define TPETRAKERNELS_CONFIG_H

/* Define if building in debug mode */
/* #undef HAVE_TPETRAKERNELS_DEBUG */

/* Define this macro if the quadmath TPL is enabled */
/* #undef HAVE_TPETRAKERNELS_QUADMATH */

/*
 * "Optimization level" for computational kernels in this subpackage.
 * The higher the level, the more code variants get generated, and
 * thus the longer the compile times.  However, more code variants
 * mean both better performance overall, and more uniform performance
 * for corner cases.
 */
#define KOKKOSLINALG_OPT_LEVEL 1

#endif // TPETRAKERNELS_CONFIG_H
