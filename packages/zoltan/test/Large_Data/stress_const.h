#include <zoltan.h>
#include <zz_util_const.h>

#ifndef STRESS_CONST_H
#define STRESS_CONST_H

/*
 * We need an invalid value for a ZOLTAN_ID_TYPE
 */

#undef ZOLTAN_ID_INVALID

#ifdef ZOLTAN_ID_TYPE_LONG
#define ZOLTAN_ID_INVALID LONG_MAX
#endif

#ifdef ZOLTAN_ID_TYPE_LONG_LONG
#define ZOLTAN_ID_INVALID LLONG_MAX
#endif

#ifdef ZOLTAN_ID_TYPE_UINT
#define ZOLTAN_ID_INVALID UINT_MAX
#endif

#ifdef ZOLTAN_ID_TYPE_INT
#define ZOLTAN_ID_INVALID INT_MAX
#endif

#ifndef ZOLTAN_ID_INVALID
#define ZOLTAN_ID_INVALID UINT_MAX
#endif

#endif
