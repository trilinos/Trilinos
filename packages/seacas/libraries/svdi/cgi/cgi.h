/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* cgi.h - header file to define useful stuff for cgi
 * 27 Oct 1989 - last date modified
 * Pat McGee, jpm@lanl.gov
 */
/*
 * 26 Oct 1989,TLC: brought up to date with Pat's data structure changes
 * 27 Oct 1989,TLC: changed function id to _FN for UNICOS problems
 */

/******************************************************************************/
/*                                                                            */
/*      type and constant declarations                                        */
/*                                                                            */
/******************************************************************************/

#ifndef CGI_H
#define CGI_H

#include "fortyp.h"
#include "stdtyp.h"

/* define true and false */
#define FALSE 0
#define TRUE 1

/* macros for min and max */
#define max(A, B) ((A) > (B) ? (A) : (B))
#define min(A, B) ((A) < (B) ? (A) : (B))

#define MAX_DEVICES 3
#define MAX_IN_PARAMS 20
#define MAX_SURFACES 4
#define MAX_ERROR_CLASS 9

/* function identifiers */
#define ACTIVATE_FN -1
#define DEACTIVATE_FN -2
#define CI_FN 1
#define CT_FN 2
#define CXDFAC_FN 3
#define CPDS_FN 4
#define CENDPG_FN 5
#define CBC_FN 6
#define CVDCX_FN 7
#define CV_FN 8
#define CCL_FN 9
#define CDSCL_FN 10
#define CDQERR_FN 11
#define CERHCT_FN 12
#define CCIXP_FN 13
#define CESC_FN 14
#define CQID_FN 15
#define CQD_FN 16
#define CLF_FN 17
#define CLPR_FN 18
#define CQSP_FN 19
#define CLESC_FN 20
#define CQP_FN 21
#define CQCL_FN 22
#define CPL_FN 23
#define CDJPL_FN 24
#define CPM_FN 25
#define CTX_FN 26
#define CPG_FN 27
#define CCA_FN 28
#define CPXA_FN 29
#define CLNT_FN 30
#define CLNW_FN 31
#define CLNC_FN 32
#define CMKT_FN 33
#define CMKS_FN 34
#define CMKC_FN 35
#define CTXP_FN 36
#define CTXC_FN 37
#define CCHH_FN 38
#define CCHO_FN 39
#define CIS_FN 40
#define CFLC_FN 41
#define CCSM_FN 42
#define CCT_FN 43
#define CGTXX_FN 44
#define CQPRL_FN 45
#define CQLN_FN 46
#define CQLNT_FN 47
#define CQSLW_FN 48
#define CQMK_FN 49
#define CQMKT_FN 50
#define CQSMS_FN 51
#define CQCHH_FN 52
#define CQFL_FN 53
#define CQC_FN 54
#define CQLNA_FN 55
#define CQMKA_FN 56
#define CQTXA_FN 57
#define CQFLA_FN 58
#define CQCTE_FN 59
#define CILI_FN 60
#define CRQLC_FN 61
/* end of function identifiers */

#define MAX_FN_ID 61

/* escape identifiers */
/* ..moved these to cgidef.h */

/* device: contains a list of pointers to active surface defn's; for use by
 * multiple CGI device routines (CGIACT, CGIDAC, CGION, CGIOFF), not by the
 * regular CGI routines (CI, CT, etc.)
 */
typedef struct
{
  void (*device_fn)();           /* which device is this list for */
  short     num_active_surfaces; /* how many active */
  short     num_on_surfaces;     /* how many have output on */
  anything *statelist[MAX_SURFACES];
  /* within surfaces, all surfaces which have output on are */
  /* found at the beginning of the list, followed by all */
  /* surfaces which have output off.  num_on_surfaces */
  /* counts the number of surfaces which have output on, and */
  /* therefore points to the first surface which has output off.*/
} device_struct;

typedef struct
{                      /* error list elements */
  f_integer err_class; /* error class */
  f_integer err_num;   /* error number */
  f_integer func_id;   /* function id */
} error_report;

typedef struct
{
  float x;
  float y;
} point;

typedef struct
{
  int x;
  int y;
} dpoint;

typedef struct
{
  utiny r;
  utiny g;
  utiny b;
} rgb;

typedef int cenum; /* CGI enumeration type */

typedef enum csm_enum { dc_csm, ic_csm } csm_enum;

void cgi_def_ini();

/* end cgi.h */
#endif
