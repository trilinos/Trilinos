/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* sdcgi - standard single device cgi routines
 * 27 Oct 1989 - last date modified
 * %Z% %M% version %I% of %G% %U%
 * Pat McGee, jpm@lanl.gov
 */
/* History:
 * 27 Oct 1989,TLC: changed function id names to _FN
 */

#include "sdcgi.h"
#include "cgi.h"
#include "fortyp.h"
#include "stdtyp.h"
#include <stdio.h>

/* ifdefc.h - ifdef file for cgi routines
 * This file is used to define the system dependent ways C is
 * called from FORTRAN.  Underscores are used by default.
 *
 * SUN DEC/ULTRIX ALLIANT : C routines must have underscores
 * SGI CONVEX             : C routines must have underscores
 *
 * VAX HP IBM/aix         : C routines do not have underscores
 *
 * CRAY/UNICOS            : C routines must be capitalized,
 *                            and no underscores
 *
 * This file also defines the system dependent macro "f2cchar".
 *
 * This include file is used by SDCGI.C and VDICGI.C
 */

#if defined(ADDC_)
#endif
#if !defined(CRA) && !defined(ADDC_) && !defined(COUGAR)
#define ci_ ci
#define ct_ ct
#define cxdfac_ cxdfac
#define cpds_ cpds
#define cendpg_ cendpg
#define cbc_ cbc
#define cvdcx_ cvdcx
#define cv_ cv
#define ccl_ ccl
#define cdqerr_ cdqerr
#define cerhct_ cerhct
#define ccixp_ ccixp
#define cesc1_ cesc1
#define cesc2_ cesc2
#define cqid_ cqid
#define cqd_ cqd
#define clf_ clf
#define clpr_ clpr
#define cqsp_ cqsp
#define clesc_ clesc
#define cqp_ cqp
#define cpl_ cpl
#define cdjpl_ cdjpl
#define cdscl_ cdscl
#define cqcl_ cqcl
#define cpm_ cpm
#define ctx1_ ctx1
#define ctx2_ ctx2
#define cpg_ cpg
#define cca_ cca
#define clnt_ clnt
#define clnw_ clnw
#define clnc_ clnc
#define cmkt_ cmkt
#define cmkc_ cmkc
#define ctxp_ ctxp
#define ctxc_ ctxc
#define cchh_ cchh
#define ccho_ ccho
#define cis_ cis
#define cflc_ cflc
#define ccsm_ ccsm
#define cct_ cct
#define cgtxx1_ cgtxx1
#define cgtxx2_ cgtxx2
#define cqprl_ cqprl
#define cqln_ cqln
#define cqlnt_ cqlnt
#define cqchh1_ cqchh1
#define cqchh2_ cqchh2
#define cqfl_ cqfl
#define cqc_ cqc
#define cqlna_ cqlna
#define cqtxa_ cqtxa
#define cqcte_ cqcte
#define cili_ cili
#define crqlc_ crqlc
#define cpxa_ cpxa
#endif
#if defined(CRA)
#define ci_ CI
#define ct_ CT
#define cxdfac_ CXDFAC
#define cpds_ CPDS
#define cendpg_ CENDPG
#define cbc_ CBC
#define cvdcx_ CVDCX
#define cv_ CV
#define ccl_ CCL
#define cdqerr_ CDQERR
#define cerhct_ CERHCT
#define ccixp_ CCIXP
#define cesc1_ CESC1
#define cesc2_ CESC2
#define cqid_ CQID
#define cqd_ CQD
#define clf_ CLF
#define clpr_ CLPR
#define cqsp_ CQSP
#define clesc_ CLESC
#define cqp_ CQP
#define cpl_ CPL
#define cdjpl_ CDJPL
#define cdscl_ CDSCL
#define cqcl_ CQCL
#define cpm_ CPM
#define ctx1_ CTX1
#define ctx2_ CTX2
#define cpg_ CPG
#define cca_ CCA
#define clnt_ CLNT
#define clnw_ CLNW
#define clnc_ CLNC
#define cmkt_ CMKT
#define cmkc_ CMKC
#define ctxp_ CTXP
#define ctxc_ CTXC
#define cchh_ CCHH
#define ccho_ CCHO
#define cis_ CIS
#define cflc_ CFLC
#define ccsm_ CCSM
#define cct_ CCT
#define cgtxx1_ CGTXX1
#define cgtxx2_ CGTXX2
#define cqprl_ CQPRL
#define cqln_ CQLN
#define cqlnt_ CQLNT
#define cqchh1_ CQCHH1
#define cqchh2_ CQCHH2
#define cqfl_ CQFL
#define cqc_ CQC
#define cqlna_ CQLNA
#define cqtxa_ CQTXA
#define cqcte_ CQCTE
#define cili_ CILI
#define crqlc_ CRQLC
#define cpxa_ CPXA
#endif
#ifdef Build64
#if defined(ADDC_)
#define cesc1_ cesc14_
#define ctx1_ ctx14_
#define cgtxx1_ cgtxx14_
#define cqchh1_ cqchh14_
#else
#define cesc1 cesc14
#define ctx1 ctx14
#define cgtxx1 cgtxx14
#define cqchh1 cqchh14
#endif
#endif

/* f2cchar macro definition */
#if !defined(CRA) && !defined(ardent)
#define f2cchar(fptr) (fptr)
#endif /* default, CRA and ardent are exceptions */
#if defined(CRA) && !defined(ardent)
#include <fortran.h>
#define f2cchar(fptr) (_fcdtocp((_fcd)fptr))
#endif
#if defined(ardent)
struct FortranStr
{
  char *str;
  int   len;
  char  id;
};

typedef struct FortranStr Fortran_Str; /* Make the declarations in the */
                                       /* source more tractable. */
#define f2cchar(fptr) (fptr->str)
#endif /* ardent */

/* end ifdefc.h */

/******************************************************************************/
/******************************************************************************/
/*                                                                            */
/*      Standard CGI functions                                                */
/*                                                                            */
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*      Global variables                                                      */
/*                                                                            */
/******************************************************************************/
/* these are also used in mdcgi.c, where they are externs */
device_struct devices[MAX_DEVICES];     /* list of currently initialized */
anything *    in_params[MAX_IN_PARAMS]; /* params sent to driver */
short         num_devices = 0;          /* how many items in devices*/
anything *    sol_surf    = NULL;       /* current solicitation surface, */

/******************************************************************************/
/*                                                                            */
/*      ci - initialize CGI                                                   */
/*                                                                            */
/******************************************************************************/
void ci_(f_integer *pds)
{
  short func_id = CI_FN;
  short i;

  if (num_devices == 0) {
    /* find the default device and initialize it */
    cgi_def_ini();
  } /* end if no initialized devices */

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)pds;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end CI */

/******************************************************************************/
/*                                                                            */
/*      ct - terminate CGI                                                    */
/*                                                                            */
/******************************************************************************/
void ct_(void)
{
  short func_id = CT_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ct */

/******************************************************************************/
/*                                                                            */
/*      cxdfac - execute deferred actions                                     */
/*                                                                            */
/******************************************************************************/
void cxdfac_(void)
{
  short func_id = CXDFAC_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cxdfac */

/******************************************************************************/
/*                                                                            */
/*      cpds - prepare drawing surface                                        */
/*                                                                            */
/******************************************************************************/
void cpds_(f_integer *clear)
{
  short func_id = CPDS_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)clear;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cpds */

/******************************************************************************/
/*                                                                            */
/*      cendpg - end page                                                     */
/*                                                                            */
/******************************************************************************/
void cendpg_(void)
{
  short func_id = CENDPG_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cendpg */

/******************************************************************************/
/*                                                                            */
/*      cbc - set background color                                            */
/*                                                                            */
/******************************************************************************/
void cbc_(f_integer *red, f_integer *green, f_integer *blue)
{
  short func_id = CBC_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)red;
  in_params[2] = (anything *)green;
  in_params[3] = (anything *)blue;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cbc */

/******************************************************************************/
/*                                                                            */
/*      cvdcx - VDC Extent                                                    */
/*                                                                            */
/******************************************************************************/
void cvdcx_(f_real *x1, f_real *y1, f_real *x2, f_real *y2)
{
  short func_id = CVDCX_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)x1;
  in_params[2] = (anything *)y1;
  in_params[3] = (anything *)x2;
  in_params[4] = (anything *)y2;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cvdcx */

/******************************************************************************/
/*                                                                            */
/*      cv - Device Viewport                                                  */
/*                                                                            */
/******************************************************************************/
void cv_(f_real *x1, f_real *y1, f_real *x2, f_real *y2)
{
  short func_id = CV_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)x1;
  in_params[2] = (anything *)y1;
  in_params[3] = (anything *)x2;
  in_params[4] = (anything *)y2;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cv */

/******************************************************************************/
/*                                                                            */
/*      ccl - Clip Indicator                                                  */
/*                                                                            */
/******************************************************************************/
void ccl_(f_integer *clipi)
{
  short func_id = CCL_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)clipi;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ccl */

/******************************************************************************/
/*                                                                            */
/*      cdscl - Drawing Surface Clip Indicator                                */
/*                                                                            */
/******************************************************************************/
void cdscl_(f_integer *clipi)
{
  short func_id = CDSCL_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)clipi;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cdsl */

/******************************************************************************/
/*                                                                            */
/*      cdqerr - dequeue error reports                                        */
/*                                                                            */
/******************************************************************************/
void cdqerr_(f_integer *nreq, f_integer *vstat, f_integer *nrem, f_integer *nret, f_integer *errcl,
             f_integer *errnm, f_integer *funcid)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CDQERR_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)nreq;
  in_params[2] = (anything *)vstat;
  in_params[3] = (anything *)nrem;
  in_params[4] = (anything *)nret;
  in_params[5] = (anything *)errcl;
  in_params[6] = (anything *)errnm;
  in_params[7] = (anything *)funcid;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cdqerr */

/******************************************************************************/
/*                                                                            */
/*      cerhct - error handling control                                       */
/*                                                                            */
/******************************************************************************/
void cerhct_(f_integer *n, f_integer *erclas, f_integer *hflag)
{
  short func_id = CERHCT_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)n;
  in_params[2] = (anything *)erclas;
  in_params[3] = (anything *)hflag;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cerht */

/******************************************************************************/
/*                                                                            */
/*      ccixp - color index precision                                         */
/*                                                                            */
/******************************************************************************/
void ccixp_(f_integer *cip)
{
  short func_id = CCIXP_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)cip;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ccixp */

/******************************************************************************/
/*                                                                            */
/*      cesc - escape function                                                */
/*                                                                            */
/******************************************************************************/
void cesc1_(f_integer *funcid, f_integer *ldr, char *data, f_integer *drec_size)
{
  char *data1;
  short func_id = CESC_FN;
  short i;

  data1 = f2cchar(data); /* convert Fortran char ptr to C ptr */

  if (num_devices == 0) {
    /* find the default device and initialize it */
    cgi_def_ini();
  } /* end if no initialized devices */

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)funcid;
  in_params[2] = (anything *)ldr;
  in_params[3] = (anything *)data1;
  in_params[4] = (anything *)drec_size;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cesc */

/******************************************************************************/
/*                                                                            */
/*      cqid - inquire device identification                                  */
/*                                                                            */
/******************************************************************************/
void cqid_(f_integer *maxchr, f_integer *vstat, f_integer *dclass, char *devid)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQID_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)maxchr;
  in_params[2] = (anything *)vstat;
  in_params[3] = (anything *)dclass;
  in_params[4] = (anything *)devid;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqid */

/******************************************************************************/
/*                                                                            */
/*      cqd - inquire device description                                      */
/*                                                                            */
/******************************************************************************/
void cqd_(f_integer *vstat, f_integer *hscopy, f_integer *disp, f_integer *bcolor, f_integer *dynbc,
          f_integer *dynvdm, f_integer *dx1, f_integer *dy1, f_integer *dx2, f_integer *dy2,
          f_real *width, f_real *height, f_integer *pixloc)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQD_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)vstat;
  in_params[2]  = (anything *)hscopy;
  in_params[3]  = (anything *)disp;
  in_params[4]  = (anything *)bcolor;
  in_params[5]  = (anything *)dynbc;
  in_params[6]  = (anything *)dynvdm;
  in_params[7]  = (anything *)dx1;
  in_params[8]  = (anything *)dy1;
  in_params[9]  = (anything *)dx2;
  in_params[10] = (anything *)dy2;
  in_params[11] = (anything *)width;
  in_params[12] = (anything *)height;
  in_params[13] = (anything *)pixloc;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqd */

/******************************************************************************/
/*                                                                            */
/*      clf - Lookup Function Support                                         */
/*                                                                            */
/******************************************************************************/
void clf_(f_integer *n, f_integer *funccl, f_integer *funcid, f_integer *vstat, f_integer *supprt)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CLF_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)n;
  in_params[2] = (anything *)funccl;
  in_params[3] = (anything *)funcid;
  in_params[4] = (anything *)vstat;
  in_params[5] = (anything *)supprt;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end clf */

/******************************************************************************/
/*                                                                            */
/*      clpr - Lookup Profile Support                                         */
/*                                                                            */
/******************************************************************************/
void clpr_(f_integer *n, char *profid, f_integer *profid_size, f_integer *vstat, f_integer *supprt)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CLPR_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)n;
  in_params[2] = (anything *)profid;
  in_params[3] = (anything *)profid_size;
  in_params[4] = (anything *)vstat;
  in_params[5] = (anything *)supprt;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end clpr */

/******************************************************************************/
/*                                                                            */
/*      cqsp - inquire supported precisions                                   */
/*                                                                            */
/******************************************************************************/
void cqsp_(f_integer *vstat, f_integer *nvip, f_integer *vip, f_integer *nvrp, f_integer *vrfmt,
           f_integer *vrexp, f_integer *vrfrac, f_integer *nip, f_integer *ip, f_integer *nrp,
           f_integer *rfmt, f_integer *rexp, f_integer *rfrac, f_integer *nixp, f_integer *ixp,
           f_integer *ncp, f_integer *cp, f_integer *ncixp, f_integer *cixp)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQSP_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)vstat;
  in_params[2]  = (anything *)nvip;
  in_params[3]  = (anything *)vip;
  in_params[4]  = (anything *)nvrp;
  in_params[5]  = (anything *)vrfmt;
  in_params[6]  = (anything *)vrexp;
  in_params[7]  = (anything *)vrfrac;
  in_params[8]  = (anything *)nip;
  in_params[9]  = (anything *)ip;
  in_params[10] = (anything *)nrp;
  in_params[11] = (anything *)rfmt;
  in_params[12] = (anything *)rexp;
  in_params[13] = (anything *)rfrac;
  in_params[14] = (anything *)nixp;
  in_params[15] = (anything *)ixp;
  in_params[16] = (anything *)ncp;
  in_params[17] = (anything *)cp;
  in_params[18] = (anything *)ncixp;
  in_params[19] = (anything *)cixp;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqsp */

/******************************************************************************/
/*                                                                            */
/*      clesc - Lookup Escape Support                                         */
/*                                                                            */
/******************************************************************************/
void clesc_(f_integer *n, f_integer *escid, f_integer *vstat, f_integer *supprt)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CLESC_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)n;
  in_params[2] = (anything *)escid;
  in_params[3] = (anything *)vstat;
  in_params[4] = (anything *)supprt;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end clesc  */

/******************************************************************************/
/*                                                                            */
/*      cqp - inquire current precisions                                      */
/*                                                                            */
/******************************************************************************/
void cqp_(f_integer *vstat, f_integer *vip, f_integer *vrfmt, f_integer *vrexp, f_integer *vrfrac,
          f_integer *ip, f_integer *rfmt, f_integer *rexp, f_integer *rfrac, f_integer *ixp,
          f_integer *cp, f_integer *cixp)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQP_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)vstat;
  in_params[2]  = (anything *)vip;
  in_params[3]  = (anything *)vrfmt;
  in_params[4]  = (anything *)vrexp;
  in_params[5]  = (anything *)vrfrac;
  in_params[6]  = (anything *)ip;
  in_params[7]  = (anything *)rfmt;
  in_params[8]  = (anything *)rexp;
  in_params[9]  = (anything *)rfrac;
  in_params[10] = (anything *)ixp;
  in_params[11] = (anything *)cp;
  in_params[12] = (anything *)cixp;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqp */

/******************************************************************************/
/*                                                                            */
/*      cqcl - inquire clipping                               */
/*                                                                            */
/******************************************************************************/
void cqcl_(f_integer *vstat, f_integer *clip1, f_integer *clipr, f_integer *sclip1,
           f_integer *sclipr)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQCL_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)vstat;
  in_params[2] = (anything *)clip1;
  in_params[3] = (anything *)clipr;
  in_params[4] = (anything *)sclip1;
  in_params[5] = (anything *)sclipr;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqcl */

/******************************************************************************/
/*                                                                            */
/*      cpl - polyline                                                        */
/*                                                                            */
/******************************************************************************/
void cpl_(f_integer *np, f_real *px, f_real *py)
{
  short func_id = CPL_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)np;
  in_params[2] = (anything *)px;
  in_params[3] = (anything *)py;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cpl */

/******************************************************************************/
/*                                                                            */
/*      cdjpl - disjoint polyline                                             */
/*                                                                            */
/******************************************************************************/
void cdjpl_(f_integer *np, f_real *px, f_real *py)
{
  short func_id = CDJPL_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)np;
  in_params[2] = (anything *)px;
  in_params[3] = (anything *)py;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cdjpl */

/******************************************************************************/
/*                                                                            */
/*      cpm - polymarker                                                      */
/*                                                                            */
/******************************************************************************/
void cpm_(f_integer *np, f_real *px, f_real *py)
{
  short func_id = CPM_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)np;
  in_params[2] = (anything *)px;
  in_params[3] = (anything *)py;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cpm */

/******************************************************************************/
/*                                                                            */
/*      ctx - text                                                            */
/*                                                                            */
/******************************************************************************/
void ctx1_(f_real *x, f_real *y, f_integer *flag, char *text, f_integer *text_size)
{
  char *text1;
  short func_id = CTX_FN;
  short i;

  text1 = f2cchar(text); /* convert Fortran char ptr to C ptr */

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)x;
  in_params[2] = (anything *)y;
  in_params[3] = (anything *)flag;
  in_params[4] = (anything *)text1;
  in_params[5] = (anything *)text_size;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ctx */

/******************************************************************************/
/*                                                                            */
/*      cpg - polygon                                                          */
/*                                                                            */
/******************************************************************************/
void cpg_(f_integer *np, f_real *px, f_real *py)
{
  short func_id = CPG_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)np;
  in_params[2] = (anything *)px;
  in_params[3] = (anything *)py;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cpg */

/******************************************************************************/
/*                                                                            */
/*      cca - cell array                                                      */
/*                                                                            */
/******************************************************************************/
void cca_(f_real *x1, f_real *y1, f_real *x2, f_real *y2, f_real *x3, f_real *y3, f_integer *nx,
          f_integer *ny, f_integer *lcp, f_integer *cells)
{
  short func_id = CCA_FN;
  short i;

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)x1;
  in_params[2]  = (anything *)y1;
  in_params[3]  = (anything *)x2;
  in_params[4]  = (anything *)y2;
  in_params[5]  = (anything *)x3;
  in_params[6]  = (anything *)y3;
  in_params[7]  = (anything *)nx;
  in_params[8]  = (anything *)ny;
  in_params[9]  = (anything *)lcp;
  in_params[10] = (anything *)cells;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cca */

/******************************************************************************/
/*                                                                            */
/*      clnt - line type                                                      */
/*                                                                            */
/******************************************************************************/
void clnt_(f_integer *lntyp)
{
  short func_id = CLNT_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)lntyp;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end clnt */

/******************************************************************************/
/*                                                                            */
/*      clnw - line width                                                     */
/*                                                                            */
/******************************************************************************/
void clnw_(f_real *lnwid)
{
  short func_id = CLNW_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)lnwid;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end clnw */

/******************************************************************************/
/*                                                                            */
/*      clnc - line color                                                     */
/*                                                                            */
/******************************************************************************/
void clnc_(f_integer *lnclr)
{
  short func_id = CLNC_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)lnclr;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end clnc */

/******************************************************************************/
/*                                                                            */
/*      cmkt - marker type                                                    */
/*                                                                            */
/******************************************************************************/
void cmkt_(f_integer *mktyp)
{
  short func_id = CMKT_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)mktyp;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cmkt*/

/******************************************************************************/
/*                                                                            */
/*      cmkc - marker color                                                   */
/*                                                                            */
/******************************************************************************/
void cmkc_(f_integer *mkclr)
{
  short func_id = CMKC_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)mkclr;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cmkc*/

/******************************************************************************/
/*                                                                            */
/*      ctxp - text precision                                         */
/*                                                                            */
/******************************************************************************/
void ctxp_(f_integer *txp)
{
  short func_id = CTXP_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)txp;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ctxp */

/******************************************************************************/
/*                                                                            */
/*      ctxc - text color                                                     */
/*                                                                            */
/******************************************************************************/
void ctxc_(f_integer *txclr)
{
  short func_id = CTXC_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)txclr;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ctxc */

/******************************************************************************/
/*                                                                            */
/*      cchh - character height                                               */
/*                                                                            */
/******************************************************************************/
void cchh_(f_real *chhit)
{
  short func_id = CCHH_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)chhit;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cchh */

/******************************************************************************/
/*                                                                            */
/*      ccho - character orientation                             */
/*                                                                            */
/******************************************************************************/
void ccho_(f_integer *xup, f_integer *yup, f_integer *xbase, f_integer *ybase)
{
  short func_id = CCHO_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)xup;
  in_params[2] = (anything *)yup;
  in_params[3] = (anything *)xbase;
  in_params[4] = (anything *)ybase;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ccho */

/******************************************************************************/
/*                                                                            */
/*      cis - interior style                                                  */
/*                                                                            */
/******************************************************************************/
void cis_(f_integer *instyl)
{
  short func_id = CIS_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)instyl;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cis*/

/******************************************************************************/
/*                                                                            */
/*      cflc - fill colour specifier                                          */
/*                                                                            */
/******************************************************************************/
void cflc_(f_integer *fclr)
{
  short func_id = CFLC_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)fclr;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cflc*/

/******************************************************************************/
/*                                                                            */
/*      ccsm - color selection mode                                           */
/*                                                                            */
/******************************************************************************/
void ccsm_(f_integer *csmode)
{
  short func_id = CCSM_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)csmode;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end ccsm */

/******************************************************************************/
/*                                                                            */
/*      cct - color table                                                     */
/*                                                                            */
/******************************************************************************/
void cct_(f_integer *starti, f_integer *nclrs, f_integer *clrs)
{
  short func_id = CCT_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)starti;
  in_params[2] = (anything *)nclrs;
  in_params[3] = (anything *)clrs;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cct */

/******************************************************************************/
/*                                                                            */
/*      cgtxx - get text extent                                               */
/*                                                                            */
/******************************************************************************/
void cgtxx1_(f_real *x, f_real *y, char *string, f_integer *vstat, f_integer *vconc, f_real *xconc,
             f_real *yconc, f_real *x1, f_real *y1, f_real *x2, f_real *y2, f_real *x3, f_real *y3,
             f_real *x4, f_real *y4, f_integer *string_size)
{
  char *    string1;
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CGTXX_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  string1 = f2cchar(string); /* convert Fortran char ptr to C ptr */

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)x;
  in_params[2]  = (anything *)y;
  in_params[3]  = (anything *)string1;
  in_params[4]  = (anything *)vstat;
  in_params[5]  = (anything *)vconc;
  in_params[6]  = (anything *)xconc;
  in_params[7]  = (anything *)yconc;
  in_params[8]  = (anything *)x1;
  in_params[9]  = (anything *)y1;
  in_params[10] = (anything *)x2;
  in_params[11] = (anything *)y2;
  in_params[12] = (anything *)x3;
  in_params[13] = (anything *)y3;
  in_params[14] = (anything *)x4;
  in_params[15] = (anything *)y4;
  in_params[16] = (anything *)string_size;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cgtxx */

/******************************************************************************/
/*                                                                            */
/*      cqprl - inquire primitive support levels                              */
/*                                                                            */
/******************************************************************************/
void cqprl_(f_integer *vstat, f_integer *maxpl, f_integer *maxdpl, f_integer *maxpg,
            f_integer *maxpgs, f_integer *maxpm, f_integer *maxcf, f_integer *maxchr,
            f_integer *maxcel, f_integer *celfil, f_integer *celaln, f_integer *comptx,
            f_integer *clofig)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQPRL_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)vstat;
  in_params[2]  = (anything *)maxpl;
  in_params[3]  = (anything *)maxdpl;
  in_params[4]  = (anything *)maxpg;
  in_params[5]  = (anything *)maxpgs;
  in_params[6]  = (anything *)maxpm;
  in_params[7]  = (anything *)maxcf;
  in_params[8]  = (anything *)maxchr;
  in_params[9]  = (anything *)maxcel;
  in_params[10] = (anything *)celfil;
  in_params[11] = (anything *)celaln;
  in_params[12] = (anything *)comptx;
  in_params[13] = (anything *)clofig;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqprl */

/******************************************************************************/
/*                                                                            */
/*      cqln - inquire line capability                                        */
/*                                                                            */
/******************************************************************************/
void cqln_(f_integer *vstat, f_integer *npdefb, f_integer *nsetb, f_integer *maxbi,
           f_integer *dynmod, f_integer *nomwid, f_integer *minwid, f_integer *maxwid)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQLN_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)vstat;
  in_params[2] = (anything *)npdefb;
  in_params[3] = (anything *)nsetb;
  in_params[4] = (anything *)maxbi;
  in_params[5] = (anything *)dynmod;
  in_params[6] = (anything *)nomwid;
  in_params[7] = (anything *)minwid;
  in_params[8] = (anything *)maxwid;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqln*/

/******************************************************************************/
/*                                                                            */
/*      cqlnt - inquire list of available line types                          */
/*                                                                            */
/******************************************************************************/
void cqlnt_(f_integer *nreq, f_integer *first, f_integer *vstat, f_integer *ntotal,
            f_integer *nlist, f_integer *lntyp)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQLNT_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)nreq;
  in_params[2] = (anything *)first;
  in_params[3] = (anything *)vstat;
  in_params[4] = (anything *)ntotal;
  in_params[5] = (anything *)nlist;
  in_params[6] = (anything *)lntyp;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqlnt*/

/******************************************************************************/
/*                                                                            */
/*      cqchh - inquire list of available character heights                   */
/*                                                                            */
/******************************************************************************/
void cqchh1_(char *font, f_integer *txp, f_integer *nreq, f_integer *first, f_integer *vstat,
             f_integer *ntotal, f_integer *nlist, f_integer *chhit, f_integer *font_size)
{
  char *    font1;
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQCHH_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  font1 = f2cchar(font); /* convert Fortran char ptr to C ptr */

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)font1;
  in_params[2] = (anything *)txp;
  in_params[3] = (anything *)nreq;
  in_params[4] = (anything *)first;
  in_params[5] = (anything *)vstat;
  in_params[6] = (anything *)ntotal;
  in_params[7] = (anything *)nlist;
  in_params[8] = (anything *)chhit;
  in_params[9] = (anything *)font_size;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqchh */

/******************************************************************************/
/*                                                                            */
/*      cqfl - inquire fill capability                                       */
/*                                                                            */
/******************************************************************************/
void cqfl_(f_integer *vstat, f_integer *npdefb, f_integer *nsetb, f_integer *maxbi,
           f_integer *dynmod, f_integer *ninsty, f_integer *instyl, f_integer *npdefp,
           f_integer *nsetp, f_integer *maxpi, f_integer *pdiv, f_integer *maxpx, f_integer *maxpy,
           f_integer *ptrans)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQFL_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)vstat;
  in_params[2]  = (anything *)npdefb;
  in_params[3]  = (anything *)nsetb;
  in_params[4]  = (anything *)maxbi;
  in_params[5]  = (anything *)dynmod;
  in_params[6]  = (anything *)ninsty;
  in_params[7]  = (anything *)instyl;
  in_params[8]  = (anything *)npdefp;
  in_params[9]  = (anything *)nsetp;
  in_params[10] = (anything *)maxpi;
  in_params[11] = (anything *)pdiv;
  in_params[12] = (anything *)maxpx;
  in_params[13] = (anything *)maxpy;
  in_params[14] = (anything *)ptrans;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqfl */

/******************************************************************************/
/*                                                                            */
/*      cqc - inquire colour capabilities                                     */
/*                                                                            */
/******************************************************************************/
void cqc_(f_integer *vstat, f_integer *nsimul, f_integer *navail, f_integer *nint, f_integer *cmode,
          f_integer *dynmod, f_integer *overit, f_integer *monoc)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQC_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)vstat;
  in_params[2] = (anything *)nsimul;
  in_params[3] = (anything *)navail;
  in_params[4] = (anything *)nint;
  in_params[5] = (anything *)cmode;
  in_params[6] = (anything *)dynmod;
  in_params[7] = (anything *)overit;
  in_params[8] = (anything *)monoc;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqc */

/******************************************************************************/
/*                                                                            */
/*      cqlna - inquire line attributes                                       */
/*                                                                            */
/******************************************************************************/
void cqlna_(f_integer *vstat, f_integer *lnbi, f_integer *lntyp, f_integer *lwmode,
            f_integer *lnwid, f_integer *csmode, f_integer *lnclr, f_integer *lcmode)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQLNA_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)vstat;
  in_params[2] = (anything *)lnbi;
  in_params[3] = (anything *)lntyp;
  in_params[4] = (anything *)lwmode;
  in_params[5] = (anything *)lnwid;
  in_params[6] = (anything *)csmode;
  in_params[7] = (anything *)lnclr;
  in_params[8] = (anything *)lcmode;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqlna */

/******************************************************************************/
/*                                                                            */
/*      cqtxa - inquire text attributes                                       */
/*                                                                            */
/******************************************************************************/
void cqtxa_(f_integer *vstat, f_integer *txbi, f_integer *fonti, f_integer *fontp, f_real *chexp,
            f_real *chspac, f_integer *csmode, f_integer *txclr[3], f_real *chhit,
            f_real *orient[4], f_integer *txpath, f_integer *horal, f_real *contha,
            f_integer *veral, f_real *contva, f_integer *chsi, f_integer *achsi)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQTXA_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0]  = (anything *)&func_id;
  in_params[1]  = (anything *)vstat;
  in_params[2]  = (anything *)txbi;
  in_params[3]  = (anything *)fonti;
  in_params[4]  = (anything *)fontp;
  in_params[5]  = (anything *)chexp;
  in_params[6]  = (anything *)chspac;
  in_params[7]  = (anything *)csmode;
  in_params[8]  = (anything *)txclr;
  in_params[9]  = (anything *)chhit;
  in_params[10] = (anything *)orient;
  in_params[11] = (anything *)txpath;
  in_params[12] = (anything *)horal;
  in_params[13] = (anything *)contha;
  in_params[14] = (anything *)veral;
  in_params[15] = (anything *)contva;
  in_params[16] = (anything *)chsi;
  in_params[17] = (anything *)achsi;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqtxa*/

/******************************************************************************/
/*                                                                            */
/*      cqcte - inquire list of colour table entries                          */
/*                                                                            */
/******************************************************************************/
void cqcte_(f_integer *nreq, f_integer *first, f_integer *vstat, f_integer *ntotal,
            f_integer *nlist, f_integer *colors)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CQCTE_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)nreq;
  in_params[2] = (anything *)first;
  in_params[3] = (anything *)vstat;
  in_params[4] = (anything *)ntotal;
  in_params[5] = (anything *)nlist;
  in_params[6] = (anything *)colors;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqcte */

/******************************************************************************/
/*                                                                            */
/*      cili - initialize logical input device                                */
/*                                                                            */
/******************************************************************************/
void cili_(f_integer *iclass, f_integer *idev)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CILI_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)iclass;
  in_params[2] = (anything *)idev;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cili*/

/******************************************************************************/
/*                                                                            */
/*      crqlc - request locator                                               */
/*                                                                            */
/******************************************************************************/
void crqlc_(f_integer *idev, f_real *timeout, f_integer *vstat, f_integer *rstat, f_integer *mvalid,
            f_integer *triggr, f_real *x, f_real *y)
{
  short     dev;       /* which device to look at now */
  short     dev_found; /* which device was it found on */
  short     func_id = CRQLC_FN;
  short     surf;           /* which surface on device */
  short     surf_found = 0; /* which active_surface was found */
  anything *temp_surface[1];

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)idev;
  in_params[2] = (anything *)timeout;
  in_params[3] = (anything *)vstat;
  in_params[4] = (anything *)rstat;
  in_params[5] = (anything *)mvalid;
  in_params[6] = (anything *)triggr;
  in_params[7] = (anything *)x;
  in_params[8] = (anything *)y;

  /* search active devices for this surface */
  dev_found = -1;
  for (dev = 0; (dev < num_devices) && (dev_found == -1); ++dev) {
    for (surf = 0; surf < devices[dev].num_active_surfaces; ++surf) {
      if (sol_surf == devices[dev].statelist[surf]) {
        dev_found  = dev;
        surf_found = surf;
        break;
      } /* end if found on list */
    }   /* end for */
  }     /* end for all devices */

  if (dev_found != -1) {
    temp_surface[0] = devices[dev_found].statelist[surf_found];
    (*devices[dev_found].device_fn)(in_params, 1, temp_surface);
  } /* end if surface not found */

} /* end cqrlc */

/******************************************************************************/
/*                                                                            */
/*      cpxa - pixel array                                                    */
/*                                                                            */
/******************************************************************************/
void cpxa_(f_real *x, f_real *y, f_integer *nx, f_integer *ny, f_integer *xscal, f_integer *yscal,
           f_integer *xdir, f_integer *ydir, f_integer *pxclrs)
{
  short func_id = CPXA_FN;
  short i;

  /* set up parameter array */
  in_params[0] = (anything *)&func_id;
  in_params[1] = (anything *)x;
  in_params[2] = (anything *)y;
  in_params[3] = (anything *)nx;
  in_params[4] = (anything *)ny;
  in_params[5] = (anything *)xscal;
  in_params[6] = (anything *)yscal;
  in_params[7] = (anything *)xdir;
  in_params[8] = (anything *)ydir;
  in_params[9] = (anything *)pxclrs;

  for (i = 0; i < num_devices; ++i) {
    (*devices[i].device_fn)(in_params, devices[i].num_on_surfaces, devices[i].statelist);
  } /* end for all devices */
} /* end cpxa */
