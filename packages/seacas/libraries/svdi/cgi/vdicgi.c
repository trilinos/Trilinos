/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 SVDI  CGI driver
    Implemented in C, callable by Fortran
    Debbie Campbell   JAN 1989
    Modifications: (mods marked by mod###)
      vbinq vdwait       13-jul-90   kcole
      vdstcs             15-aug-90    kcole
 */

#include "cgidef.h"
#include <stdio.h>
#include <string.h>
#include <time.h>
#define NO_GET_DEVID_CHAR
#include "devid.h"
#undef NO_GET_DEVID_CHAR

static float get_devid_num(char *string)
{
  int i;

  for (i = 0; i < MAX_DEVID; i++) {
    if (*(device_values[i].devid_char) == string[0] &&
        *(device_values[i].devid_char + 1) == string[1] &&
        *(device_values[i].devid_char + 2) == string[2]) {
      return (device_values[i].devid_num);
    }
  }

  /*  return a zero if there is no character string match */
  return (0.);
}

/* ifdefv.h - ifdef file for svdi routines
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
 * This include file is used by VDICGI.C
 */

#if defined(ADDC_)
#endif
#if !defined(CRA) && !defined(ADDC_) && !defined(COUGAR)
#define vdinit_ vdinit
#define vdterm_ vdterm
#define vdfram_ vdfram
#define vdiqdc_ vdiqdc
#define vdnwpg_ vdnwpg
#define vdbell_ vdbell
#define vdwait_ vdwait
#define vdbufl_ vdbufl
#define vdstco_ vdstco
#define vdiqco_ vdiqco
#define vdescp_ vdescp
#define vdiqes_ vdiqes
#define vdiqnd_ vdiqnd
#define vdmova_ vdmova
#define vdlina_ vdlina
#define vdpnta_ vdpnta
#define vdtext_ vdtext
#define vdpoly_ vdpoly
#define vdiqcp_ vdiqcp
#define vdiqos_ vdiqos
#define vdstos_ vdstos
#define vdstfc_ vdstfc
#define vdstbc_ vdstbc
#define vdstin_ vdstin
#define vdstls_ vdstls
#define vdstlw_ vdstlw
#define vdstcs_ vdstcs
#define vdaabu_ vdaabu
#define vdaloc_ vdaloc
#define vdabgl_ vdabgl
#define vdakgl_ vdakgl
#define vdstla_ vdstla
#define vdloge_ vdloge
#define vdmoni_ vdmoni
#define vbpkg _vbpkg
#define vberrh_ vberrh
#define vbiqpk_ vbiqpk
#define vbiqdv_ vbiqdv
#define vbdev _vbdev
#define cdrofs_ cdrofs
#define cdrrfs_ cdrrfs
#define cdrwfs_ cdrwfs
#define cdrcfs_ cdrcfs
#define cdroab_ cdroab
#define bgpbuf_ bgpbuf
#endif
#if defined(CRA)
#define vdinit_ VDINIT
#define vdterm_ VDTERM
#define vdfram_ VDFRAM
#define vdiqdc_ VDIQDC
#define vdnwpg_ VDNWPG
#define vdbell_ VDBELL
#define vdwait_ VDWAIT
#define vdbufl_ VDBUFL
#define vdstco_ VDSTCO
#define vdiqco_ VDIQCO
#define vdescp_ VDESCP
#define vdiqes_ VDIQES
#define vdiqnd_ VDIQND
#define vdmova_ VDMOVA
#define vdlina_ VDLINA
#define vdpnta_ VDPNTA
#define vdtext_ VDTEXT
#define vdpoly_ VDPOLY
#define vdiqcp_ VDIQCP
#define vdiqos_ VDIQOS
#define vdstos_ VDSTOS
#define vdstfc_ VDSTFC
#define vdstbc_ VDSTBC
#define vdstin_ VDSTIN
#define vdstls_ VDSTLS
#define vdstlw_ VDSTLW
#define vdstcs_ VDSTCS
#define vdaabu_ VDAABU
#define vdaloc_ VDALOC
#define vdabgl_ VDABGL
#define vdakgl_ VDAKGL
#define vdstla_ VDSTLA
#define vdloge_ VDLOGE
#define vdmoni_ VDMONI
#define vbpkg_ VBPKG
#define vberrh_ VBERRH
#define vbiqpk_ VBIQPK
#define vbiqdv_ VBIQDV
#define vbdev_ VBDEV
#define cdrofs_ CDROFS
#define cdrrfs_ CDRRFS
#define cdrwfs_ CDRWFS
#define cdrcfs_ CDRCFS
#define cdroab_ CDROAB
#define bgpbuf_ BGPBUF
#endif
#ifdef Build64
#if !defined(ADDC_)
#define vdinit vdinit4
#define vdiqdc vdiqdc4
#define vdstco vdstco4
#define vdiqco vdiqco4
#define vdescp vdescp4
#define vdiqes vdiqes4
#define vdiqnd vdiqnd4
#define vdmova vdmova4
#define vdlina vdlina4
#define vdpnta vdpnta4
#define vdtext vdtext4
#define vdpoly vdpoly4
#define vdiqcp vdiqcp4
#define vdiqos vdiqos4
#define vdstfc vdstfc4
#define vdstbc vdstbc4
#define vdstin vdstin4
#define vdstls vdstls4
#define vdstlw vdstlw4
#define vdstcs vdstcs4
#define vdaabu vdaabu4
#define vdaloc vdaloc4
#define vdabgl vdabgl4
#define vdakgl vdakgl4
#define vdstla vdstla4
#define vdstos vdstos4
#define vdfram vdfram4
#else
#define vdinit_ vdinit4_
#define vdiqdc_ vdiqdc4_
#define vdstco_ vdstco4_
#define vdiqco_ vdiqco4_
#define vdescp_ vdescp4_
#define vdiqes_ vdiqes4_
#define vdiqnd_ vdiqnd4_
#define vdmova_ vdmova4_
#define vdlina_ vdlina4_
#define vdpnta_ vdpnta4_
#define vdtext_ vdtext4_
#define vdpoly_ vdpoly4_
#define vdiqcp_ vdiqcp4_
#define vdiqos_ vdiqos4_
#define vdstfc_ vdstfc4_
#define vdstbc_ vdstbc4_
#define vdstin_ vdstin4_
#define vdstls_ vdstls4_
#define vdstlw_ vdstlw4_
#define vdstcs_ vdstcs4_
#define vdaabu_ vdaabu4_
#define vdaloc_ vdaloc4_
#define vdabgl_ vdabgl4_
#define vdakgl_ vdakgl4_
#define vdstla_ vdstla4_
#define vdstos_ vdstos4_
#define vdfram_ vdfram4_
#endif
#endif

/* end ifdefv.h */
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
#if !defined(ADDC_)
#define cesc cesc4
#define ctx ctx4
#define cgtxx cgtxx4
#define cqchh cqchh4
#define cesc2 cesc24
#define ctx2 ctx24
#define cgtxx2 cgtxx24
#define cqchh2 cqchh24
#else
#define cesc_ cesc4_
#define ctx_ ctx4_
#define cgtxx_ cgtxx4_
#define cqchh_ cqchh4_
#define cesc2_ cesc24_
#define ctx2_ ctx24_
#define cgtxx2_ cgtxx24_
#define cqchh2_ cqchh24_
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

#include "sdcgi.h"

/* current position */
static float xcp = 0.;
static float ycp = 0.;

/* internal buffer for polylines, polygons */
#define VBUF_SIZE 2000
static float vlist_x[VBUF_SIZE];
static float vlist_y[VBUF_SIZE];
static int   nvert = 0;

/* transformation parameters*/
static float scale;

#define MAX_VECTOR 7
static float vector[MAX_VECTOR] = {7., 0., 1., 0., .01, 0., 0.};

#define MAX_DEV_CAP 33
static float dev_cap[MAX_DEV_CAP] = {0., 1., 1.,    255.,   0.,     15., 2., 1., 1., 1., 1.,
                                     1., 1., 0.,    32767., 32767., 0.,  0., 1., 1., 1., 0.,
                                     0., 3., 2000., 1.,     0.,     0.,  0., 1., 1., 1., 1.};
static float linewidth_nominal;

/*  Maximum NDC Values */
static float ndc_xmax = 1.;
static float ndc_ymax = 1.;
static int   alpha_mode;

static float color_scale;

/* Logical Color Table */
static int init_colors[24] = {0, 0, 0,   255, 0, 0,   0, 255, 0,   255, 255, 0,
                              0, 0, 255, 255, 0, 255, 0, 255, 255, 255, 255, 255};

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

/* macro to convert ascii(integer) to char (note: machine dependent) */
#define a_to_c(ain) ((char)ain) /* for ascii machine */

/* macros which map ndc into CGI coords. */
#define map_x(xin) ((float)(scale * (xin)))
#define map_y(yin) ((float)(scale * (yin)))

/* macros which map CGI coords. into ndc */
#define ndc_map_x(xin) ((float)((xin) / scale))
#define ndc_map_y(yin) ((float)((yin) / scale))

#define min(p1, p2) ((p1) < (p2) ? p1 : p2)
#define max(p1, p2) ((p1) > (p2) ? p1 : p2)

void vdicgi_errh(char errmsg[])
{
  int   temp1, temp2;
  float temp3;

  temp1 = XEAGMD;
  temp2 = 1;
  temp3 = 0.; /*  alpha  */
  cesc2_(&temp1, &temp2, &temp3);
  alpha_mode = TRUE;
  fprintf(stderr, " %s\n", errmsg);
}

void vbinq(void)
/* --  THIS IS WHERE CHANGES TO DEV_CAP OCCUR -- */
{
  int   vstat, hscopy, disp, bcolor, dynbc, dynvdm, dx1, dy1, dx2, dy2, pixloc;
  int   maxpl, maxdpl, maxpg, maxpgs, maxpm, maxcf, maxchr, maxcel;
  int   celfil, celaln, comptx, clofig, dclass;
  int   npdefb, nsetb, maxbi, dynmod, nomwid, minwid, maxwid;
  int   nreq, first, ntotal, lntyp[6], nlist;
  int   nsimul, navail, nint, cmode, overit, monoc, txp, chhit, i;
  float width, height;
  char  devid[4];

  /*  INQUIRE DEVICE IDENTIFICATION */
  maxchr = 3;
  cqid_(&maxchr, &vstat, &dclass, devid);
  if (vstat == CVAL) {
    dev_cap[22] = get_devid_num(devid);
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqid ");
  }

  /*  INQUIRE DEVICE DESCRIPTION */
  cqd_(&vstat, &hscopy, &disp, &bcolor, &dynbc, &dynvdm, &dx1, &dy1, &dx2, &dy2, &width, &height,
       &pixloc);
  if (vstat == CVAL) {
    dev_cap[0] = (float)hscopy; /* Erasibility */
    dev_cap[1] = (float)disp;   /* Scan type (vector, raster) */
    /*     the next two may be wrong, want view surface instead of device surface
           : also note, the next 2 items may change if surface is window--later */
    if (!((dx1 == 0) && (dx2 == 0))) { /* X dimension view surface */
      dev_cap[14] = (float)(dx2 - dx1);
    }
    if (!((dy1 == 0) && (dy2 == 0))) { /* Y dimension view surface */
      dev_cap[15] = (float)(dy2 - dy1);
    }
    dev_cap[16] = width;  /* X dimension physical units */
    dev_cap[17] = height; /* Y dimension physical units */
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqd ");
  }

  /*    input inquiries to zeros - for non-interactive devices */
  /*    mod### set these based on dclass not hscopy            */
  if (dclass == 0) /* dclass = dev(13) value */
  {
    for (i = 7; i < 13; i++) {
      dev_cap[i] = 0.;
    }
  }

  /*  INQUIRE PRIMITIVE SUPPORT LEVELS */
  cqprl_(&vstat, &maxpl, &maxdpl, &maxpg, &maxpgs, &maxpm, &maxcf, &maxchr, &maxcel, &celfil,
         &celaln, &comptx, &clofig);
  if (vstat == CVAL) {
    if (maxpg == -1) {
      maxpg = VBUF_SIZE;
    }
    dev_cap[24] = min(maxpg, VBUF_SIZE); /* Maximum polygon points */
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqprl ");
  }

  /*  INQUIRE LINE CAPABILITY */
  cqln_(&vstat, &npdefb, &nsetb, &maxbi, &dynmod, &nomwid, &minwid, &maxwid);
  if (vstat == CVAL) {
    dev_cap[18] = (float)minwid; /* Minimum line width */

    /*  Due to the fact that cgi device coordinates have to be integer, the
        linewidth numbers may come back zero, enforce a minimum for nominal
        since we use it later in vdstlw for calculations */
    linewidth_nominal = (nomwid == 0) /* Nominal line width */
                            ? .001
                            : (float)nomwid / dev_cap[14];
    dev_cap[30] = (float)maxwid; /* Maximum line width */

    /* Could inquire marker capability, but this will be close enough*/
    dev_cap[19] = dev_cap[18]; /* Minimum pointsize  */
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqln ");
  }

  /*  INQUIRE LIST OF AVAILABLE LINE TYPES */
  nreq  = 6;
  first = 1;
  cqlnt_(&nreq, &first, &vstat, &ntotal, &nlist, lntyp);
  if (vstat == CVAL) {
    /* map CGI to VDI linestyles */
    /* cgi linestyles: 1 - solid 2 - dash 3 - dot,
     *                 4 - dashdot 5 - dash dot dot
     * vdi linestyles: 0 - solid 1 - dotted 2 - dot dash
     *                 3 - short dash 4 - long dash 5 - medium dash
     */
    dev_cap[5] = 0.;
    for (i = 0; i < nlist; i++) {
      switch (lntyp[i]) {
      case 1: break;
      case 2: dev_cap[5] += 16.; break;
      case 3: dev_cap[5] += 1.; break;
      case 4: dev_cap[5] += 2.; break;
      case 5: dev_cap[5] += 4.; break;
      default: break;
      } /* end switch */
    }
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqlnt ");
  }

  /*  INQUIRE COLOR CAPABILITIES */
  cqc_(&vstat, &nsimul, &navail, &nint, &cmode, &dynmod, &overit, &monoc);
  if (vstat == CVAL) {
    dev_cap[3]  = (float)nsimul - 2; /* Simultaneous colors    */
    dev_cap[26] = (float)navail - 2; /* Available colors       */
    dev_cap[2]  = (float)nint;       /* Available intensities   */
    dev_cap[31] = 1.;                /* Color or monochrome    */
    if (monoc == CYES) {
      dev_cap[31] = 0.;
    }
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqc ");
  }

  /*  INQUIRE LIST OF AVAILABLE CHARACTER HEIGHTS */
  txp   = CSTRNG;
  nreq  = 1;
  first = 1;
  cqchh2_(&txp, &nreq, &first, &vstat, &ntotal, &nlist, &chhit);
  if (vstat == CVAL) {
    dev_cap[6] = (nlist == 0) ? 0. : (float)ntotal;
  }
  else {
    vdicgi_errh(" SVDI Shell (VBINQ) invalid inquire for cqchh ");
  }
}

void vdinit_(float *aspect, int *justif)
{
  float asp;
  int   just, temp, temp2, vstat, vconc;
  float xconc, yconc, x1, y1, x2, y2, x3, y3, x4, y4, temp_xcp, temp_ycp;
  float scaled_ndc_xmax, scaled_ndc_ymax;
  float rtemp = 0.0;
  xconc       = 0.0;
  yconc       = 0.0;
  x1          = 0.0;
  y1          = 0.0;
  x2          = 0.0;
  y2          = 0.0;
  x3          = 0.0;
  y3          = 0.0;
  x4          = 0.0;
  y4          = 0.0;

  asp  = *aspect;
  just = *justif;

  if (asp < 0.) {
    vdicgi_errh(" SVDI Shell (VDINIT) Error Number 721 Severity 5: ");
    asp = 0.;
  }

  if (just < 0 || just > 9) {
    vdicgi_errh(" SVDI Shell (VDINIT) Error Number 720 Severity 5: ");
    just = 0;
  }

  /*  Initialize CGI         */
  temp       = CACT;
  alpha_mode = FALSE;
  ci_(&temp);

  /*  Inquire everything you always wanted to know about */
  vbinq();

  /*  Turn off clip indicators */
  temp = CDCOFF;
  cdscl_(&temp);
  temp = COFF;
  ccl_(&temp);

  /*  Set up proper scaling to take advantage of whole device (not just square) */
  if (asp == 0.) {
    asp = dev_cap[14] / dev_cap[15];
  }
  if (asp > 1.) {
    ndc_xmax = 1.;
    ndc_ymax = 1. / asp;
  }
  else {
    ndc_xmax = asp;
    ndc_ymax = 1.;
  }
  scale           = 32767.;
  scaled_ndc_xmax = map_x(ndc_xmax);
  scaled_ndc_ymax = map_y(ndc_ymax);
  rtemp           = 0.0;
  cvdcx_(&rtemp, &rtemp, &scaled_ndc_xmax, &scaled_ndc_ymax);

  /*  Set color mode to index, and set color index precision to 8 bits  */
  temp = CINDEX;
  ccsm_(&temp);
  temp = 8;
  ccixp_(&temp);
  color_scale = 255.;

  /*  set up the standard 8 colors in indices 2 - 9 (0 reserved for background,
      1 reserved for default foreground) */
  temp  = 2;
  temp2 = 8;
  cct_(&temp, &temp2, init_colors);

  /*  Set default marker type to dot  */
  temp = 1;
  cmkt_(&temp);

  /*  Set default interior style to solid */
  temp = CSOLID;
  cis_(&temp);

  /*  Inquire what the default character size is - use cgtxx instead of
      cqtxa because need both x and y size (may have to adjust for inter
      character/line spacing later   */
  temp_xcp = map_x(xcp);
  temp_ycp = map_y(ycp);
  cgtxx2_(&temp_xcp, &temp_ycp, &vstat, &vconc, &xconc, &yconc, &x1, &y1, &x2, &y2, &x3, &y3, &x4,
          &y4);
  if (vstat == CVAL) {
    vector[5] = ndc_map_x(x2 - x1);
    vector[6] = ndc_map_y(y4 - y1);
  }
  else {
    vdicgi_errh(" SVDI Shell (VDINIT) inquire error from cgtxx ");
  }

  /* Initialize locator device */
  temp  = CLOCAT;
  temp2 = 1;
  cili_(&temp, &temp2);
}

void vflush(void)
/*  flush polyline buffer */
{
  if (alpha_mode) {
    int   temp1 = XEAGMD;
    int   temp2 = 1;
    float temp3 = 1.; /*  graphics  */
    cesc2_(&temp1, &temp2, &temp3);
    alpha_mode = FALSE;
  }

  if (nvert > 0) {
    cpl_(&nvert, vlist_x, vlist_y);
    nvert = 0;
  }
}

void vdterm_()
{
  vflush();
  ct_();
}

void vdiqdc_(int *index, float *value)
{
  if (*index < 1 || *index > MAX_DEV_CAP) {
    vdicgi_errh(" SVDI Shell (VDIQDC) Error Number 726, Severity Code 5 ");
    return;
  }
  {
    *value = dev_cap[*index - 1];
  }
}

void vdnwpg_()
{
  int temp, nreq, first, vstat, ntotal, colors[3], nlist;

  vflush(); /*  flush polyline buffer  */

  /* execute deferred actions */
  cxdfac_();

  /*  inquire the rgb values of color_index. Set index 0 to those values   */

  nreq = 1;

  /*  the 3 below is because you asking for the 3rd index -- which happens
      to be logical index 0 */

  first = vector[1] + 3;
  cqcte_(&nreq, &first, &vstat, &ntotal, &nlist, colors);
  if (vstat == CVAL) {
    cbc_(&colors[0], &colors[1], &colors[2]);
  }
  else {
    vdicgi_errh(" SVDI Shell (VDNWPG) invalid inquire for cqcte ");
  }

  /*  prepare drawing surface - do background color   */
  temp = CCONDC;
  cpds_(&temp);
}

void vdbell_() {}

void vdwait_()
{
  int   temp, vstat, rstat, mvalid, trigger;
  float xpos, ypos, timeout;

  /*  make picture current - flush output buffers */
  vflush(); /*  flush polyline buffer  */
  cxdfac_();

  /*  do read waiting for viewing for interactive devices only */
  /*  mod### changed from 0(erasibility) to 12(input)          */
  if (dev_cap[12] != 0.) {
    temp    = 1;
    timeout = 1.;
    crqlc_(&temp, &timeout, &vstat, &rstat, &mvalid, &trigger, &xpos, &ypos);
    if (vstat != CVAL) {
      vdicgi_errh(" SVDI Shell (VDWAIT) invalid request for crqlc ");
    }
  }
}

void vdbufl_()
{
  int   temp1, temp2;
  float temp3;

  vflush(); /*  flush polyline buffer  */
  cxdfac_();

  /*   set the device into alpha mode - 0. */
  temp1 = XEAGMD;
  temp2 = 1;
  temp3 = 0.; /*  alpha  */
  cesc2_(&temp1, &temp2, &temp3);
  alpha_mode = TRUE;
}

void vdstco_(int *num, int index_array[], float color_array[][3], int *color_mod)
{
  int i              = 0;
  int start_index    = 0;
  int hld_colors_ptr = 0;
  int first          = 0;
  int count          = 0;
  int hld_colors[768];

  if (*num < 1 || *num > dev_cap[3]) {
    vdicgi_errh(" SVDI Shell (VDSTCO) Error Number 723, Severity Code 5 ");
    return;
  }

  if (*color_mod != 0 && *color_mod != 1) {
    vdicgi_errh(" SVDI Shell (VDSTCO) Error Number 725, Severity Code 5 ");
    return;
  }

  vflush(); /*  flush polyline buffer  */
  first = TRUE;

  for (i = 0; i < *num; i++) {
    if (index_array[i] < 0 || index_array[i] > 255) {
      vdicgi_errh(" SVDI Shell (VDSTCO) Error Number 724, Severity Code 5 ");

      /*  Set RGB color                                               */
    }
    else if (*color_mod == 0) {
      if (color_array[i][0] < 0. || color_array[i][0] > 1. || color_array[i][1] < 0. ||
          color_array[i][1] > 1. || color_array[i][2] < 0. || color_array[i][2] > 1.) {
        vdicgi_errh(" SVDI Shell (VDSTCO) Error Number 727, Severity Code 5 ");
      }
      else
      /*  Check to see if indices are consecutive (buffer up if they are)  */
      {

        if (first) {

          /*  set color table from index requested + 2 (0 reserved for background) */

          start_index    = index_array[i] + 2;
          hld_colors_ptr = 0;
          first          = FALSE;
        }

        hld_colors[hld_colors_ptr++] = (int)(color_array[i][0] * color_scale);
        hld_colors[hld_colors_ptr++] = (int)(color_array[i][1] * color_scale);
        hld_colors[hld_colors_ptr++] = (int)(color_array[i][2] * color_scale);

        if (i < *num - 1) {
          if (index_array[i] != (index_array[i + 1] - 1)) {
            count = hld_colors_ptr / 3;
            cct_(&start_index, &count, hld_colors);
            first = TRUE;
          }
        }
      }
    }

    /*  Set HLS color (HLS not supported)                           */
    else {
      vdicgi_errh(" SVDI Shell (VDSTCO) HLS option being phased out - ignored ");
      vdicgi_errh(" Contact Computer Graphics Group - Div. 2644 ");
    }
  }

  count = hld_colors_ptr / 3;
  if (count > 0) {
    cct_(&start_index, &count, hld_colors);
  }
}

void vdiqco_(int *num, int index_array[], float color_array[][3], int *color_mod)
{
  int   i, ntotal, colors[3], first, nreq, vstat, nlist;
  int   temp1, temp2;
  float temp3;

  if (*num < 1 || *num > dev_cap[3]) {
    vdicgi_errh(" SVDI Shell (VDIQCO) Error Number 723, Severity Code 5 ");
    return;
  }

  if (*color_mod != 0 && *color_mod != 1) {
    vdicgi_errh(" SVDI Shell (VDIQCO) Error Number 725, Severity Code 5 ");
    return;
  }

  if (alpha_mode) {
    temp1 = XEAGMD;
    temp2 = 1;
    temp3 = 1.; /* "graphics" */
    cesc2_(&temp1, &temp2, &temp3);
    alpha_mode = FALSE;
  }

  for (i = 0; i < *num; i++) {
    if (index_array[i] < 0 || index_array[i] > 255) {
      vdicgi_errh(" SVDI Shell (VDIQCO) Error Number 724, Severity Code 5 ");
    }
    else {
      if (*color_mod == 1) {
        vdicgi_errh(" HLS option being phased out - returning RGB values ");
        vdicgi_errh(" Contact Computer Graphics Group - Div. 2644 ");
      }

      /*  find out from cgi what the color table is */
      nreq  = 1;
      first = index_array[i] + 3;
      cqcte_(&nreq, &first, &vstat, &ntotal, &nlist, colors);
      if (vstat == CVAL) {
        color_array[i][0] = (float)colors[0] / color_scale;
        color_array[i][1] = (float)colors[1] / color_scale;
        color_array[i][2] = (float)colors[2] / color_scale;
      }
      else {
        vdicgi_errh(" SVDI Shell (VDIQCO) invalid inquire for cqcte ");
        color_array[i][0] = -1.;
      }
    }
  }
}

void vdescp_(int *escape_code, int *n, float args)
{
  if (*n < 0) {
    vdicgi_errh(" SVDI Shell (VDESCP) Error Number 802, Severity Code 5 ");
    return;
  }
}

void vdiqes_(int *escape_code, int *support) { *support = 0; }

void vdiqnd_(float *x_ndc, float *y_ndc)
{
  *x_ndc = ndc_xmax;
  *y_ndc = ndc_ymax;
}

void vdmova_(float *x, float *y)
{
  vflush();
  vlist_x[nvert] = map_x(*x);
  vlist_y[nvert] = map_y(*y);
  nvert++;

  xcp = *x;
  ycp = *y;
}

void vdlina_(float *x, float *y)
{

  /*  if a move hasn't already been done, doit it now  */
  if (nvert == 0) {
    vlist_x[nvert] = map_x(xcp);
    vlist_y[nvert] = map_y(ycp);
    nvert++;
  }

  /*  append to polyline buffer (if full, flush it first)  */
  if (nvert >= VBUF_SIZE) {
    vflush();
    nvert = 0; /* vflush sets this, but static analysis still warns. */
  }
  vlist_x[nvert] = map_x(*x);
  vlist_y[nvert] = map_y(*y);
  nvert++;

  xcp = *x;
  ycp = *y;
}

void vdpnta_(float *x, float *y)
{
  float new_x, new_y;
  int   temp;

  new_x = map_x(*x);
  new_y = map_y(*y);
  temp  = 1;

  vflush();
  vlist_x[nvert] = new_x;
  vlist_y[nvert] = new_y;
  nvert++;

  cpm_(&temp, &new_x, &new_y);

  xcp = *x;
  ycp = *y;
}

void vdtext_(int *length, int char_array[])
{
  int   i, lenout, len, temp1, temp2;
  float dx, dy, temp_xcp, temp_ycp, temp3;
  int   strout[150];

  len = *length;

  if (len < 1) {
    vdicgi_errh(" SVDI Shell (VDTEXT) Error Number 212, Severity Code 5 ");
    return;
  }

  if (len > 136) {
    vdicgi_errh(" SVDI Shell (VDTEXT) Error Number 213, Severity Code 5 ");
    len = 136;
  }

  if (alpha_mode) {
    temp1 = XEAGMD;
    temp2 = 1;
    temp3 = 1.; /*  graphics  */
    cesc2_(&temp1, &temp2, &temp3);
    alpha_mode = FALSE;
  }

  lenout = 0; /*count characters in string output buffer "strout" */

  for (i = 0; i < len; i++) {
    if (char_array[i] < 32 || char_array[i] > 126) {
      switch (char_array[i]) {
      case 8:
        dx = -vector[6];
        dy = 0.;
        break;
      case 10:
        dx = 0.;
        dy = -vector[5];
        break;
      case 13:
        dx = -xcp;
        dy = 0.;
        break;
      default:
        dx = 0.;
        dy = 0.;
        vdicgi_errh(" SVDI Shell (VDTEXT) Error Number 208, Severity Code 5 ");
        break;
      }
      /*  Stuff to send, finish the string   */

      if (lenout != 0) {
        temp_xcp = map_x(xcp);
        temp_ycp = map_y(ycp);
        ctx2_(&temp_xcp, &temp_ycp, strout, &lenout);
        xcp += lenout * vector[6];
        lenout = 0;
      }
      xcp += dx;
      ycp += dy;
      vdmova_(&xcp, &ycp);
    }

    else {

      strout[lenout++] = char_array[i];
    }
  }

  /*  All done, get rid of them         */
  if (lenout != 0) {
    temp_xcp = map_x(xcp);
    temp_ycp = map_y(ycp);
    ctx2_(&temp_xcp, &temp_ycp, strout, &lenout);
    xcp += lenout * vector[6];
    lenout = 0;
  }
}

void vdpoly_(float x_array[], float y_array[], int *npts)
{
  int i, iend;

  vflush();

  if (*npts > (int)dev_cap[24]) {
    vdicgi_errh(" SVDI Shell (VDPOLY) exceeded max points for - device limit ");
    iend = (int)dev_cap[24];
  }

  else if (*npts > VBUF_SIZE) {
    vdicgi_errh(" SVDI Shell (VDPOLY) exceeded max points - internal buffer limit ");
    iend = VBUF_SIZE;
  }
  else {
    iend = *npts;
  }

  for (i = 0; i < iend; i++) {
    vlist_x[nvert] = map_x(x_array[i]);
    vlist_y[nvert] = map_y(y_array[i]);
    nvert++;
  }

  cpg_(&nvert, vlist_x, vlist_y);
  nvert = 0;

  xcp = x_array[0];
  ycp = y_array[0];
}

void vdiqcp_(float *x, float *y)
{
  *x = xcp;
  *y = ycp;
}

void vdiqos_(float attr_array[])
{
  int i;

  for (i = 0; i < MAX_VECTOR; i++) {
    attr_array[i] = vector[i];
  }
}

void vdstfc_(int *color_index)
{
  int new_color;

  if (*color_index < 0 || *color_index > 255) {
    vdicgi_errh(" SVDI Shell (VDSTFC) Error Number 724 Severity 5: ");
    return;
  }

  vector[0] = (float)*color_index;

  /*   map foreground color into indices 2-n since 0 is reserved for background,
       and 1 is default foreground */
  new_color = *color_index + 2;

  vflush();

  clnc_(&new_color);
  cmkc_(&new_color);
  ctxc_(&new_color);
  cflc_(&new_color);
}

void vdstbc_(int *color_index)
{
  if (*color_index < 0 || *color_index > 255) {
    vdicgi_errh(" SVDI Shell (VDSTBC) Error Number 724 Severity 5: ");
    return;
  }
  vector[1] = (float)*color_index;
}

void vdstin_(float *intensity)
{
  if (*intensity < 0. || *intensity > 1.) {
    vdicgi_errh(" SVDI Shell (VDSTIN) Error Number 401 Severity 5: ");
    return;
  }
}

void vdstls_(int *line_style)
{
  int temp;

  if (*line_style < 0 || *line_style > 5) {
    vdicgi_errh(" SVDI Shell (VDSTLS) Error Number 401 Severity 5: ");
    return;
  }

  vector[3] = (float)*line_style;
  vflush();
  switch (*line_style) {
  /*  cgi linestyles are: 1 - solid, 2 - dash, 3 - dot, 4 - dashdot 5 - dash dot dot
      vdi linestyles are: 0 - solid, 1 - dotted, 2 - dot dash
      3 - short dash, 4 - long dash, 5 - medium dash    */
  case 0:
    temp = 1;
    clnt_(&temp);
    break;
  case 1:
    temp = 3;
    clnt_(&temp);
    break;
  case 2:
    temp = 4;
    clnt_(&temp);
    break;
  case 3:
    temp = 2;
    clnt_(&temp);
    break;
  case 4:
    temp = 5;
    clnt_(&temp);
    break;
  case 5:
    temp = 2;
    clnt_(&temp);
    break;
  default: vector[3] = 0.; break;
  }
}

void vdstlw_(float *line_width)
{
  float new_size;

  if (*line_width < 0. || *line_width > 1.) {
    vdicgi_errh(" SVDI Shell (VDSTLW) Error Number 401 Severity 5: ");
    return;
  }

  vector[4] = *line_width;

  /* cgi default state is scaled linewidth - not expressed in VDC */
  new_size = *line_width / linewidth_nominal;
  vflush();
  clnw_(&new_size);
}

void vdstcs_(float *y_size)
{
  float new_size, x, y, xconc, yconc, x1, y1, x2, y2, x3, y3, x4, y4;
  float temp3;
  int   vstat, vconc, temp1, temp2;

  xconc = 0.0;
  yconc = 0.0;
  x1    = 0.0;
  y1    = 0.0;
  x2    = 0.0;
  y2    = 0.0;
  x3    = 0.0;
  y3    = 0.0;
  x4    = 0.0;
  y4    = 0.0;

  new_size = map_x(*y_size);

  if (alpha_mode) {
    temp1 = XEAGMD;
    temp2 = 1;
    temp3 = 1.; /*  graphics  */
    cesc2_(&temp1, &temp2, &temp3);
    alpha_mode = FALSE;
  }

  cchh_(&new_size);

  /*  inquire what the realized size is and stuff it here */
  x = map_x(xcp);
  y = map_y(ycp);
  cgtxx2_(&x, &y, &vstat, &vconc, &xconc, &yconc, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4);
  if (vstat == CVAL) {
    /* mod## y value being stored incorrectyl; swapped x & y */
    vector[5] = ndc_map_y(y4 - y1);
    vector[6] = ndc_map_x(x2 - x1);
  }
  else {
    vdicgi_errh(" SVDI Shell (VDSTCS) inquire error from cgtxx ");
  }
}

void vdaabu_(int *button)
{
  int   temp, vstat, rstat, mvalid, trigger;
  float xpos, ypos, timeout;

  if (dev_cap[0] != 0.) {
    vflush(); /*  flush polyline buffer  */
    temp    = 1;
    timeout = 1.;
    crqlc_(&temp, &timeout, &vstat, &rstat, &mvalid, &trigger, &xpos, &ypos);
    if (vstat == CVAL) {
      *button = trigger;
    }
    else {
      vdicgi_errh(" SVDI Shell (VDAABU) invalid request for crqlc ");
    }
  }
}

void vdaloc_(float *x, float *y)
{
  int   temp, vstat, rstat, mvalid, trigger;
  float xpos, ypos, timeout;

  if (dev_cap[0] != 0.) {
    vflush(); /*  flush polyline buffer  */
    temp    = 1;
    timeout = 1.;
    crqlc_(&temp, &timeout, &vstat, &rstat, &mvalid, &trigger, &xpos, &ypos);
    if (vstat == CVAL) {
      *x = ndc_map_x(xpos);
      *y = ndc_map_y(ypos);
    }
    else {
      vdicgi_errh(" SVDI Shell (VDALOC) invalid request for crqlc ");
    }
  }
}

void vdabgl_(int *button, float *x, float *y)
{
  int   temp, vstat, rstat, mvalid, trigger;
  float xpos, ypos, timeout;

  if (dev_cap[0] != 0.) {
    vflush(); /*  flush polyline buffer  */
    temp    = 1;
    timeout = 1.;
    crqlc_(&temp, &timeout, &vstat, &rstat, &mvalid, &trigger, &xpos, &ypos);
    if (vstat == CVAL) {
      *x      = ndc_map_x(xpos);
      *y      = ndc_map_y(ypos);
      *button = trigger;
    }
    else {
      vdicgi_errh(" SVDI Shell (VDABGL) invalid request for crqlc ");
    }
  }
}

void vdakgl_(int *charac, float *x, float *y)
{
  int   temp, vstat, rstat, mvalid, trigger;
  float xpos, ypos, timeout;

  if (dev_cap[0] != 0.) {
    vflush(); /*  flush polyline buffer  */
    temp    = 1;
    timeout = 1.;
    crqlc_(&temp, &timeout, &vstat, &rstat, &mvalid, &trigger, &xpos, &ypos);
    if (vstat == CVAL) {
      *x      = ndc_map_x(xpos);
      *y      = ndc_map_y(ypos);
      *charac = trigger;
    }
    else {
      vdicgi_errh(" SVDI Shell (VDAKGL) invalid request for crqlc ");
    }
  }
}

void vdstla_(float *x, float *y)
{
  int   temp, vstat, rstat, mvalid, trigger;
  float xpos, ypos, timeout;

  if (dev_cap[0] != 0.) {
    vflush(); /*  flush polyline buffer  */
    temp    = 1;
    timeout = 1.;
    crqlc_(&temp, &timeout, &vstat, &rstat, &mvalid, &trigger, &xpos, &ypos);
    if (vstat == CVAL) {
      *x = ndc_map_x(xpos);
      *y = ndc_map_y(ypos);
    }
    else {
      vdicgi_errh(" SVDI Shell (VDSTLA) invalid request for crqlc ");
    }
  }
}

void vdstos_(float attr_array[])
{
  int temp;

  temp = (int)attr_array[0];
  vdstfc_(&temp);
  temp = (int)attr_array[1];
  vdstbc_(&temp);
  vdstin_(attr_array + 2);
  temp = (int)attr_array[3];
  vdstls_(&temp);
  vdstlw_(attr_array + 4);
  vdstcs_(attr_array + 5);
}

void vdfram_(int itype) {}

void vbiqpk_() {}
void vberrh_() {}
