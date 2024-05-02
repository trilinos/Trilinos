/*
 * Copyright(C) 1999-2020, 2022, 2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 * SUN DEC/ULTRIX ALLIANT : C routines must have underscores
 * SGI CONVEX             : C routines must have underscores
 *
 * VAX HP IBM/aix         : C routines do not have underscores
 *
 * CRAY/UNICOS            : C routines must be capitalized,
 *                            and no underscores
 */

#if !defined(CRA) && !defined(ADDC_) && !defined(COUGAR)
#define vdmova wmetmv
#define vdmoni wmetmo
#define vdgnam wmetgn
#define vbiqdv wmetiv
#define vbiqpk wmetqp
#define vdlina wmetln
#define vdtext wmettx
#define vdpnta wmetpt
#define vdpoly wmetpy
#define vdiqcp wmetcp
#define vdstos wmetos
#define vdiqos wmetio
#define vdstfc wmetfc
#define vdstbc wmetbc
#define vdstin wmetin
#define vdstls wmetls
#define vdstlw wmetlw
#define vdstcs wmetcs
#define vdaabu wmetbu
#define vdaloc wmetlo
#define vdabgl wmetbl
#define vdakgl wmetkl
#define vdstla wmetla
#define vdinit wmetnt
#define vdfram wmetfr
#define vdterm wmettr
#define vdiqdc wmetdc
#define vdnwpg wmetpg
#define vdbell wmetbe
#define vdwait wmetwt
#define vdbufl wmetfl
#define vdstco wmetco
#define vdiqco wmetic
#define vdescp wmetes
#define vdiqes wmetie
#define vdiqnd wmetid
#define vimova wmetim
#define vilina wmetil
#define vipnta wmetip
#define vitext wmetix
#define viinit wmetii
#define viterm wmetit
#define vinwpg wmetig
#define vcjob  vcjob
#define vberrh wmeter
#define vdloge wmetle
#define cdrwfs wmetwf
#define cdrrfs wmetrf
#define cdrofs wmetof
#define cdrof3 wmeto3
#define cdrcfs wmetcf
#define cdroff wmetff
#define cdroab wmetab
#define bgpbuf wmetbf
#define qmsbuf wmetqm
#define qmsbu1 wmetbf
#define ddcbuf wmetbf
#define h75buf wmetbf
#define btkbuf wmetbf
#define nmtbuf wmetbf
#define vbimbf wmetib
#define vbpkg  wmetpk
#define vbdev  wmetdv
#define vdiqrs wmetqr
#define vdstmp wmetmp
#define vdstrs wmetrs
#define vdstrv wmetrv
#define vdbrgb wmetbg
#define vdfrgb wmetfg
#define vdpixl wmetpx
#define vdpixi wmetpi
#define vdrpix wmetrp
#define vdrpxi wmetri
#define vdrscl wmetrl
#define vdiqci wmetci
#define vbstmp wmet01
#define vifram wmet02
#define vcndcm wmet03
#define vcattr wmet04
#define vbini1 wmet05
#define vb2hls wmet06
#define vb2rgb wmet07
#define vccolt wmet08
#define vccrps wmet09
#define vcscal wmet10
#define vcddim wmet11
#define vipoly wmet12
#define vbout  wmet13
#define cgixxx cgimet
#endif
#if defined(CRA)
#define vdmova WMETMV
#define vdmoni WMETMO
#define vdgnam WMETGN
#define vbiqdv WMETIV
#define vbiqpk WMETQP
#define vdlina WMETLN
#define vdtext WMETTX
#define vdpnta WMETPT
#define vdpoly WMETPY
#define vdiqcp WMETCP
#define vdstos WMETOS
#define vdiqos WMETIO
#define vdstfc WMETFC
#define vdstbc WMETBC
#define vdstin WMETIN
#define vdstls WMETLS
#define vdstlw WMETLW
#define vdstcs WMETCS
#define vdaabu WMETBU
#define vdaloc WMETLO
#define vdabgl WMETBL
#define vdakgl WMETKL
#define vdstla WMETLA
#define vdinit WMETNT
#define vdfram WMETFR
#define vdterm WMETTR
#define vdiqdc WMETDC
#define vdnwpg WMETPG
#define vdbell WMETBE
#define vdwait WMETWT
#define vdbufl WMETFL
#define vdstco WMETCO
#define vdiqco WMETIC
#define vdescp WMETES
#define vdiqes WMETIE
#define vdiqnd WMETID
#define vimova WMETIM
#define vilina WMETIL
#define vipnta WMETIP
#define vitext WMETIX
#define viinit WMETII
#define viterm WMETIT
#define vinwpg WMETIG
#define cdrcom CDRCOM
#define vcjob  VCJOB
#define vconod VCONOD
#define vberrh WMETER
#define vdloge WMETLE
#define cdrwfs WMETWF
#define cdrrfs WMETRF
#define cdrofs WMETOF
#define cdrof3 WMETO3
#define cdrcfs WMETCF
#define cdroff WMETFF
#define cdroab WMETAB
#define bgpbuf WMETBF
#define qmsbuf WMETQM
#define qmsbu1 WMETBF
#define ddcbuf WMETBF
#define h75buf WMETBF
#define btkbuf WMETBF
#define nmtbuf WMETBF
#define vbimbf WMETIB
#define vbpkg  WMETPK
#define vbdev  WMETDV
#define vdiqrs WMETQR
#define vdstmp WMETMP
#define vdstrs WMETRS
#define vdstrv WMETRV
#define vdbrgb WMETBG
#define vdfrgb WMETFG
#define vdpixl WMETPX
#define vdpixi WMETPI
#define vdrpix WMETRP
#define vdrpxi WMETRI
#define vdrscl WMETRL
#define vdiqci WMETCI
#define vbstmp WMET01
#define vifram WMET02
#define vcndcm WMET03
#define vcattr WMET04
#define vbini1 WMET05
#define vb2hls WMET06
#define vb2rgb WMET07
#define vccolt WMET08
#define vccrps WMET09
#define vcscal WMET10
#define vcddim WMET11
#define vipoly WMET12
#define vbout  WMET13
#define wmetzz WMETZZ
#define cgixxx CGIMET
#endif
#if defined(ADDC_) || defined(COUGAR)
#define vdmova wmetmv_
#define vdmoni wmetmo_
#define vdgnam wmetgn_
#define vbiqdv wmetiv_
#define vbiqpk wmetqp_
#define vdlina wmetln_
#define vdtext wmettx_
#define vdpnta wmetpt_
#define vdpoly wmetpy_
#define vdiqcp wmetcp_
#define vdstos wmetos_
#define vdiqos wmetio_
#define vdstfc wmetfc_
#define vdstbc wmetbc_
#define vdstin wmetin_
#define vdstls wmetls_
#define vdstlw wmetlw_
#define vdstcs wmetcs_
#define vdaabu wmetbu_
#define vdaloc wmetlo_
#define vdabgl wmetbl_
#define vdakgl wmetkl_
#define vdstla wmetla_
#define vdinit wmetnt_
#define vdfram wmetfr_
#define vdterm wmettr_
#define vdiqdc wmetdc_
#define vdnwpg wmetpg_
#define vdbell wmetbe_
#define vdwait wmetwt_
#define vdbufl wmetfl_
#define vdstco wmetco_
#define vdiqco wmetic_
#define vdescp wmetes_
#define vdiqes wmetie_
#define vdiqnd wmetid_
#define vimova wmetim_
#define vilina wmetil_
#define vipnta wmetip_
#define vitext wmetix_
#define viinit wmetii_
#define viterm wmetit_
#define vinwpg wmetig_
#define cdrcom cdrcom_
#define vcjob  vcjob_
#define vconod vconod_
#define vberrh wmeter_
#define vdloge wmetle_
#define cdrwfs wmetwf_
#define cdrrfs wmetrf_
#define cdrofs wmetof_
#define cdrof3 wmeto3_
#define cdrcfs wmetcf_
#define cdroff wmetff_
#define cdroab wmetab_
#define bgpbuf wmetbf_
#define qmsbuf wmetqm_
#define qmsbu1 wmetbf_
#define ddcbuf wmetbf_
#define h75buf wmetbf_
#define btkbuf wmetbf_
#define nmtbuf wmetbf_
#define vbimbf wmetib_
#define vbpkg  wmetpk_
#define vbdev  wmetdv_
#define vdiqrs wmetqr_
#define vdstmp wmetmp_
#define vdstrs wmetrs_
#define vdstrv wmetrv_
#define vdbrgb wmetbg_
#define vdfrgb wmetfg_
#define vdpixl wmetpx_
#define vdpixi wmetpi_
#define vdrpix wmetrp_
#define vdrpxi wmetri_
#define vdrscl wmetrl_
#define vdiqci wmetci_
#define vbstmp wmet01_
#define vifram wmet02_
#define vcndcm wmet03_
#define vcattr wmet04_
#define vbini1 wmet05_
#define vb2hls wmet06_
#define vb2rgb wmet07_
#define vccolt wmet08_
#define vccrps wmet09_
#define vcscal wmet10_
#define vcddim wmet11_
#define vipoly wmet12_
#define vbout  wmet13_
#define wmetbf wmetbf_
#define wmetzz wmetzz_
#define cgixxx cgimet_
#endif
#if defined(CONVEX)
#define cdrcom_ _cdrcom_
#define cdrcm2_ _cdrcm2_
#define vcjob_  _vcjob_
#define vconod_ _vconod_
#define cdrunx_ _cdrunx_
#endif

/* cgimet.h - header file to define device dependent stuff for driver
 *   Metafile (met)
 * Sandia National Laboratories, Div 2634
 * Sun Nov 19 12:02:43 MST 1989 - last date modified
 */

#include "svdi.h"

/******************************************************************************/
/*                                                                            */
/*      type and constant declarations                                        */
/*                                                                            */
/******************************************************************************/

#define MAX_DEVICE_SURFACES 6   /* num of surfaces device supports */
                                /* ...set to 1 for interactive device */
#define ERROR_LIST_SIZE  10     /* number of errors stored */
#define COLOR_TABLE_SIZE 256    /* max color table for SVDI is 256 */
#define MAX_TEXT_LENGTH  136    /* max text length for SVDI is 136 */
#define MAX_POLYLINE     -1     /* -1 = no limit */
#define MAX_DJ_POLYLINE  -1     /* -1 = no limit */
#define MAX_POLYMARKER   -1     /* -1 = no limit */
#define MAX_POLYGON      2000   /* uses min of this and SVDI limit */
#define MAX_ARRAY        1      /* for pixel array and cell array */
                                /* ...set to 1 for non raster device */
#define BUFFER_SIZE 1440        /* for batch device buffer routines */
                                /* ...set to 1 for interactive device */
#define DEFAULT_OUTFILE_NAME "" /* not used, driver builds the name */

/* end cgimet.h */
#include "data_def.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#if defined(CRA) || (defined(SUN) && (defined(SYSV) || defined(SVR4)))
#include <fcntl.h>
#include <sys/types.h>
#endif
#include "cdrcom.h"
#include "cgi.h"
#include "cgidef.h"
#include "devid.h"
#include "fortyp.h"
#include "stdtyp.h"
#include <sys/file.h>

/* ------------------------------------------------------------*/
/* >> FUNCTION PROTOTYPES                                      */
/*-------------------------------------------------------------*/
/* cgi/vdi routines */
static void xactivate(anything **surf_list);
static void xdeactivate(anything **surf_list);
static void xci(anything **params, int num_surfaces, anything **surf_list);
static void xct(int num_surfaces, anything **surf_list);
static void xcxdfac(int num_surfaces, anything **surf_list);
static void xcpds(anything **params, int num_surfaces, anything **surf_list);
static void xcendpg(int num_surfaces, anything **surf_list);
static void xcbc(anything **params, int num_surfaces, anything **surf_list);
static void xcvdcx(anything **params, int num_surfaces, anything **surf_list);
static void xcv(anything **params, int num_surfaces, anything **surf_list);
static void xccl(anything **params, int num_surfaces, anything **surf_list);
static void xcdscl(anything **params, int num_surfaces, anything **surf_list);
static void xcdqerr(anything **params, anything **surf_list);
static void xcerhct(anything **params, int num_surfaces, anything **surf_list);
static void xccixp(anything **params, int num_surfaces, anything **surf_list);
static void xcesc(anything **params, int num_surfaces, anything **surf_list);
static void xcqid(anything **params, anything **surf_list);
static void xcqd(anything **params, anything **surf_list);
static void xclf(anything **params, anything **surf_list);
static void xclpr(anything **params, anything **surf_list);
static void xcqsp(anything **params, anything **surf_list);
static void xclesc(anything **params, anything **surf_list);
static void xcqp(anything **params, anything **surf_list);
static void xcqcl(anything **params, anything **surf_list);
static void xcpl(anything **params, int num_surfaces, anything **surf_list);
static void xcdjpl(anything **params, int num_surfaces, anything **surf_list);
static void xcpm(anything **params, int num_surfaces, anything **surf_list);
static void xctx(anything **params, int num_surfaces, anything **surf_list);
static void xcpg(anything **params, int num_surfaces, anything **surf_list);
static void xcca(anything **params, int num_surfaces, anything **surf_list);
static void xcpxa(anything **params, int num_surfaces, anything **surf_list);
static void xclnt(anything **params, int num_surfaces, anything **surf_list);
static void xclnw(anything **params, int num_surfaces, anything **surf_list);
static void xclnc(anything **params, int num_surfaces, anything **surf_list);
static void xcmkt(int num_surfaces, anything **surf_list);
static void xcmks(anything **params, int num_surfaces, anything **surf_list);
static void xcmkc(anything **params, int num_surfaces, anything **surf_list);
static void xctxp(anything **params, int num_surfaces, anything **surf_list);
static void xctxc(anything **params, int num_surfaces, anything **surf_list);
static void xcchh(anything **params, int num_surfaces, anything **surf_list);
static void xccho(anything **params, int num_surfaces, anything **surf_list);
static void xcis(anything **params, int num_surfaces, anything **surf_list);
static void xcflc(anything **params, int num_surfaces, anything **surf_list);
static void xccsm(anything **params, int num_surfaces, anything **surf_list);
static void xcct(anything **params, int num_surfaces, anything **surf_list);
static void xcgtxx(anything **params, anything **surf_list);
static void xcqprl(anything **params, anything **surf_list);
static void xcqln(anything **params, anything **surf_list);
static void xcqlnt(anything **params, anything **surf_list);
static void xcqslw(anything **params, anything **surf_list);
static void xcqmk(anything **params, anything **surf_list);
static void xcqmkt(anything **params, anything **surf_list);
static void xcqsms(anything **params, anything **surf_list);
static void xcqchh(anything **params, anything **surf_list);
static void xcqfl(anything **params, anything **surf_list);
static void xcqc(anything **params, anything **surf_list);
static void xcqlna(anything **params, anything **surf_list);
static void xcqmka(anything **params, anything **surf_list);
static void xcqtxa(anything **params, anything **surf_list);
static void xcqfla(anything **params, anything **surf_list);
static void xcqcte(anything **params, anything **surf_list);
static void xcili(anything **params, anything **surf_list);
static void xcrqlc(anything **params, anything **surf_list);
/* utility routines */
static void  init_state(surf_statelist *surf_state);
static void  set_dev_descrip(void);
static void  reset_vdi(surf_statelist *surf_state);
static void  set_mapping(surf_statelist *surf_state);
static void  set_clipping(surf_statelist *cur_state);
static void  set_foreground_color(surf_statelist *surf_state, int *colors);
static void  set_background_color(surf_statelist *surf_state, int *colors);
static void  report_error(surf_statelist *surf_state, int e_class, int e_num, int f_id);
static void  gettoken(int *index_p, int *numrecs_p, char *data_p, int max_chars, char *outtoken_p);
static int   poly_clip(point *cmin, point *cmax, float *vx, float *vy, int vlen, float *xout,
                       float *yout, int *lenout);
static int   inside_bnd(point *v, point *bmin, point *bmax, int bound_num);
static point intersect_bnd(point *p1, point *p2, point *bmin, point *bmax, int bound_num);
/* cdr routines */
/* -- not used with interactive devices */
void cdrofs(int *);
void cdrof3(int *, int *);
void cdroff(int *, int *, int *, int *);
void cdrrfs(int *, int *, char *, int *);
void cdrwfs(int *, int *, char *, int *);
void cdrcfs(int *, int *);
void cdroab(int *, int *);

static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

/*-------------------------------------------------------------*/
/* >> GLOBAL VARIABLE DECLARATIONS                             */
/*-------------------------------------------------------------*/

/* variables used to maintain state list, free list */
static int first_free_state = -1;
static int last_alloc_state = -1;

/* declare one statelist per surface */
static surf_statelist surf_states[MAX_DEVICE_SURFACES];

/* declare one description table per device */
static dev_descrip_table dev_descrip;

/* variable to keep track of current statelist */
/* -- this is global so CDR can know which state to look at */
static surf_statelist *cur_state;

/* These arrays are global because of size. Used in cell and pixel array */
static float rarray[MAX_ARRAY], garray[MAX_ARRAY], barray[MAX_ARRAY];
static int   varray[MAX_ARRAY];

/*-------------------------------------------------------------*/
/* >> CGI/SVDI DEVICE DRIVER ROUTINES                          */
/*-------------------------------------------------------------*/

void cgixxx(anything *params[], int num_surfaces, anything *surf_list[])
{
  /* Main entry point */
  switch (*(short *)params[0]) {
  case ACTIVATE_FN: xactivate(surf_list); break;
  case DEACTIVATE_FN: xdeactivate(surf_list); break;
  case CI_FN: xci(params, num_surfaces, surf_list); break;
  case CT_FN: xct(num_surfaces, surf_list); break;
  case CXDFAC_FN: xcxdfac(num_surfaces, surf_list); break;
  case CPDS_FN: xcpds(params, num_surfaces, surf_list); break;
  case CENDPG_FN: xcendpg(num_surfaces, surf_list); break;
  case CBC_FN: xcbc(params, num_surfaces, surf_list); break;
  case CVDCX_FN: xcvdcx(params, num_surfaces, surf_list); break;
  case CV_FN: xcv(params, num_surfaces, surf_list); break;
  case CCL_FN: xccl(params, num_surfaces, surf_list); break;
  case CDSCL_FN: xcdscl(params, num_surfaces, surf_list); break;
  case CDQERR_FN: xcdqerr(params, surf_list); break;
  case CERHCT_FN: xcerhct(params, num_surfaces, surf_list); break;
  case CCIXP_FN: xccixp(params, num_surfaces, surf_list); break;
  case CESC_FN: xcesc(params, num_surfaces, surf_list); break;
  case CQID_FN: xcqid(params, surf_list); break;
  case CQD_FN: xcqd(params, surf_list); break;
  case CLF_FN: xclf(params, surf_list); break;
#ifndef sgi
  case CLPR_FN: xclpr(params, surf_list); break;
#endif
  case CQSP_FN: xcqsp(params, surf_list); break;
  case CLESC_FN: xclesc(params, surf_list); break;
  case CQP_FN: xcqp(params, surf_list); break;
  case CQCL_FN: xcqcl(params, surf_list); break;
  case CPL_FN: xcpl(params, num_surfaces, surf_list); break;
  case CDJPL_FN: xcdjpl(params, num_surfaces, surf_list); break;
  case CPM_FN: xcpm(params, num_surfaces, surf_list); break;
  case CTX_FN: xctx(params, num_surfaces, surf_list); break;
  case CPG_FN: xcpg(params, num_surfaces, surf_list); break;
  case CCA_FN: xcca(params, num_surfaces, surf_list); break;
  case CPXA_FN: xcpxa(params, num_surfaces, surf_list); break;
  case CLNT_FN: xclnt(params, num_surfaces, surf_list); break;
  case CLNW_FN: xclnw(params, num_surfaces, surf_list); break;
  case CLNC_FN: xclnc(params, num_surfaces, surf_list); break;
  case CMKT_FN: xcmkt(num_surfaces, surf_list); break;
  case CMKS_FN: xcmks(params, num_surfaces, surf_list); break;
  case CMKC_FN: xcmkc(params, num_surfaces, surf_list); break;
  case CTXP_FN: xctxp(params, num_surfaces, surf_list); break;
  case CTXC_FN: xctxc(params, num_surfaces, surf_list); break;
  case CCHH_FN: xcchh(params, num_surfaces, surf_list); break;
  case CCHO_FN: xccho(params, num_surfaces, surf_list); break;
  case CIS_FN: xcis(params, num_surfaces, surf_list); break;
  case CFLC_FN: xcflc(params, num_surfaces, surf_list); break;
  case CCSM_FN: xccsm(params, num_surfaces, surf_list); break;
  case CCT_FN: xcct(params, num_surfaces, surf_list); break;
  case CGTXX_FN: xcgtxx(params, surf_list); break;
  case CQPRL_FN: xcqprl(params, surf_list); break;
  case CQLN_FN: xcqln(params, surf_list); break;
  case CQLNT_FN: xcqlnt(params, surf_list); break;
  case CQSLW_FN: xcqslw(params, surf_list); break;
  case CQMK_FN: xcqmk(params, surf_list); break;
  case CQMKT_FN: xcqmkt(params, surf_list); break;
  case CQSMS_FN: xcqsms(params, surf_list); break;
  case CQCHH_FN: xcqchh(params, surf_list); break;
  case CQFL_FN: xcqfl(params, surf_list); break;
  case CQC_FN: xcqc(params, surf_list); break;
  case CQLNA_FN: xcqlna(params, surf_list); break;
  case CQMKA_FN: xcqmka(params, surf_list); break;
  case CQTXA_FN: xcqtxa(params, surf_list); break;
  case CQFLA_FN: xcqfla(params, surf_list); break;
  case CQCTE_FN: xcqcte(params, surf_list); break;
  case CILI_FN: xcili(params, surf_list); break;
  case CRQLC_FN: xcrqlc(params, surf_list); break;
  } /* end switch */
} /* end cgixxx */

/* ACTIVATE */
static void xactivate(anything **surf_list)
{
  int new_state; /* index of new state to use */

  if (first_free_state >= 0) { /* there is free one, use it */
    new_state        = first_free_state;
    first_free_state = surf_states[new_state].next_free_state;
  }
  else if (last_alloc_state < MAX_DEVICE_SURFACES - 1) {
    new_state = last_alloc_state + 1; /* can alloc one */
    ++last_alloc_state;
  }
  else { /* none available */
    new_state = -1;
  } /* end if */

  if (new_state < 0) { /* problem, report error & go away */
    surf_list[0] = NULL;
    return;
  } /* end if problem */

  /* got new state, init it */
  surf_states[new_state].cgi_inited = CNO;
  surf_states[new_state].this_index = new_state; /* save for dealloc */
  /* -- for batch devices */
  copy_string(surf_states[new_state].filename, DEFAULT_OUTFILE_NAME, 100);
  surf_states[new_state].file_d = -1;

  /* set surface state_list pointer to point to this state list */
  surf_list[0] = (anything *)&(surf_states[new_state]);

} /* end xactivate */

/* DEACTIVATE */
static void xdeactivate(anything **surf_list)
{
  /* deallocate memory allocated for cgi state list */

  int index; /* which state is this surface linked to */

  index = ((surf_statelist *)surf_list[0])->this_index;
  if ((index >= 0) && (index < MAX_DEVICE_SURFACES)) {
    surf_states[index].next_free_state = first_free_state;
    first_free_state                   = index;
  } /* end if */

  /* unlink cgi state list from surface */
  surf_list[0] = NULL;

} /* end xdeactivate */

/* INITIALIZE */
static void xci(anything **params, int num_surfaces, anything **surf_list)
{
  int   i; /* index for loop on surf_list */
  float zero  = 0.;
  int   izero = 0;

  for (i = 0; i < num_surfaces; ++i) {

    /* cur_state is global - points to current statelist */
    cur_state = (surf_statelist *)surf_list[i];

    /* initialize/reinit SVDI and/or set device descript table */
    /*  -- this needs to happen before calling init_state */

    if (cur_state->cgi_inited == CNO) { /* this surface has never been
                                          initialized */
      vdinit(&zero, &izero);

      /* set device descript table for this device only once */
      if (dev_descrip.set_flag == FALSE) {
        set_dev_descrip();
        dev_descrip.set_flag = TRUE;
      } /* end if set_flag */
    } /* end if surface not initialized */

    else { /* this surface has been initialized once already */

      /* reset SVDI attributes */
      reset_vdi(cur_state);

      /* if prepare drawing surface, issue a new page */
      if (*(int *)params[1] == CYES) {
        vdnwpg();
      }
    } /* end else surface initialized */

    /* initialize/reinitialize the default device state */
    init_state(cur_state);

    /* set cgi_inited to initialized */
    cur_state->cgi_inited = CYES;

  } /* end for each surface */
} /* end xci */

/* TERMINATE */
static void xct(int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* terminate and set state to uninitialized */
    vdterm();
    cur_state->cgi_inited = CNO;

  } /* end for each surface */
} /* end xct */

/* EXECUTE DEFERRED ACTIONS */
static void xcxdfac(int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    vdbufl();

  } /* end for each surface */
} /* end xcxdfac */

/* PREPARE DRAWING SURFACE */
static void xcpds(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* new page if page has been marked or parameter is forced */
    if (cur_state->pic_dirty == CDIRTY || *(int *)params[1] == CFORCC) {
      vdnwpg();
      cur_state->pic_dirty = CCLEAN;
    }

  } /* end for each surface */
} /* end xcpds */

/* END PAGE */
static void xcendpg(int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* perform an implicit execute deferred actions */
    vdbufl();

    /* for hardcopy, call new page only if page has been marked */
    if (dev_descrip.copy_class == CHARD) {
      if (cur_state->pic_dirty == CDIRTY) {
        vdnwpg();
        cur_state->pic_dirty = CCLEAN;
      } /* end if dirty page */
    }

  } /* end for each surface */
} /* end xcendpg */

/* BACKGROUND COLOR */
static void xcbc(anything **params, int num_surfaces, anything **surf_list)
{
  int i;           /* index for loop on surfaces */
  int bg_color[3]; /* temporary array */

  bg_color[0] = *(int *)params[1];
  bg_color[1] = *(int *)params[2];
  bg_color[2] = *(int *)params[3];

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    set_background_color(cur_state, bg_color);

  } /* end for each surface */
} /* end xcbc */

/* VDC EXTENT */
static void xcvdcx(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* store the values in the state list */
    /* default integer precision is 16 - check range */
    if (*(float *)params[1] > 32767.0f || *(float *)params[2] > 32767.0f ||
        *(float *)params[3] > 32767.0f || *(float *)params[4] > 32767.0f ||
        *(float *)params[1] < -32767.0f || *(float *)params[2] < -32767.0f ||
        *(float *)params[3] < -32767.0f || *(float *)params[4] < -32767.0f) {
      /* error 1:-103 VDC Extent out of range.  Function ignored */
      report_error(cur_state, 1, -103, *(short *)params[0]);
      return;
    } /* end if error */

    /* no error */
    cur_state->vdc1.x = *(float *)params[1];
    cur_state->vdc1.y = *(float *)params[2];
    cur_state->vdc2.x = *(float *)params[3];
    cur_state->vdc2.y = *(float *)params[4];
    /* end else no error */

    /* VDC to viewport (in NDC) mapping */
    /* -- this routine sets effective viewport in VDC and
       NDC (forcing isotropic mapping) and the mapping scale */
    set_mapping(cur_state);

    /* determine new effective clipping rectangle */
    /* ..First figure out the min and max VDC points.  Since
     * the user cannot change the clip rectangle, the min and
     * max clip rectangle points are known. (min = clip_rect1 = [0,0];
     * max = clip_rect2 = [32767.,32767.])
     */
    /* need to check for degenerate case */
    cur_state->eff_clip_rect1.x =
        max(min(cur_state->vdc1.x, cur_state->vdc2.x), cur_state->clip_rect1.x);
    cur_state->eff_clip_rect1.y =
        max(min(cur_state->vdc1.y, cur_state->vdc2.y), cur_state->clip_rect1.y);
    cur_state->eff_clip_rect2.x =
        min(max(cur_state->vdc1.x, cur_state->vdc2.x), cur_state->clip_rect2.x);
    cur_state->eff_clip_rect2.y =
        min(max(cur_state->vdc1.y, cur_state->vdc2.y), cur_state->clip_rect2.y);

    /* set clipping region */
    set_clipping(cur_state);

  } /* end for each surface */
} /* end xcvdcx */

/* DEVICE VIEWPORT */
static void xcv(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* viewport  - default viewport specification mode is FRACTION OF
     * DRAWING SURFACE.  Thus, the viewport (0,0),(1,1) corresponds to
     * the lower and upper corners, respectively, of the default
     * device viewport, which is the largest unrotated rectangular area
     * visible on the drawing surface. Numbers outside the range
     * [0.0...1.0] may be specified.
     */

    /* store values in state list */
    cur_state->vp1.x = *(float *)params[1];
    cur_state->vp1.y = *(float *)params[2];
    cur_state->vp2.x = *(float *)params[3];
    cur_state->vp2.y = *(float *)params[4];

    /* do VDC to NDC mapping */
    /* -- this routine sets effective viewports in VDC and
       NDC (forcing isotropic mapping) and the mapping scale */
    set_mapping(cur_state);

    /* set clipping region */
    set_clipping(cur_state);

  } /* end for each surface */
} /* end xcv */

/* CLIP INDICATOR */
static void xccl(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* check for a legal enumerated value, store in state table
       and set the clipping region */
    if (*(int *)params[1] == CON || *(int *)params[1] == COFF) {
      cur_state->clip_indicator = *(int *)params[1];
      set_clipping(cur_state);
    }
    else {
      /* error 1:-101 Enumerated parameter out of range. Function ignored */
      report_error(cur_state, 1, -101, *(short *)params[0]);
    }

  } /* end for each surface */
} /* end xccl */

/* DRAWING SURFACE CLIP INDICATOR */
static void xcdscl(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* check for a legal enumerated value, store in state table
       and set the clipping region */
    if (*(int *)params[1] == CDCOFF || *(int *)params[1] == CDCREC || *(int *)params[1] == CVPORT) {
      cur_state->ds_clip_indicator = *(int *)params[1];
      set_clipping(cur_state);
    }
    else {
      /* error 1:-101 Enumerated parameter out of range. Function ignored */
      report_error(cur_state, 1, -101, *(short *)params[0]);
    }

  } /* end for each surface */
} /* end xcdscl */

/* DEQUEUE ERROR REPORTS */
static void xcdqerr(anything **params, anything **surf_list)
{
  int i;       /* loop variable */
  int max_ret; /* max number of errors to return */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* find the number of errors to be returned */
  max_ret = (*(int *)params[1] < cur_state->err_count) ? *(int *)params[1] : cur_state->err_count;

  for (i = 0; i < max_ret; i++) {

    ((int *)params[5])[i] = cur_state->err_queue[cur_state->err_head_ptr].err_class;
    ((int *)params[6])[i] = cur_state->err_queue[cur_state->err_head_ptr].err_num;
    ((int *)params[7])[i] = cur_state->err_queue[cur_state->err_head_ptr].func_id;

    /* update the head of list ptr */
    cur_state->err_head_ptr++;
    if (cur_state->err_head_ptr == ERROR_LIST_SIZE) {
      cur_state->err_head_ptr = 0;
    }

    /* decrement the count of errors */
    cur_state->err_count--;

  } /* end for i */

  /* number of reports remaining in the queue */
  *(int *)params[3] = cur_state->err_count;

  /* number of reports returned */
  *(int *)params[4] = i;

  /* validity flag */
  *(int *)params[2] = CVAL;

} /* end xcdqerr */

/* ERROR HANDLING CONTROL */
static void xcerhct(anything **params, int num_surfaces, anything **surf_list)
{
  int i, j; /* loop indices */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* store error control level in state list */
    for (j = 0; j < *(int *)params[1]; j++) {

      /* check enumeration types */
      if (((int *)params[3])[j] != CEHON && ((int *)params[3])[j] != CEHROF &&
          ((int *)params[3])[j] != CEHDOF) {
        /* error 1:-101 Enumerated parameter out of range. Function ignored */
        report_error(cur_state, 1, -101, *(short *)params[0]);

        /* check for legal index */
      }
      else if (((int *)params[2])[j] < 1 || ((int *)params[2])[j] > MAX_ERROR_CLASS) {
        /* error 1:-102 Index out of range. Function ignored */
        report_error(cur_state, 1, -102, *(short *)params[0]);
      }
      else {
        cur_state->err_flag[((int *)params[2])[j] - 1] = ((int *)params[3])[j];
      }
    } /* end for j */

  } /* end for each surface */
} /* end xcerhct */

/* COLOR INDEX PRECISION */
static void xccixp(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  /* I'm not doing anything with color index precision now */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    if (*(int *)params[1] == 8 || *(int *)params[1] == 16 || *(int *)params[1] == 24 ||
        *(int *)params[1] == 32) {
      cur_state->color_index_prec = *(int *)params[1];
    }
    else {
      /* error 3:204 Selected precision not supported */
      report_error(cur_state, 3, 204, *(short *)params[0]);
    }

  } /* end for each surface */
} /* end xccixp */

/* ESCAPE */
static void xcesc(anything **params, int num_surfaces, anything **surf_list)
{
  int   i;          /* index for loop on surfaces */
  int   j;          /* loop index */
  int   first;      /* logical; TRUE = first time through */
  char  data[80];   /* temp for parsing data record */
  char  c;          /* single character */
  int   vdi_esc;    /* used for SVDI escape call */
  int   count;      /* used for SVDI escape call */
  float args[10];   /* used for SVDI escape call */
  int   tokenindex; /* index in data record of first char */
                    /*     after end of token*/
  first = TRUE;

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];

    switch (*(int *)params[1]) {

    case XEMFNM: /* file name */

      if (first) {
        tokenindex = 0; /* first char of data record */

        /* params[4] is an added parameter - contains the string length */
        gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);

        first = FALSE;
      } /* end if first */

      /* this function must be called before CI */
      /* ...error message if not?? check for legal file name?? */
      if (cur_state->cgi_inited != CYES) {
        copy_string(cur_state->filename, data, 100);
      }

      break;

    case XEAGMD: /* set terminal to alpha/graphics */

      if (first) {
        tokenindex = 0; /* first char of data record */

        /* params[4] is an added parameter - contains the string length */
        gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);

        /* SVDI takes care of 'graphics', check for 'alpha' */
        j = 0;
        while ((c = data[j++])) {
          if (isalpha(c) && isupper(c)) {
            data[j - 1] = tolower(c);
          }
        }

        first = FALSE;
      } /* end if first */

      if (strncmp(data, "alpha", 5) == 0) {
        vdbufl();
      }

      break;

    case XESVDI: /* SVDI escape */

      if (first) {
        tokenindex = 0; /* first char of data record */

        /* first token is SVDI escape number */
        gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);
        vdi_esc = strtol(data, NULL, 10);

        /* second token is number of parameters to SVDI escape */
        gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);
        count = strtol(data, NULL, 10);

        /* the rest of the tokens contain the arguments */
        for (j = 0; j < count && j < 10; j++) {
          gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);
          args[j] = (float)strtod(data, NULL);
        }

        first = FALSE;
      } /* end if first */

      vdescp(&vdi_esc, &count, args);

      break;

    case XEPGSZ: /* Page size */

      if (first) {
        tokenindex = 0; /* first char of data record */

        /* first token is page size in x */
        gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);
        args[0] = (float)strtod(data, NULL);

        /* second token is page size in y */
        gettoken(&tokenindex, (int *)params[2], (char *)params[3], 80, data);
        args[1] = (float)strtod(data, NULL);

        /* page size svdi number is 1400 */
        vdi_esc = 1400;
        count   = 2;

        first = FALSE;
      } /* end if first */

      vdescp(&vdi_esc, &count, args);

      break;

    case XEBELL: /* Ring the bell */ vdbell(); break;

    default:
      /* error 2:001 Request unsupported feature. Function ignored */
      report_error(cur_state, 2, 001, *(short *)params[0]);
      break;

    } /* end switch */
  } /* end for i */
} /* end xcesc */

/* INQUIRE DEVICE IDENTIFICATION */
static void xcqid(anything **params, anything **surf_list)
{
  char      *cgi_devid;      /* pointer to character device id */
  int        maxchr;         /* max number of chars in device id */
  int        qdc_index = 23; /* index for inquiries to vdiqdc */
  float      value;          /* value returned by vdiqdc */
  static int set = FALSE;    /* flag whether values have been set yet */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* device id */
    vdiqdc(&qdc_index, &value);
    cgi_devid = get_devid_char(value);
    if (cgi_devid != NULL) {
      copy_string(dev_descrip.dev_id, cgi_devid, 4);
    }

  } /* end if not set */

  /* response validity */
  *(int *)params[2] = CVAL;

  /* return device class - this is set in set_dev_descrip() */
  *(int *)params[3] = dev_descrip.dev_class;

  /* return device id */
  maxchr = (*(int *)params[1] > 3) ? 3 : *(int *)params[1];
  copy_string((char *)params[4], dev_descrip.dev_id, maxchr);
  *((char *)params[4] + maxchr) = '\0';

} /* end xcqid */

/* INQUIRE DEVICE DESCRIPTION */
static void xcqd(anything **params, anything **surf_list)
{
  int        qdc_index;   /* index for inquiries to vdiqdc  */
  float      value;       /* value returned by vdiqdc */
  static int set = FALSE; /* flag whether values have been set */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* display type */
    qdc_index = 2;
    vdiqdc(&qdc_index, &value);
    if (value == 0.0f) {
      dev_descrip.display_type = CVECT;
    }
    else if (value == 1.0f) {
      dev_descrip.display_type = CRAST;
    }
    else {
      dev_descrip.display_type = COTHER;
    }

    /* background color capability */
    qdc_index = 32;
    vdiqdc(&qdc_index, &value);
    dev_descrip.bcolor_cap = (value == 1.0f) ? CYES : CNO;

    /* dynamic modification for background color */
    dev_descrip.dynamic_mod_bg = CIRG;

    /* dynamic modification for VDC to device mapping */
    dev_descrip.dynamic_mod_map = CIRG;

    /* width and height in millimeters */
    /* ...SVDI returns 0 if it doesn't know. CGI does this also */
    qdc_index = 17;
    vdiqdc(&qdc_index, &dev_descrip.draw_surf_width);
    qdc_index = 18;
    vdiqdc(&qdc_index, &dev_descrip.draw_surf_height);

    /* pixel location relative to coordinates */
    dev_descrip.pix_loc = CPXON;

  } /* end haven't been set */

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* stuff info into return values */
  /* ...copy_class and device corners are set in set_dev_descrip() */
  *(int *)params[2]    = dev_descrip.copy_class;
  *(int *)params[3]    = dev_descrip.display_type;
  *(int *)params[4]    = dev_descrip.bcolor_cap;
  *(int *)params[5]    = dev_descrip.dynamic_mod_bg;
  *(int *)params[6]    = dev_descrip.dynamic_mod_map;
  *(int *)params[7]    = dev_descrip.dc1.x;
  *(int *)params[8]    = dev_descrip.dc1.y;
  *(int *)params[9]    = dev_descrip.dc2.x;
  *(int *)params[10]   = dev_descrip.dc2.y;
  *(float *)params[11] = dev_descrip.draw_surf_width;
  *(float *)params[12] = dev_descrip.draw_surf_height;
  *(int *)params[13]   = dev_descrip.pix_loc;

} /* end xcqd */

/* LOOKUP FUNCTION SUPPORT */
static void xclf(anything **params, anything **surf_list)
{
  int          i;                  /* index for loop on surfaces */
  int          qdc_index    = 24;  /* for inquire SVDI */
  static float poly_support = -1.; /* level of SVDI poly support */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  for (i = 0; i < *(int *)params[1]; i++) {

    switch (*(short *)params[0]) {
    case ACTIVATE_FN: ((int *)params[5])[i] = CSYES; break;
    case DEACTIVATE_FN: ((int *)params[5])[i] = CSYES; break;
    case CI_FN: ((int *)params[5])[i] = CSYES; break;
    case CT_FN: ((int *)params[5])[i] = CSYES; break;
    case CXDFAC_FN: ((int *)params[5])[i] = CSYES; break;
    case CPDS_FN: ((int *)params[5])[i] = CSYES; break;
    case CENDPG_FN: ((int *)params[5])[i] = CSYES; break;
    case CBC_FN: ((int *)params[5])[i] = CSYES; break;
    case CVDCX_FN: ((int *)params[5])[i] = CSYES; break;
    case CV_FN: ((int *)params[5])[i] = CSYES; break;
    case CCL_FN: ((int *)params[5])[i] = CSYES; break;
    case CDSCL_FN: ((int *)params[5])[i] = CSYES; break;
    case CDQERR_FN: ((int *)params[5])[i] = CSYES; break;
    case CERHCT_FN: ((int *)params[5])[i] = CSYES; break;
    case CCIXP_FN: ((int *)params[5])[i] = CSNO; break;
    case CESC_FN: ((int *)params[5])[i] = CSYES; break;
    case CQID_FN: ((int *)params[5])[i] = CSYES; break;
    case CQD_FN: ((int *)params[5])[i] = CSYES; break;
    case CLF_FN: ((int *)params[5])[i] = CSYES; break;
    case CLPR_FN: ((int *)params[5])[i] = CSYES; break;
    case CQSP_FN: ((int *)params[5])[i] = CSYES; break;
    case CLESC_FN: ((int *)params[5])[i] = CSYES; break;
    case CQP_FN: ((int *)params[5])[i] = CSYES; break;
    case CQCL_FN: ((int *)params[5])[i] = CSYES; break;
    case CPL_FN: ((int *)params[5])[i] = CSYES; break;
    case CDJPL_FN: ((int *)params[5])[i] = CSYES; break;
    case CPM_FN: ((int *)params[5])[i] = CSYES; break;
    case CTX_FN: ((int *)params[5])[i] = CSYES; break;
    case CPG_FN:
      if (poly_support == -1.0f) {
        vdiqdc(&qdc_index, &poly_support);
      }
      if (poly_support < 3.0f) {
        ((int *)params[5])[i] = CSNO;
      }
      else {
        ((int *)params[5])[i] = CSYES;
      }
      break;
    case CCA_FN: ((int *)params[5])[i] = CSYES; break;
    case CPXA_FN: ((int *)params[5])[i] = CSYES; break;
    case CLNT_FN: ((int *)params[5])[i] = CSYES; break;
    case CLNW_FN: ((int *)params[5])[i] = CSYES; break;
    case CLNC_FN: ((int *)params[5])[i] = CSYES; break;
    case CMKT_FN: ((int *)params[5])[i] = CSYES; break;
    case CMKS_FN: ((int *)params[5])[i] = CSNO; break;
    case CMKC_FN: ((int *)params[5])[i] = CSYES; break;
    case CTXP_FN: ((int *)params[5])[i] = CSNO; break;
    case CTXC_FN: ((int *)params[5])[i] = CSNO; break;
    case CCHH_FN: ((int *)params[5])[i] = CSYES; break;
    case CCHO_FN: ((int *)params[5])[i] = CSNO; break;
    case CIS_FN: ((int *)params[5])[i] = CSYES; break;
    case CFLC_FN: ((int *)params[5])[i] = CSYES; break;
    case CCSM_FN: ((int *)params[5])[i] = CSYES; break;
    case CCT_FN: ((int *)params[5])[i] = CSYES; break;
    case CGTXX_FN: ((int *)params[5])[i] = CSYES; break;
    case CQPRL_FN: ((int *)params[5])[i] = CSYES; break;
    case CQLN_FN: ((int *)params[5])[i] = CSYES; break;
    case CQLNT_FN: ((int *)params[5])[i] = CSYES; break;
    case CQSLW_FN: ((int *)params[5])[i] = CSNO; break;
    case CQMK_FN: ((int *)params[5])[i] = CSNO; break;
    case CQMKT_FN: ((int *)params[5])[i] = CSYES; break;
    case CQSMS_FN: ((int *)params[5])[i] = CSNO; break;
    case CQCHH_FN: ((int *)params[5])[i] = CSYES; break;
    case CQFL_FN: ((int *)params[5])[i] = CSYES; break;
    case CQC_FN: ((int *)params[5])[i] = CSYES; break;
    case CQLNA_FN: ((int *)params[5])[i] = CSYES; break;
    case CQMKA_FN: ((int *)params[5])[i] = CSNO; break;
    case CQTXA_FN: ((int *)params[5])[i] = CSYES; break;
    case CQFLA_FN: ((int *)params[5])[i] = CSNO; break;
    case CQCTE_FN: ((int *)params[5])[i] = CSYES; break;
    case CILI_FN: ((int *)params[5])[i] = CSYES; break;
    case CRQLC_FN: ((int *)params[5])[i] = CSYES; break;
    default: ((int *)params[5])[i] = CSNO; break;
    } /* end switch */
  }

  /* validity flag */
  *(int *)params[4] = CVAL;

} /* end xclf */

/* LOOKUP PROFILE SUPPORT */
static void xclpr(anything **params, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  for (i = 0; i < *(int *)params[1]; i++) {
    if (strncmp(((char **)params[2])[i], "1-WAY-OUTPUT", 12) != 0 ||
        strncmp(((char **)params[2])[i], "2-WAY-OUTPUT", 12) ||
        strncmp(((char **)params[2])[i], "SLATEC 1", 8) != 0) {
      ((int *)params[4])[i] = CSYES;
    }
    else {
      ((int *)params[4])[i] = CSNO;
    }
  }
  /* validity flag */
  *(int *)params[3] = CVAL;

} /* end xclpr */

/* INQUIRE SUPPORTED PRECISIONS */
static void xcqsp(anything **params, anything **surf_list)
{
  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* this stuff is basically machine dependent. I return INVAL and
     leave it at that */
  *(int *)params[1] = CINVAL;

} /* end xcqsp */

/* LOOKUP ESCAPE SUPPORT */
static void xclesc(anything **params, anything **surf_list)
{
  int i;          /* loop index */
  int esc = 1400; /* for inquire escape */
  int support;    /* for inquire escape */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  for (i = 0; i < *(int *)params[1]; i++) {

    /* return YES for supported escapes, else return NO */
    switch (((int *)params[2])[i]) {
    case XEMFNM: ((int *)params[4])[i] = CSYES; break;
    case XEAGMD: ((int *)params[4])[i] = CSYES; break;
    case XESVDI: ((int *)params[4])[i] = CSYES; break;
    case XEBELL: ((int *)params[4])[i] = CSYES; break;
    case XEPGSZ:
      vdiqes(&esc, &support);
      if (support == 1 || support == 2) {
        ((int *)params[4])[i] = CSYES;
      }
      else {
        ((int *)params[4])[i] = CSNO;
      }
      break;

    default: ((int *)params[4])[i] = CSNO; break;
    } /* end switch */
  }

  /* validity flag */
  *(int *)params[3] = CVAL;

} /* end xclesc */

/* INQUIRE CURRENT PRECISIONS */
static void xcqp(anything **params, anything **surf_list)
{
  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* validity flag */
  *(int *)params[1] = CINVAL;

} /* end xcqp */

/* INQUIRE CLIPPING */
static void xcqcl(anything **params, anything **surf_list)
{
  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* clip indicator */
  *(int *)params[2] = cur_state->clip_indicator;

  /* clip rectangle */
  ((float *)params[3])[0] = cur_state->clip_rect1.x;
  ((float *)params[3])[1] = cur_state->clip_rect1.y;
  ((float *)params[3])[2] = cur_state->clip_rect2.x;
  ((float *)params[3])[3] = cur_state->clip_rect2.y;

  /* drawing surface clip indicator */
  *(int *)params[4] = cur_state->ds_clip_indicator;

  /* drawing surface clip rectangle */
  /* Default is the whole drawing surface; routines to change the
   * drawing surface clip rectangle are not included. Since veiwport
   * specification mode is FRACTION OF DRAWING SURFACE, units for
   * rectangle are 0 - 1 */
  ((float *)params[5])[0] = 0.;
  ((float *)params[5])[1] = 0.;
  ((float *)params[5])[2] = 1.;
  ((float *)params[5])[3] = 1.;

  /* validity flag */
  *(int *)params[1] = CVAL;

} /* end xcqcl */

/* POLYLINE */
static void xcpl(anything **params, int num_surfaces, anything **surf_list)
{
  int             i, j;                /* indices for loops */
  int             np;                  /* number of points in polyline */
  float           prev_x, prev_y;      /* line end points */
  float           cur_x, cur_y;        /* line end points */
  float           save_x, save_y;      /* used for clipping */
  int             prev_code, cur_code; /* encoded endpoints - for clipping */
  int             mode = 0, done;      /* stuff used for clipping */
  static unsigned mask = ~(~0u << 1);  /* for masking off bits */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* do some error checking */
    if (MAX_POLYLINE != -1 && *(int *)params[1] > MAX_POLYLINE) {
      /* error 6:301 - Number of points too large. Max is used */
      report_error(cur_state, 6, 301, *(short *)params[0]);
      np = MAX_POLYLINE;
    }
    else {
      np = *(int *)params[1];
    }

    /* set SVDI foreground color if needed */
    set_foreground_color(cur_state, cur_state->line_color);

    /* take a short cut if no clipping is on */
    if (!cur_state->clip_on) {
      for (j = 0; j < np; j++) {

        /* map VDC to NDC */
        cur_x = *((float *)params[2] + j) * cur_state->xscale + cur_state->xoffset;
        cur_y = *((float *)params[3] + j) * cur_state->yscale + cur_state->yoffset;
        if (j == 0) {
          vdmova(&cur_x, &cur_y);
        }
        else {
          vdlina(&cur_x, &cur_y);
        }
      } /* end for j */
    } /* end no clipping on */

    else { /* clipping is on */

      /* This clipping algorithm is loosely based on the Cohen-
       * Sutherland algorithm.  This algorithm takes into account
       * the fact that points are consecutive, ie., it tries to
       * encode each point only once. Everything is done
       * inline for efficiency
       */

      /* convert VDC to NDC */
      prev_x = *((float *)params[2]) * cur_state->xscale + cur_state->xoffset;
      prev_y = *((float *)params[3]) * cur_state->yscale + cur_state->yoffset;

      /* encode the point */
      prev_code = 0;
      if ((cur_state->clipmax.y - prev_y) < 0) {
        prev_code = prev_code | (1 << 3);
      }
      else if ((prev_y - cur_state->clipmin.y) < 0) {
        prev_code = prev_code | (1 << 2);
      }

      if ((cur_state->clipmax.x - prev_x) < 0) {
        prev_code = prev_code | (1 << 1);
      }
      else if ((prev_x - cur_state->clipmin.x) < 0) {
        prev_code = prev_code | (1 << 0);
      }

      if (!prev_code) {
        vdmova(&prev_x, &prev_y);
      }

      j = 0;
      while (++j < np) {

        /* convert VDC to NDC */
        cur_x = *((float *)params[2] + j) * cur_state->xscale + cur_state->xoffset;
        cur_y = *((float *)params[3] + j) * cur_state->yscale + cur_state->yoffset;

        /* encode the endpoint */
        cur_code = 0;
        if ((cur_state->clipmax.y - cur_y) < 0) {
          cur_code = cur_code | (1 << 3);
        }
        else if ((cur_y - cur_state->clipmin.y) < 0) {
          cur_code = cur_code | (1 << 2);
        }

        if ((cur_state->clipmax.x - cur_x) < 0) {
          cur_code = cur_code | (1 << 1);
        }
        else if ((cur_x - cur_state->clipmin.x) < 0) {
          cur_code = cur_code | (1 << 0);
        }

        /* determine which case this is */
        if (!prev_code && !cur_code) { /* trivial accept */
          vdlina(&cur_x, &cur_y);
          prev_x = cur_x;
          prev_y = cur_y;
          done   = TRUE;
        } /* end trivial accept */

        else if (prev_code & cur_code) { /* trivial reject */
          prev_x    = cur_x;
          prev_y    = cur_y;
          prev_code = cur_code;
          done      = TRUE;
        }

        else {
          done = FALSE;
          if (prev_code) {
            if (cur_code) {
              mode = 2; /* both out */
            }
            else {
              mode = 1; /* prev out; cur in */
            }
          }
          else {
            mode = 3; /* prev in; cur out */
          }
        }

        while (!done) {

          switch (mode) {

          case 1: /* prev out, cur in... clip */
          case 2: /* both out, clip prev first */

            if (prev_code & mask) { /* clip at left edge */
              prev_y =
                  prev_y + (cur_y - prev_y) * (cur_state->clipmin.x - prev_x) / (cur_x - prev_x);
              prev_x = cur_state->clipmin.x;
            }
            else if ((prev_code >> 1) & mask) { /* clip at right edge */
              prev_y =
                  prev_y + (cur_y - prev_y) * (cur_state->clipmax.x - prev_x) / (cur_x - prev_x);
              prev_x = cur_state->clipmax.x;
            }

            if ((prev_code >> 2) & mask) { /* clip at bottom edge */
              prev_x =
                  prev_x + (cur_x - prev_x) * (cur_state->clipmin.y - prev_y) / (cur_y - prev_y);
              prev_y = cur_state->clipmin.y;
            }
            else if ((prev_code >> 3) & mask) { /* clip at top edge */
              prev_x =
                  prev_x + (cur_x - prev_x) * (cur_state->clipmax.y - prev_y) / (cur_y - prev_y);
              prev_y = cur_state->clipmax.y;
            }

            switch (mode) {
            case 1: /* not looping through, draw it */
              vdmova(&prev_x, &prev_y);
              vdlina(&cur_x, &cur_y);
              prev_x    = cur_x;
              prev_y    = cur_y;
              prev_code = cur_code;
              done      = TRUE;
              break;

            case 2: /* both were out, encode new point */
              prev_code = 0;
              if ((cur_state->clipmax.y - prev_y) < 0) {
                prev_code = prev_code | (1 << 3);
              }
              else if ((prev_y - cur_state->clipmin.y) < 0) {
                prev_code = prev_code | (1 << 2);
              }

              if ((cur_state->clipmax.x - prev_x) < 0) {
                prev_code = prev_code | (1 << 1);
              }
              else if ((prev_x - cur_state->clipmin.x) < 0) {
                prev_code = prev_code | (1 << 0);
              }

              if (prev_code) { /* reject it */
                prev_x    = cur_x;
                prev_y    = cur_y;
                prev_code = cur_code;
                done      = TRUE;
              }
              else { /* prev is in */
                vdmova(&prev_x, &prev_y);
                mode = 3; /* cur is still out, go through again */
              }
              break;
            } /* end inner switch */
            break;

          case 3: /* prev in, cur out...clip */
            save_x = cur_x;
            save_y = cur_y;

            if (cur_code & mask) { /* clip at left edge */
              cur_y = cur_y + (prev_y - cur_y) * (cur_state->clipmin.x - cur_x) / (prev_x - cur_x);
              cur_x = cur_state->clipmin.x;
            }
            else if ((cur_code >> 1) & mask) { /* clip at right edge */
              cur_y = cur_y + (prev_y - cur_y) * (cur_state->clipmax.x - cur_x) / (prev_x - cur_x);
              cur_x = cur_state->clipmax.x;
            }

            if ((cur_code >> 2) & mask) { /* clip at bottom edge */
              cur_x = cur_x + (prev_x - cur_x) * (cur_state->clipmin.y - cur_y) / (prev_y - cur_y);
              cur_y = cur_state->clipmin.y;
            }
            else if ((cur_code >> 3) & mask) { /* clip at top edge */
              cur_x = cur_x + (prev_x - cur_x) * (cur_state->clipmax.y - cur_y) / (prev_y - cur_y);
              cur_y = cur_state->clipmax.y;
            }

            vdlina(&cur_x, &cur_y);
            prev_x    = save_x;
            prev_y    = save_y;
            prev_code = cur_code;
            done      = TRUE;
            break;

          default: done = TRUE; break;
          } /* end switch */
        } /* end while !done */
      } /* end while j*/

    } /* end else clipping is on */

    /* flag that the page has been marked */
    cur_state->pic_dirty = CDIRTY;

  } /* end for i */
} /* end xcpl */

/* DISJOINT POLYLINE */
static void xcdjpl(anything **params, int num_surfaces, anything **surf_list)
{
  int             i, j;               /* indices for loops */
  int             np;                 /* number of points in polyline */
  float           x1, y1, x2, y2;     /* line endpoints in NDC */
  int             code1, code2;       /* encoded endpoints for clipping */
  int             mode = 0, done;     /* variables used for clipping */
  static unsigned mask = ~(~0u << 1); /* for masking off bits */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* check for errors */
    if (MAX_DJ_POLYLINE != -1 && *(int *)params[1] > MAX_DJ_POLYLINE) {
      /* error 6:301 - Number of points too large. Max is used */
      report_error(cur_state, 6, 301, *(short *)params[0]);
      np = MAX_DJ_POLYLINE;
    }
    else {
      np = *(int *)params[1];
    }

    /* set SVDI foreground color if needed */
    set_foreground_color(cur_state, cur_state->line_color);

    /* if clipping is off, just do it */
    if (!cur_state->clip_on) {

      for (j = 0; j < np - 1; j = j + 2) {

        /* map VDC to NDC */
        x1 = *((float *)params[2] + j) * cur_state->xscale + cur_state->xoffset;
        y1 = *((float *)params[3] + j) * cur_state->yscale + cur_state->yoffset;
        x2 = *((float *)params[2] + j + 1) * cur_state->xscale + cur_state->xoffset;
        y2 = *((float *)params[3] + j + 1) * cur_state->yscale + cur_state->yoffset;

        vdmova(&x1, &y1);
        vdlina(&x2, &y2);

      } /* end for j */
    } /* end if clipping is off */

    else { /* clipping is on */

      /* Clip a line.  Kind of based on Cohen-Sutherland algorithm.
       * Everything is done inline for efficiency.
       */

      for (j = 0; j < np - 1; j = j + 2) {

        /* map VDC to NDC */
        x1 = *((float *)params[2] + j) * cur_state->xscale + cur_state->xoffset;
        y1 = *((float *)params[3] + j) * cur_state->yscale + cur_state->yoffset;
        x2 = *((float *)params[2] + j + 1) * cur_state->xscale + cur_state->xoffset;
        y2 = *((float *)params[3] + j + 1) * cur_state->yscale + cur_state->yoffset;

        /* encode the two points */
        code1 = 0;
        if ((cur_state->clipmax.y - y1) < 0) {
          code1 = code1 | (1 << 3);
        }
        else if ((y1 - cur_state->clipmin.y) < 0) {
          code1 = code1 | (1 << 2);
        }

        if ((cur_state->clipmax.x - x1) < 0) {
          code1 = code1 | (1 << 1);
        }
        else if ((x1 - cur_state->clipmin.x) < 0) {
          code1 = code1 | (1 << 0);
        }

        code2 = 0;
        if ((cur_state->clipmax.y - y2) < 0) {
          code2 = code2 | (1 << 3);
        }
        else if ((y2 - cur_state->clipmin.y) < 0) {
          code2 = code2 | (1 << 2);
        }

        if ((cur_state->clipmax.x - x2) < 0) {
          code2 = code2 | (1 << 1);
        }
        else if ((x2 - cur_state->clipmin.x) < 0) {
          code2 = code2 | (1 << 0);
        }

        /* determine which case this is */
        if (!code1 && !code2) { /* trivial accept - draw it */
          vdmova(&x1, &y1);
          vdlina(&x2, &y2);
          done = TRUE;
        } /* end trivial accept */

        else if (code1 & code2) { /* trivial reject - do nothing */
          done = TRUE;
        }
        else {
          done = FALSE;
          if (code1) {
            if (code2) {
              mode = 2; /* both out */
            }
            else {
              mode = 1; /* x1,y1 out; x2,y2 in */
            }
          }
          else {
            mode = 3; /* x1,y1 in;  x2,y2 out */
          }
        }

        while (!done) {

          switch (mode) {

          case 1: /* x1,y1 out, x2,y2 in... clip x1,y1*/
          case 2: /* both out, clip x1,y1 first */

            if (code1 & mask) { /* clip at left edge */
              y1 = y1 + (y2 - y1) * (cur_state->clipmin.x - x1) / (x2 - x1);
              x1 = cur_state->clipmin.x;
            }
            else if ((code1 >> 1) & mask) { /* clip at right edge */
              y1 = y1 + (y2 - y1) * (cur_state->clipmax.x - x1) / (x2 - x1);
              x1 = cur_state->clipmax.x;
            }

            if ((code1 >> 2) & mask) { /* clip at bottom edge */
              x1 = x1 + (x2 - x1) * (cur_state->clipmin.y - y1) / (y2 - y1);
              y1 = cur_state->clipmin.y;
            }
            else if ((code1 >> 3) & mask) { /* clip at top edge */
              x1 = x1 + (x2 - x1) * (cur_state->clipmax.y - y1) / (y2 - y1);
              y1 = cur_state->clipmax.y;
            }

            switch (mode) {
            case 1: /* not looping through, draw it */
              vdmova(&x1, &y1);
              vdlina(&x2, &y2);
              done = TRUE;
              break;

            case 2: /* both were out, encode new point */
              code1 = 0;
              if ((cur_state->clipmax.y - y1) < 0) {
                code1 = code1 | (1 << 3);
              }
              else if ((y1 - cur_state->clipmin.y) < 0) {
                code1 = code1 | (1 << 2);
              }

              if ((cur_state->clipmax.x - x1) < 0) {
                code1 = code1 | (1 << 1);
              }
              else if ((x1 - cur_state->clipmin.x) < 0) {
                code1 = code1 | (1 << 0);
              }

              if (code1) {
                done = TRUE; /* new point is out, reject */
              }
              else {
                mode = 3; /* x1,y1 in; x2,y2 out; clip x2,y2 */
              }
              break;

            } /* end inner switch */

            break;

          case 3: /* x1,y1 in, x2,y2 out...clip */

            if (code2 & mask) { /* clip at left edge */
              y2 = y2 + (y1 - y2) * (cur_state->clipmin.x - x2) / (x1 - x2);
              x2 = cur_state->clipmin.x;
            }
            else if ((code2 >> 1) & mask) { /* clip at right edge */
              y2 = y2 + (y1 - y2) * (cur_state->clipmax.x - x2) / (x1 - x2);
              x2 = cur_state->clipmax.x;
            }

            if ((code2 >> 2) & mask) { /* clip at bottom edge */
              x2 = x2 + (x1 - x2) * (cur_state->clipmin.y - y2) / (y1 - y2);
              y2 = cur_state->clipmin.y;
            }

            else if ((code2 >> 3) & mask) { /* clip at top edge */
              x2 = x2 + (x1 - x2) * (cur_state->clipmax.y - y2) / (y1 - y2);
              y2 = cur_state->clipmax.y;
            }

            vdmova(&x1, &y1);
            vdlina(&x2, &y2);
            done = TRUE;
            break;

          default: done = TRUE; break;

          } /* end switch */
        } /* end while !done */
      } /* end for j */
    } /* end else clipping is on */

    /* flag that the page has been marked */
    cur_state->pic_dirty = CDIRTY;

  } /* end for i */
} /* end xcdjpl */

/* POLYMARKER */
static void xcpm(anything **params, int num_surfaces, anything **surf_list)
{
  int    i, j;         /* indices for loops */
  int    np;           /* number of points */
  float *x, *y;        /* x,y values */
  float  x_vdi, y_vdi; /* x,y values in NDC */
  int    ok;           /* flag TRUE =mark lies in clip region */

  for (i = 0; i < num_surfaces; ++i) {

    np = *(int *)params[1];
    x  = (float *)params[2];
    y  = (float *)params[3];

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* check for errors */
    if (MAX_POLYMARKER != -1 && *(int *)params[1] > MAX_POLYMARKER) {
      /* error 6:301 - Number of points too large. Max is used */
      report_error(cur_state, 6, 301, *(short *)params[0]);
      np = MAX_POLYMARKER;
    }
    else {
      np = *(int *)params[1];
    }

    /* set SVDI foreground color if needed */
    set_foreground_color(cur_state, cur_state->mark_color);

    /* set marker type */
    /* --only support "dot" for now - so nothing to set */

    /* set marker size */
    /* --only support "dot" for now - so nothing to set */

    for (j = 0; j < np; j++) {

      /* map VDC to NDC */
      x_vdi = x[j] * cur_state->xscale + cur_state->xoffset;
      y_vdi = y[j] * cur_state->yscale + cur_state->yoffset;

      ok = TRUE;
      if (cur_state->clip_on) { /* clip is on */
        /* ...clip the whole polymarker - it's in or it's out */
        ok = ((x_vdi <= cur_state->clipmax.x && x_vdi >= cur_state->clipmin.x) &&
              (y_vdi <= cur_state->clipmax.y && y_vdi >= cur_state->clipmin.y));
      }

      if (ok) {
        vdpnta(&x_vdi, &y_vdi);
      }

    } /* end for j */

    /* flag that the page has been marked */
    cur_state->pic_dirty = CDIRTY;

  } /* end for i */
} /* end xcpm */

/* TEXT */
static void xctx(anything **params, int num_surfaces, anything **surf_list)
{
  int   i, j;                 /* indices for loops */
  float x, y;                 /* x,y starting position */
  char *ch;                   /* input character string */
  char  c;                    /* 1 character */
  int   np;                   /* number of chars in text string */
  int   ich[MAX_TEXT_LENGTH]; /* array of ascii chars to pass SVDI */
  int   ok;                   /* flag TRUE =line lies in clip region */
  float temp_array[14];       /* used for SVDI inquiries */
  float char_height;          /* char height, in VDC or NDC */
  float char_width;           /* char width, in VDC or NDC */
  int   skip;                 /* nbr of chars to skip when clipping */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* error checking */
    if (*(int *)params[3] != CFINAL) {
      /* compound text not supported */
      /* error 2:001 - Request for unsupported feature. Function ignored */
      report_error(cur_state, 2, 001, *(short *)params[0]);
      break;
    }
    /* params[5] is an added parameter which is length of string */
    else if (MAX_TEXT_LENGTH != -1 && *(int *)params[5] > MAX_TEXT_LENGTH) {
      /* error 6:301 - Number of points too large. Max is used */
      report_error(cur_state, 6, 301, *(short *)params[0]);
      np = MAX_TEXT_LENGTH;
    }
    else {
      np = *(int *)params[5];
    }

    /* if there is nothing to draw, break  */
    if (np <= 0) {
      break;
    }

    /* set SVDI foreground color if needed */
    set_foreground_color(cur_state, cur_state->text_color);

    /* find out what the SVDI character height/width is */
    vdiqos(&temp_array[1]);

    /* convert stuff to NDC */
    char_width  = temp_array[7];
    char_height = temp_array[6];
    /* map VDC to NDC */
    x = *(float *)params[1] * cur_state->xscale + cur_state->xoffset;
    y = *(float *)params[2] * cur_state->yscale + cur_state->yoffset;

    ok   = TRUE;
    skip = 0;

    /* this text clipping assumes no backspaces, etc.,and only the
     * default character orientation.
     */
    if (cur_state->clip_on) { /* clipping is on, do it */

      /* if the string doesn't fit in y, reject it */
      ok = (y >= cur_state->clipmin.y && (y + char_height) <= cur_state->clipmax.y);

      if (ok) { /* ok in y, check it in x */

        /* check for trivial reject */
        ok = (x < cur_state->clipmax.x && (x + (np * char_width)) > cur_state->clipmin.x);

        if (ok) { /* it's in there, clip it */

          /* clip to left edge */
          if (x < cur_state->clipmin.x) {

            /* ...compute how many chars to skip */
            skip = ceil((cur_state->clipmin.x - x) / char_width);
            /* ...compute new x */
            x = x + ((float)skip * char_width);
            /* ...compute new number of chars */
            np = np - skip;
          }

          /* clip to right edge */
          if (x + (np * char_width) > cur_state->clipmax.x) {

            /* ...compute how many chars fit */
            np = (int)((cur_state->clipmax.x - x) / char_width);
          }
        } /* end if ok in x */
      } /* end if ok in y */
    } /* end if clip_on */

    ok = ok && np > 0; /* make sure there is still some text left */

    if (ok) { /* everything is still ok */

      /* move to position, and draw the text */
      vdmova(&x, &y);

      ch = &(((char *)params[4])[skip]);
      for (j = 0; j < np; j++) {
        c = *ch++;
        /* this works for ASCII machines */
        ich[j] = (int)c;
      }
      vdtext(&j, ich);

    } /* end if ok */

  } /* end for each surface */
} /* end xctx */

/* POLYGON */
static void xcpg(anything **params, int num_surfaces, anything **surf_list)
{
  int          i, j;              /* indices for loops */
  int          np;                /* number of points */
  float       *x, *y;             /* x,y VDC values */
  float        xnew[MAX_POLYGON]; /* clipped x,y values (VDC or NDC) */
  float        ynew[MAX_POLYGON];
  int          npnew;              /* number of points after clip */
  float        temp_array[14];     /* temp used for SVDI inquire */
  int          qdc_index;          /* for inquire SVDI */
  int          vdi_ls;             /* SVDI linestyle */
  int          ok;                 /* flag TRUE =line lies in clip region */
  static float poly_support = -1.; /* save whether SVDI does polygons */
  static float vdi_polymax  = -1.; /* store SVDI polygon max value */

  /* find out SVDI support for polygons. */
  if (poly_support == -1.0f) {
    qdc_index = 24;
    vdiqdc(&qdc_index, &poly_support);
  }

  /* get and store max polygon points from SVDI */
  if (vdi_polymax == -1.0f) {
    qdc_index = 25;
    vdiqdc(&qdc_index, &vdi_polymax);
  }

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* error checking */
    if (MAX_POLYGON != -1 && (*(int *)params[1] > MAX_POLYGON || *(int *)params[1] > vdi_polymax)) {
      /* error 6:302 - Number of points in polygon too large. Max used */
      report_error(cur_state, 6, 302, *(short *)params[0]);
      np = (MAX_POLYGON < vdi_polymax) ? MAX_POLYGON : vdi_polymax;
    }
    else {
      np = *(int *)params[1];
    }

    x = (float *)params[2];
    y = (float *)params[3];

    /* set SVDI foreground color if needed */
    set_foreground_color(cur_state, cur_state->fill_color);

    /* map VDC to NDC */
    for (j = 0; j < np; j++) {

      xnew[j] = *x++ * cur_state->xscale + cur_state->xoffset;
      ynew[j] = *y++ * cur_state->yscale + cur_state->yoffset;

    } /* end for j */

    ok = TRUE;
    if (cur_state->clip_on) { /* clip is on */

      /* poly clip returns new x,y values and a count */
      /* this routine should be done inline for efficiency */
      ok = poly_clip(&cur_state->clipmin, &cur_state->clipmax, xnew, ynew, np, xnew, ynew, &npnew);
    }
    else { /* no clpping */
      npnew = np;
    }

    if (ok) { /* everything's ok, do the polygon */

      /* if hollow or if no polygon support, draw the border lines.
         if solid, call vdpoly */
      if ((cur_state->interior_style == CHOLLO) || /* hallow OR   */
          (poly_support < 3.0f)) {                 /*no poly supp */

        vdiqos(&temp_array[1]); /* make sure line style is solid */
        if (temp_array[4] != 0) {
          vdi_ls = 0;
          vdstls(&vdi_ls);
        }

        vdmova(&xnew[0], &ynew[0]);
        for (j = 1; j < npnew; j++) {
          vdlina(&xnew[j], &ynew[j]);
        }
        vdlina(&xnew[0], &ynew[0]);
        if (temp_array[4] != 0) { /* set line style back */
          vdi_ls = (int)temp_array[4];
          vdstls(&vdi_ls);
        }
      } /* end if hollow OR no poly support */
      else { /* solid polygon */
        vdpoly(xnew, ynew, &npnew);
      }

      /* flag that the page has been marked */
      cur_state->pic_dirty = CDIRTY;

    } /*end if ok */

  } /* end for i */
} /* end xcpg */

/* CELL ARRAY */
static void xcca(anything **params, int num_surfaces, anything **surf_list)
{
  int   i;                       /* index for loop on surfaces */
  int   j, k;                    /* loop indices */
  int   index, index_inc, count; /* array indices */
  int   nx, ny;                  /* number of cells in x,y */
  int   nx1, ny1;                /* number of cells in x,y after clipping */
  int  *cells;                   /* color values */
  int   ix, iy, yinc;            /* SVDI logical raster coordinates */
  float x1, x2, y1, y2;          /* corners of rectangle in NDC */
  float xmin, xmax, ymin, ymax;  /* SVDI raster viewport (NDC) */
  int   xdir, ydir;              /* x,y directions on the device */
  int   imap;
  float xcell, ycell; /* area per cell */
  int   ok;           /* flag TRUE =line lies in clip region */
  int   skipx, skipy; /* number of cells clipped in xmin,ymin */

  /* if there isn't anything to draw, return */
  if (*(int *)params[7] == 0 || *(int *)params[8] == 0) {
    return;
  }

  /* if cell array has negative dimension, use absolute values */
  nx = nx1 = abs(*(int *)params[7]);
  ny = ny1 = abs(*(int *)params[8]);

  cells = (int *)params[10];

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* if color selection mode is in error, ignore function */
    if (cur_state->csm == -1) {
      break;
    }

    /* for now, only allow 4 cases (same as pixel array),
     * where nx runs only horizontal, and ny runs vertical.
     * this routine needs to be fixed to allow all 8 cases
     * Also, this routine does not allow skewing.
     */
    if (*(float *)params[2] != *(float *)params[6] || *(float *)params[3] != *(float *)params[5]) {
      /* error 1:-106 Illegal cell array. Function ignored */
      report_error(cur_state, 1, -106, *(short *)params[0]);
      return;
    }

    /* map corners to NDC */
    x1 = *(float *)params[1] * cur_state->xscale + cur_state->xoffset;
    x2 = *(float *)params[3] * cur_state->xscale + cur_state->xoffset;

    y1 = *(float *)params[2] * cur_state->yscale + cur_state->yoffset;
    y2 = *(float *)params[4] * cur_state->yscale + cur_state->yoffset;

    /* figure out the raster viewport limits, and figure
     * out x and y directions based on mappings and
     * rectangle orientation.
     */
    if (x1 < x2) {
      xmin = x1;
      xmax = x2;
      xdir = CINCR;
    }
    else {
      xmin = x2;
      xmax = x1;
      xdir = CDECR;
    }

    if (y1 < y2) {
      ymin = y1;
      ymax = y2;
      ydir = CINCR;
    }
    else {
      ymin = y2;
      ymax = y1;
      ydir = CDECR;
    }

    ok    = TRUE;
    skipx = skipy = 0;

    if (cur_state->clip_on) { /* clipping is on, do it */

      /* trivial reject */
      ok = ((xmin < cur_state->clipmax.x && xmax > cur_state->clipmin.x) &&
            (ymin < cur_state->clipmax.y && ymax > cur_state->clipmin.y));

      if (ok) { /* it's in there, clip it */

        /* clip in x first */
        xcell = (xmax - xmin) / nx; /* area per cell */

        /* clip at left edge */
        if (xmin < cur_state->clipmin.x) {

          /* ...compute how many cells to skip */
          skipx = ceil((cur_state->clipmin.x - xmin) / xcell);
          /* ...compute new xmin */
          xmin = xmin + ((float)skipx * xcell);
          /* ...compute new nx */
          nx1 = nx - skipx;
        } /* end clip at left edge */

        /* clip at right edge */
        if (xmax > cur_state->clipmax.x) {

          /* ...compute new nx - how much will fit */
          nx1 = (int)((cur_state->clipmax.x - xmin) / xcell);
          /* ...compute new xmax */
          xmax = xmin + (nx1 * xcell);
        } /* end clip at right edge */

        /* now clip in y */
        ycell = (ymax - ymin) / ny; /* area per cell */

        /* clip at bottom edge */
        if (ymin < cur_state->clipmin.y) {

          /* ...compute how many cells to skip */
          skipy = ceil((cur_state->clipmin.y - ymin) / ycell);
          /* ...compute new ymin */
          ymin = ymin + ((float)skipy * ycell);
          /* ...compute new ny */
          ny1 = ny - skipy;
        } /* end clip at bottom edge */

        /* clip at top edge */
        if (ymax > cur_state->clipmax.y) {

          /* ...compute new ny - how much will fit */
          ny1 = (int)((cur_state->clipmax.y - ymin) / ycell);
          /* ...compute new ymax */
          ymax = ymin + (ny1 * ycell);
        } /* end clip at top edge */

      } /* end if ok */
    } /* end if clip is on */

    /* make sure there is still something left */
    ok = ok && (nx1 > 0 && ny1 > 0);

    if (ok) { /* everything is still ok */

      /* set the raster space */
      vdstrs(&nx1, &ny1);

      imap = 3;
      vbstmp(&imap);
      imap = 5;
      vbstmp(&imap);

      /* set the raster viewport */
      vdstrv(&xmin, &xmax, &ymin, &ymax);

      /* direct color or indexed */
      if (cur_state->csm == CINDEX) { /* indexed color */

        /* Indexing into the cell array is based on clipping info
         * and the orientation of the rectangle. SVDI+raster puts
         * [0,0] at the top left corner.
         */

        /* There are 4 different cases based on the orientation
         * of the rectangle. ( really, there are 8, but I'm not
         * supporting the other 4 right now ).  If x is increasing,
         * the most natural case, the cell array will be written
         * out one raster line at a time, either top to bottom or
         * bottom to top depending on the y direction. If x is
         * decreasing, the cell array will have to be copied and
         * reordered, so it will be written out in chunks.  In this
         * case, it will always be written out top to bottom.
         */

        if (xdir == CINCR) { /* x increasing */

          ix = 0;
          /* compute iy and y increment based on y direction */
          if (ydir == CINCR) { /* y increasing */
            iy   = ny1 - 1;
            yinc = -1;
          }
          else {
            iy   = 0;
            yinc = 1;
          }

          /* This special case is for efficiency, since is will probably
           * be most of the cases: If no clipping is on and [0,0] is at
           * the top left (like SVDI+raster) draw the whole cell array.
           */
          if (!cur_state->clip_on && iy == 0) {
            count = nx * ny;
            vdpixi(&ix, &iy, &cells[0], &count);
          } /* end special case */

          else { /* special case doesn't apply, draw one line at a time */

            /* starting index based on clipping info */
            index = (ny - ny1 - skipy) * nx + skipx;

            /* output one raster line at a time */
            for (k = 0; k < ny1; k++) {
              vdpixi(&ix, &iy, &cells[index], &nx1);
              index = index + nx; /* compute next index */
              iy    = iy + yinc;  /* compute next raster line */
            } /* end for k */
          } /* end else no special case */
        } /* end x increasing */

        else { /* x decreasing */

          /* always write top to bottom */
          ix = iy = count = 0;

          /* starting index and index increment are based on
             clipping info and the y direction*/

          if (ydir == CINCR) { /* y increasing */
            index     = (ny1 + skipy) * nx - skipx - 1;
            index_inc = -nx + nx1;
            yinc      = -1;
          }
          else { /* y decreasing */
            index     = (ny - ny1 - skipy) * nx + (nx - 1) - skipx;
            index_inc = nx + nx1;
            yinc      = 1;
          }

          /* reorder the cell colors */
          for (k = 0; k < ny1; k++) {
            for (j = 0; j < nx1; j++) {
              assert(count < MAX_ARRAY);
              varray[count++] = cells[index--];
            }

            /* if another row won't fit, or done, output varray */
            if (count + nx1 > MAX_ARRAY || k == ny1 - 1) {
              vdpixi(&ix, &iy, varray, &count);
              count = 0;                     /* reset varray counter */
              iy    = iy + (yinc * (k + 1)); /* update iy */
            }

            /* compute next index */
            index = index + index_inc;

          } /* end for k */
        } /* end x decreasing */
      } /* end if indexed color */

      else { /*  direct color */

        /* Indexing into the cells array is based on clipping info,
         * and the orientation of the cell array.  In SVDI+raster,
         * the top left corner is [0,0]. For direct color, the cell
         * array is always drawn from top to bottom. If the cell array
         * is oriented with [0,0] at the the bottom left corner, the
         * way the cells array is indexed takes care of the flipping.
         */

        ix = iy = count = 0;

        if (xdir == CINCR) { /* x increasing */

          /* starting index and index increment are based on
             clipping info and the y direction */
          if (ydir == CINCR) { /* y increasing */
            index     = ((ny1 + skipy - 1) * nx * 3) + (skipx * 3);
            index_inc = -(nx + nx1) * 3;
          }
          else { /* y decreasing */
            index     = (ny - ny1 - skipy) * nx * 3 + (skipx * 3);
            index_inc = (nx - nx1) * 3;
          }

          /* store as much info as possible before calling vdpixl */
          for (k = 0; k < ny1; k++) {
            for (j = 0; j < nx1; j++) {
              assert(count < MAX_ARRAY);
              rarray[count] = (float)cells[index++] / 255.0f;
              garray[count] = (float)cells[index++] / 255.0f;
              barray[count] = (float)cells[index++] / 255.0f;
              count++;
            } /* end for j */

            /* if another row won't fit, or done, output rgb arrays */
            if (count + nx1 > MAX_ARRAY || k == ny1 - 1) {
              vdpixl(&ix, &iy, rarray, garray, barray, &count);
              count = 0;     /* reset rgb array counter */
              iy    = k + 1; /* update iy */
            }

            /* compute next index into the cells array */
            index = index + index_inc;

          } /* end for k */
        } /* end if x increasing */

        else { /* x decreasing */

          /* starting index and index increment are based on
             clipping info and the y direction */
          if (ydir == CINCR) { /* y increasing */
            index     = (ny1 + skipy) * nx * 3 - (skipx * 3) - 3;
            index_inc = (-nx + nx1) * 3;
          }
          else { /* y decreasing */
            index     = (ny - ny1 - skipy) * nx * 3 + ((nx - 1) * 3) - (skipx * 3);
            index_inc = (nx + nx1) * 3;
          }

          /* store as much info as possible before calling vdpixl */
          for (k = 0; k < ny1; k++) {
            for (j = 0; j < nx1; j++) {
              assert(count < MAX_ARRAY);
              rarray[count] = (float)cells[index] / 255.0f;
              garray[count] = (float)cells[index + 1] / 255.0f;
              barray[count] = (float)cells[index + 2] / 255.0f;
              index         = index - 3;
              count++;
            } /* end for j */

            /* if another row won't fit, or done, output rgb arrays */
            if (count + nx1 > MAX_ARRAY || k == ny1 - 1) {
              vdpixl(&ix, &iy, rarray, garray, barray, &count);
              count = 0;     /* reset rgb array counter */
              iy    = k + 1; /* update iy */
            }

            /* compute next index into the cells array */
            index = index + index_inc;

          } /* end for k */
        } /* end else x decreasing */
      } /* end else direct color */

      /* flag that the page has been marked */
      cur_state->pic_dirty = CDIRTY;

    } /* end if ok */
  } /* end for each surface */
} /* end xcca */

/* PIXEL ARRAY */
static void xcpxa(anything **params, int num_surfaces, anything **surf_list)
{
  int          i;                      /* index for loop on surfaces */
  int          j, k;                   /* loop indices */
  int         *pxclrs;                 /* pixel color values */
  int          nx, ny;                 /* size of pixel color array */
  int          nx1, ny1;               /* number of pixels after clipping */
  int          index, index_inc;       /* index into pxclrs array */
  int          count;                  /* index into temp pixel arrays */
  int          repx, repy;             /* replication in x and y */
  int          xdir, ydir;             /* x and y direction on the DEVICE */
  float        x_orig, y_orig;         /* x and y origin */
  float        x_dist, y_dist;         /* distance pixels run from the origin */
  int          ix, iy, yinc;           /* SVDI logical raster coordinates */
  float        xmin, xmax, ymin, ymax; /* SVDI raster viewport after clipping */
  int          ok;                     /* flag TRUE =line lies in clip region */
  int          skipx, skipy;           /* count pixels skipped in xmin,ymin */
  int          imap;
  int          qrs_index;     /* index for SVDI inquiries */
  float        temp_array[2]; /* temporary array for SVDI inquiries */
  float        x_pxl_scale;   /* for converting pixels to NDC/VDC space */
  float        y_pxl_scale;
  static float x_pixels = -1.; /* number of pixels in full NDC space */
  static float y_pixels = -1.;
  float        zero     = 0.0;

  /* if nx,ny,repx,or repy are non positive, nothing is drawn */
  if (*(int *)params[3] < 0 || *(int *)params[4] < 0 || *(int *)params[5] < 0 ||
      *(int *)params[6] < 0) {
    return;
  }

  /* set up initial values */
  nx = nx1 = *(int *)params[3];
  ny = ny1 = *(int *)params[4];
  repx     = *(int *)params[5];
  repy     = *(int *)params[6];
  pxclrs   = (int *)params[9];

  /* find (and save) the number of pixels in full NDC space */
  if (x_pixels == -1.0f && y_pixels == -1.0f) {
    vdstrv(&zero, &dev_descrip.xndc_max, &zero, &dev_descrip.yndc_max);
    qrs_index = 2;
    vdiqrs(&qrs_index, temp_array);
    x_pixels = temp_array[0];
    y_pixels = temp_array[1];
  } /* end find number of pixels */

  /* loop through surfaces */
  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* if color selection mode is in error, ignore function */
    if (cur_state->csm == -1) {
      break;
    }

    /* check enumerated types */
    if ((*(int *)params[7] != CINCR && *(int *)params[7] != CDECR) ||
        (*(int *)params[8] != CINCR && *(int *)params[8] != CDECR)) {
      /* error 1:-101 enumerated parameter out of range. Function ignored */
      report_error(cur_state, 1, -101, *(short *)params[0]);
      break;
    }

    /* figure out what x direction and y direction are in terms
     * of the device, since x direction and y direction are dependent
     * on the VDC to device mapping. If xscale is positive, then
     * VDC and VP are either both increasing or both decreasing and
     * x direction on the device is the same as x direction in VDC.
     * If xscale is negative, one is increasing and one is decreasing,
     * and x direction on the device is opposite of x direction in
     * VDC. Same holds for y (where "increasing" means bottom to top)
     */
    if (cur_state->xscale > 0) {
      xdir = *(int *)params[7];
    }
    else {
      xdir = (*(int *)params[7] == CINCR) ? CDECR : CINCR;
    }

    if (cur_state->yscale > 0) {
      ydir = *(int *)params[8];
    }
    else {
      ydir = (*(int *)params[8] == CINCR) ? CDECR : CINCR;
    }

    /* map to NDC */

    /* origin */
    x_orig = *(float *)params[1] * cur_state->xscale + cur_state->xoffset;
    y_orig = *(float *)params[2] * cur_state->yscale + cur_state->yoffset;

    /* scale to convert pixels to NDC space */
    x_pxl_scale = dev_descrip.xndc_max / x_pixels;
    y_pxl_scale = dev_descrip.yndc_max / y_pixels;

    /* figure out raster viewport mins and maxes */
    x_dist = nx * repx * x_pxl_scale;
    y_dist = ny * repy * y_pxl_scale;

    /* ...x is increasing or decreasing */
    if (xdir == CINCR) {
      xmin = x_orig;
      xmax = x_orig + x_dist;
    }
    else {
      xmin = x_orig - x_dist;
      xmax = x_orig;
    }
    /* y is increasing or decreasing */
    if (ydir == CINCR) {
      ymin = y_orig;
      ymax = y_orig + y_dist;
    }
    else {
      ymin = y_orig - y_dist;
      ymax = y_orig;
    }

    ok    = TRUE;
    skipx = skipy = 0;

    if (cur_state->clip_on) { /* do the clipping */

      /* trivial reject */
      ok = ((xmin < cur_state->clipmax.x && xmax > cur_state->clipmin.x) &&
            (ymin < cur_state->clipmax.y && ymax > cur_state->clipmin.y));

      if (ok) { /* it's in there, clip it */

        /* clip in x first */
        if (xmin < cur_state->clipmin.x) { /* ...left edge */

          /* ...compute how many pixels to skip */
          skipx = ceil((cur_state->clipmin.x - xmin) / (x_pxl_scale * repx));
          /* ...compute new xmin */
          xmin = xmin + ((float)skipx * x_pxl_scale * repx);
          /* ...compute new nx */
          nx1 = nx - skipx;
        } /* end clip left edge */

        if (xmax > cur_state->clipmax.x) { /* ...right edge */

          /* ...compute new nx */
          nx1 = (int)((cur_state->clipmax.x - xmin) / (x_pxl_scale * repx));
          /* ...compute new xmax */
          xmax = xmin + (nx1 * x_pxl_scale * repx);
        } /* end clip right edge */

        /* now clip in y */
        if (ymin < cur_state->clipmin.y) { /* ...bottom edge */

          /* ...compute how many pixels to skip */
          skipy = ceil((cur_state->clipmin.y - ymin) / (y_pxl_scale * repy));
          /* ...compute new ymin */
          ymin = ymin + ((float)skipy * y_pxl_scale * repy);
          /* ...compute new ny */
          ny1 = ny - skipy;
        } /* end bottom left edge */

        if (ymax > cur_state->clipmax.y) { /* ...top edge */

          /* ...compute new ny */
          ny1 = (int)((cur_state->clipmax.y - ymin) / (y_pxl_scale * repy));
          /* ...comput new ymax */
          ymax = ymin + (ny1 * y_pxl_scale * repy);
        } /* end clip top edge */

      } /* end if ok */
    } /* end if clip_on */

    /* make sure there is still something to draw */
    ok = ok && (nx1 > 0 && ny1 > 0);

    if (ok) { /* everything is still ok */

      /* set the raster viewport */
      vdstrv(&xmin, &xmax, &ymin, &ymax);

      imap = 2;
      vbstmp(&imap);
      imap = 5;
      vbstmp(&imap);

      /* set the raster space to number of device pixels in viewport */
      vdstrs(&nx1, &ny1);

      /* direct color or indexed */
      if (cur_state->csm == CINDEX) { /* indexed color */

        /* Indexing into the pixel array is based on clipping info
         * and the orientation of the pixel array, which is
         * determined by whether x and y are increasing or
         * decreasing in VDC from the origin.  SVDI+raster puts
         * [0,0] at the top left corner.
         */

        /* There are 4 different cases based on the x and y
         * directions.  If x is increasing, the most natural
         * case, the pixel array will be written out one raster
         * line at a time, either top to bottom or bottom to
         * top depending on the y direction. If x is decreasing,
         * the pixel array will have to be copied and reordered,
         * so it will be written out in chunks.  In this case,
         * it will always be written out top to bottom.
         */

        if (xdir == CINCR) { /* x increasing */

          ix = 0;
          /* compute iy and y increment based on y direction */
          if (ydir == CINCR) { /* y increasing */
            iy   = ny1 - 1;
            yinc = -1;
          }
          else { /* y decreasing */
            iy   = 0;
            yinc = 1;
          }

          /* This special case is for efficiency, since is will probably
           * be most of the cases: If no clipping is on and [0,0] is at
           * the top left (like SVDI+raster) draw the whole pixel array.
           */
          if (!cur_state->clip_on && iy == 0) {
            count = nx * ny;
            vdpixi(&ix, &iy, &pxclrs[0], &count);
          } /* end special case */

          else { /* special case doesn't apply, draw one line at a time */

            /* starting index based on clipping info */
            index = (ny - ny1 - skipy) * nx + skipx;

            /* output one raster line at a time */
            for (k = 0; k < ny1; k++) {
              vdpixi(&ix, &iy, &pxclrs[index], &nx1);
              index = index + nx; /* compute next index */
              iy    = iy + yinc;  /* compute next raster line */
            } /* end for k */
          } /* end else no special case */
        } /* end x increasing */

        else { /* x decreasing */

          /* always write top to bottom */
          ix = iy = count = 0;

          /* starting index and index increment are based on
             clipping info and the y direction*/

          if (ydir == CINCR) { /* y increasing */
            index     = (ny1 + skipy) * nx - skipx - 1;
            index_inc = -nx + nx1;
            yinc      = -1;
          }
          else { /* y decreasing */
            index     = (ny - ny1 - skipy) * nx + (nx - 1) - skipx;
            index_inc = nx + nx1;
            yinc      = 1;
          }

          /* reorder the pixel colors */
          for (k = 0; k < ny1; k++) {
            for (j = 0; j < nx1; j++) {
              assert(count < MAX_ARRAY);
              varray[count++] = pxclrs[index--];
            }

            /* if another row won't fit, or done, output varray */
            if (count + nx1 > MAX_ARRAY || k == ny1 - 1) {
              vdpixi(&ix, &iy, varray, &count);
              count = 0;                     /* reset varray counter */
              iy    = iy + (yinc * (k + 1)); /* update iy */
            }

            /* compute next index */
            index = index + index_inc;

          } /* end for k */
        } /* end x decreasing */
      } /* end if indexed color */

      else { /*  direct color */

        /* For direct color, the color values always have to be
         * copied, so the raster lines are written out in chucks.
         * The raster lines will always be written top to bottom;
         * indexing and reordering of the color values will be
         * based on x and y direction.
         */

        ix = iy = count = 0;

        if (xdir == CINCR) { /* x increasing */

          /* starting index and index increment are based on
             clipping info and the y direction */
          if (ydir == CINCR) { /* y increasing */
            index     = ((ny1 + skipy - 1) * nx * 3) + (skipx * 3);
            index_inc = -(nx + nx1) * 3;
          }
          else { /* y decreasing */
            index     = (ny - ny1 - skipy) * nx * 3 + (skipx * 3);
            index_inc = (nx - nx1) * 3;
          }

          /* store as much info as possible before calling vdpixl */
          for (k = 0; k < ny1; k++) {
            for (j = 0; j < nx1; j++) {
              assert(count < MAX_ARRAY);
              rarray[count] = (float)pxclrs[index++] / 255.0f;
              garray[count] = (float)pxclrs[index++] / 255.0f;
              barray[count] = (float)pxclrs[index++] / 255.0f;
              count++;
            } /* end for j */

            /* if another row won't fit, or done, output rgb arrays */
            if (count + nx1 > MAX_ARRAY || k == ny1 - 1) {
              vdpixl(&ix, &iy, rarray, garray, barray, &count);
              count = 0;     /* reset rgb array counter */
              iy    = k + 1; /* update iy */
            }

            /* compute next index into the pxclrs array */
            index = index + index_inc;

          } /* end for k */
        } /* end if x increasing */

        else { /* x decreasing */

          /* starting index and index increment are based on
             clipping info and the y direction */
          if (ydir == CINCR) { /* y increasing */
            index     = (ny1 + skipy) * nx * 3 - (skipx * 3) - 3;
            index_inc = (-nx + nx1) * 3;
          }
          else { /* y decreasing */
            index     = (ny - ny1 - skipy) * nx * 3 + ((nx - 1) * 3) - (skipx * 3);
            index_inc = (nx + nx1) * 3;
          }

          /* store as much info as possible before calling vdpixl */
          for (k = 0; k < ny1; k++) {
            for (j = 0; j < nx1; j++) {
              assert(count < MAX_ARRAY);
              rarray[count] = (float)pxclrs[index] / 255.0f;
              garray[count] = (float)pxclrs[index + 1] / 255.0f;
              barray[count] = (float)pxclrs[index + 2] / 255.0f;
              index         = index - 3;
              count++;
            } /* end for j */

            /* if another row won't fit, or done, output rgb arrays */
            if (count + nx1 > MAX_ARRAY || k == ny1 - 1) {
              vdpixl(&ix, &iy, rarray, garray, barray, &count);
              count = 0;     /* reset rgb array counter */
              iy    = k + 1; /* update iy */
            }

            /* compute next index into the pxclrs array */
            index = index + index_inc;

          } /* end for k */
        } /* end else x decreasing */
      } /* end else direct color */

      /* flag that the page has been marked */
      cur_state->pic_dirty = CDIRTY;

    } /* end if ok */
  } /* end for each surface */
} /* end xcpxa */

/* LINE TYPE */
static void xclnt(anything **params, int num_surfaces, anything **surf_list)
{

  int i;         /* index for loop on surfaces */
  int vdi_style; /* SVDI linestyle */

  /* CGI line types: */
  /* 1: solid - supported
   * 2: dash - supported
   * 3: dot - supported
   * 4: dash dot - supported
   * 5: dash dot dot - not supported
   */

  /* SVDI linestyles: */
  /* 0: solid
   * 1: dotted
   * 2: dot dash
   * 3: short dash
   * 4: long dash
   * 5: medium dash
   */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* let SVDI catch any out of range indices */
    cur_state->line_type = *(int *)params[1];

    /* make sure it really needs to be set */
    if (cur_state->line_type != cur_state->vdi_attrib.line_type) {

      cur_state->vdi_attrib.line_type = cur_state->line_type;

      /* map to SVDI styles */
      switch (cur_state->line_type) {
      case 1: vdi_style = 0; break;
      case 2: vdi_style = 5; break;
      case 3: vdi_style = 1; break;
      case 4: vdi_style = 2; break;
      case 5: vdi_style = 0; break; /* not supported, set to solid */
      default: vdi_style = cur_state->line_type; break;
      } /* end switch */

      vdstls(&vdi_style);
    } /* end it needs to be set */

  } /* end for each surface */
} /* end xclnt */

/* LINE WIDTH */
static void xclnw(anything **params, int num_surfaces, anything **surf_list)
{

  int   i;              /* index for loop on surfaces */
  float vdi_lw;         /* SVDI linewidth */
  float temp_array[14]; /* use for SVDI inquiry */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* make sure it needs to be updated */
    if (*(float *)params[1] != cur_state->vdi_attrib.line_width) {

      /* save the original in vdi_attrib array */
      cur_state->vdi_attrib.line_width = *(float *)params[1];

      /* convert to NDC and set the linewidth */
      vdi_lw = dev_descrip.linewidth_nominal * *(float *)params[1] * fabs(cur_state->xscale);
      vdstlw(&vdi_lw);

      /* find out what it really got set at */
      vdiqos(&temp_array[1]);

      /* map it back to (scaled) VDC and store it in state table */
      cur_state->line_width =
          (temp_array[5] * fabs(cur_state->xscale) * 100) / dev_descrip.linewidth_nominal;

    } /* end it needs to be updated */

  } /* end for each surface */
} /* end xclnw */

/* LINE COLOR */
static void xclnc(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* if color selection mode is in error, ignore function */
    if (cur_state->csm == -1) {
      break;
    }

    if (cur_state->csm == CINDEX) { /* indexed mode */

      /* let SVDI catch out of range indices */
      cur_state->line_csm      = CINDEX;
      cur_state->line_color[0] = *(int *)params[1];
    } /* end indexed mode */

    else { /* direct color mode */

      /* let SVDI catch illegal color values */
      cur_state->line_csm      = CDRECT;
      cur_state->line_color[0] = *((int *)params[1] + 0);
      cur_state->line_color[1] = *((int *)params[1] + 1);
      cur_state->line_color[2] = *((int *)params[1] + 2);
    } /* end direct color mode */

  } /* end for each surface */
} /* end xclnc */

/* MARKER TYPE */
static void xcmkt(int num_surfaces, anything **surf_list)
{
  /* CGI marker types: */
  /* 1: dot - supported
   * 2: plus - not supported
   * 3: asterisk - not supported
   * 4: circle - not supported
   * 5: cross - not supported
   */

  /* SVDI only supports a dot, which is the default */

  /* An unsupported index is mapped to a supported value, such
   * as the default, without error. So can just ignore this
   * function.
   */

} /* end xcmkt */

/* MARKER SIZE */
static void xcmks(anything **params, int num_surfaces, anything **surf_list)
{
  /* marker size is mapped to nearest avail marker size on
   * the device. SVDI only supports 1 marker size, which is
   * the default size, so this function not supported.
   */
  /* error 4:001 - Function not supported. Function ignored */
  report_error(cur_state, 4, 001, *(short *)params[0]);

} /* end xcmks */

/* MARKER COLOR */
static void xcmkc(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* if color selection mode is in error, ignore function */
    if (cur_state->csm == -1) {
      break;
    }

    if (cur_state->csm == CINDEX) { /* indexed mode */

      /* let SVDI catch out of range indices */
      cur_state->mark_csm      = CINDEX;
      cur_state->mark_color[0] = *(int *)params[1];
    } /* end indexed mode  */

    else { /* direct color mode */

      /* let SVDI catch illegal color values */
      cur_state->mark_csm      = CDRECT;
      cur_state->mark_color[0] = *((int *)params[1] + 0);
      cur_state->mark_color[1] = *((int *)params[1] + 1);
      cur_state->mark_color[2] = *((int *)params[1] + 2);
    } /* end direct color mode */

  } /* end for each surface */
} /* end xcmkc */

/* TEXT PRECISION */
static void xctxp(anything **params, int num_surfaces, anything **surf_list)
{
  /* error 4:001 - Function not supported. Function ignored */
  report_error(cur_state, 4, 001, *(short *)params[0]);

} /* end xctxp */

/* TEXT COLOR */
static void xctxc(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* if color selection mode is in error, ignore function */
    if (cur_state->csm == -1) {
      break;
    }

    if (cur_state->csm == CINDEX) { /* indexed mode */

      /* let SVDI catch out of range indices */
      cur_state->text_csm      = CINDEX;
      cur_state->text_color[0] = *(int *)params[1];
    } /* end indexed mode */

    else { /* direct color mode */

      /* let SVDI catch illegal color values */
      cur_state->text_csm      = CDRECT;
      cur_state->text_color[0] = *((int *)params[1] + 0);
      cur_state->text_color[1] = *((int *)params[1] + 1);
      cur_state->text_color[2] = *((int *)params[1] + 2);
    } /* end direct mode */

  } /* end for each surface */
} /* end xctxc */

/* CHARACTER HEIGHT */
static void xcchh(anything **params, int num_surfaces, anything **surf_list)
{
  int   i;              /* index for loop on surfaces */
  float vdi_cs;         /* SVDI character size */
  float temp_array[14]; /* used for SVDI inquire */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* set SVDI character size if needed */
    if (*(float *)params[1] != cur_state->vdi_attrib.char_x) {

      /* let SVDI catch illegal character height */
      cur_state->vdi_attrib.char_x = *(float *)params[1];

      /* map to NDC */
      vdi_cs = *(float *)params[1] * fabs(cur_state->xscale);
      vdstcs(&vdi_cs);

      /* inquire to find out what it really got set to */
      vdiqos(&temp_array[1]);

      /* map back to VDC and store in state table */
      cur_state->char_height = temp_array[6] / fabs(cur_state->xscale);

    } /* end if char size needs to be set */

  } /* end for each surface */
} /* end xcchh */

/* CHARACTER ORIENTATION */
static void xccho(anything **params, int num_surfaces, anything **surf_list)
{
  /* error 4:001 - Function not supported. Function ignored */
  report_error(cur_state, 4, 001, *(short *)params[0]);

} /* end xccho */

/* INTERIOR STYLE */
static void xcis(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  /* CGI interior styles: */
  /*  hollow - supported
   *  solid -  supported
   *  pattern - not supported
   *  hatch - not supported
   *  empty - not supported
   *  bitmap - not supported
   */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* check for legal enumerated value */
    if (*(int *)params[1] == CHOLLO || *(int *)params[1] == CSOLID) {
      cur_state->interior_style = *(int *)params[1];
    }
    else {
      /* error 1:-101 Enumerated parameter out of range. Function ignored */
      report_error(cur_state, 1, -101, *(short *)params[0]);
    }

  } /* end for each surface */
} /* end xcis */

/* FILL COLOR */
static void xcflc(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* if color selection mode is in error, ignore function */
    if (cur_state->csm == -1) {
      break;
    }

    if (cur_state->csm == CINDEX) { /* indexed mode */

      /* let SVDI catch out of range indices */
      cur_state->fill_csm      = CINDEX;
      cur_state->fill_color[0] = *(int *)params[1];
    } /* end indexed color */

    else { /* direct color mode */

      /* let SVDI catch illegal color values */
      cur_state->fill_csm      = CDRECT;
      cur_state->fill_color[0] = *((int *)params[1] + 0);
      cur_state->fill_color[1] = *((int *)params[1] + 1);
      cur_state->fill_color[2] = *((int *)params[1] + 2);
    } /* end direct color mode */

  } /* end for each surface */
} /* end xcflc */

/* COLOR SELECTION MODE */
static void xccsm(anything **params, int num_surfaces, anything **surf_list)
{
  int i; /* index for loop on surfaces */

  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* check enumerated parameters */
    if (*(int *)params[1] != CDRECT && *(int *)params[1] != CINDEX) {
      /* error 1:-101 enumerated parameter out of range. Function ignored */
      report_error(cur_state, 1, -101, *(short *)params[0]);
    }
    else {

      /* if the device supports both index and direct, or the
         request if only for index, everything is honky-dory */
      if (dev_descrip.csm_avail == CCLRID || *(int *)params[1] == CINDEX) {
        cur_state->csm = *(int *)params[1];
      }
      else { /* must be an error */
        /* error 3:305 - Color selection mode not supported. All
          functions with affected parameters ignored until reset */
        report_error(cur_state, 3, 305, *(short *)params[0]);
        cur_state->csm = -1; /* flag error condition */
      } /* end else must be an error */
    } /* end else */

  } /* end for each surface */
} /* end xccsm */

/* COLOR TABLE */
static void xcct(anything **params, int num_surfaces, anything **surf_list)
{
  int   i;                   /* index for loop on surfaces */
  int   j, k;                /* indices for looping */
  int   starti;              /* starting color index */
  int   num_cols;            /* number of colors in color list */
  int   tmp_array[3];        /* temp array used for setting bg color */
  float color_array[256][3]; /* used for setting SVDI color table */
  int   num_set;             /* number of colors being set, for SVDI */
  int   maxindex;            /* max color index to set */
  int   indx_ptr;            /* for keeping track of color indices */
  int   index1, index2;      /* defines a range of indices to set */
  int   first;               /* marks first time through loop */
  int   one = 1;

  first = TRUE;

  /* loop through surfaces */
  for (i = 0; i < num_surfaces; ++i) {

    cur_state = (surf_statelist *)surf_list[i];
    if (cur_state->cgi_inited != CYES) {
      break;
    }

    /* starting color index and number of colors */
    starti   = *(int *)params[1];
    num_cols = *(int *)params[2];

    /* error checking */
    /* an out of range index in an index definition function should
       be ignored.  should an error be reported ?? */
    if (starti < 0) {
      starti = 0;
    }

    if ((starti + num_cols - 1) > dev_descrip.num_cols) {
      /* error 6:308 - Too many color specifiers. Discard extras */
      report_error(cur_state, 6, 308, *(short *)params[0]);
      num_cols = dev_descrip.num_cols - starti + 1;
    } /* end error */

    /* convert rgb to lie between 0. and 1. */
    /* ...only do this once */
    if (first) {
      k = 0;
      for (j = starti; j < starti + num_cols; j++) {
        color_array[j][0] = (float)((int *)params[3])[k++] / 255.0f;
        color_array[j][1] = (float)((int *)params[3])[k++] / 255.0f;
        color_array[j][2] = (float)((int *)params[3])[k++] / 255.0f;
      }
      first = FALSE;
    }

    /* mark that the color table has been set */
    cur_state->color_set = TRUE;

    /* update the current state color table */
    k = 0;
    for (j = starti; j < starti + num_cols; j++) {
      cur_state->color_table[j].r = ((int *)params[3])[k++];
      cur_state->color_table[j].g = ((int *)params[3])[k++];
      cur_state->color_table[j].b = ((int *)params[3])[k++];
    }

    /* set the SVDI color table */
    /* color indices are mapped so that fg and bg colors are
     * 'reasonable' colors.  If the user sets either the fg index
     * color or bg color index, map those colors so that the
     * fg and bg colors won't change
     */
    /* color index mappings:
     *   if set index 0  --> set bg_index
     *   if set index 1  --> set fg_index
     *   if set bg_index --> set 0
     *   if set fg_index --> set 1
     */

    /* set the special cases first */
    /* 1. caller setting index 0  - set the background color */
    /* --set_background_color() will set color table entry bg_index */
    if (starti == 0) {
      tmp_array[0] = cur_state->color_table[0].r;
      tmp_array[1] = cur_state->color_table[0].g;
      tmp_array[2] = cur_state->color_table[0].b;
      set_background_color(cur_state, tmp_array);
    }

    /* 2. caller setting index 1  - set index fg_index */
    if (starti <= 1 && (starti + num_cols > 1)) {
      vdstco(&one, &dev_descrip.index_array[cur_state->fg_index], &color_array[1],
             &dev_descrip.col_mode);
    }

    /* 3. caller setting index bg_index  - set index 0 */
    if (cur_state->bg_index != 0) { /* .. don't do it twice */
      if (starti <= cur_state->bg_index && (starti + num_cols > cur_state->bg_index)) {
        vdstco(&one, &dev_descrip.index_array[0], &color_array[cur_state->bg_index],
               &dev_descrip.col_mode);
      }
    }

    /* 4. caller setting index fg_index  - set index 1 */
    if (cur_state->fg_index != 1) { /* .. don't do it twice */
      if (starti <= cur_state->fg_index && (starti + num_cols > cur_state->fg_index)) {
        vdstco(&one, &dev_descrip.index_array[1], &color_array[cur_state->fg_index],
               &dev_descrip.col_mode);
      }
    }

    /* now do all the rest */
    /* ...sort for convenience */
    if (cur_state->fg_index < cur_state->bg_index) {
      index1 = cur_state->fg_index;
      index2 = cur_state->bg_index;
    }
    else {
      index1 = cur_state->bg_index;
      index2 = cur_state->fg_index;
    } /* end sort */

    maxindex = starti + num_cols - 1;
    indx_ptr = starti;
    if (indx_ptr < 2) {
      indx_ptr = 2;
    }

    /* 1. set color table between index 1 and index1 */
    if (indx_ptr <= 2) {
      num_set = min(index1 - 2, maxindex - 1);
      if (num_set > 0) {
        vdstco(&num_set, &dev_descrip.index_array[2], &color_array[2], &dev_descrip.col_mode);
        indx_ptr = index1 + 1;
      }
    }

    if (indx_ptr == index1) { /* update indx_ptr */
      indx_ptr = indx_ptr + 1;
    }

    /* 2. set color table between index1 and index2 */
    if (indx_ptr < index2) {
      num_set = min(index2 - indx_ptr, maxindex - indx_ptr + 1);
      if (num_set > 0) {
        vdstco(&num_set, &dev_descrip.index_array[indx_ptr], &color_array[indx_ptr],
               &dev_descrip.col_mode);
        indx_ptr = index2 + 1;
      }
    }

    if (indx_ptr == index2) { /* update indx_ptr */
      indx_ptr = indx_ptr + 1;
    }

    /* 3. set color table after index 2 */
    if (indx_ptr >= index2 + 1) {
      num_set = maxindex - indx_ptr + 1;
      if (num_set > 0) {
        vdstco(&num_set, &dev_descrip.index_array[indx_ptr], &color_array[indx_ptr],
               &dev_descrip.col_mode);
      }
    }
  } /* end for each surface */
} /* end xcct */

/* GET TEXT EXTENT */
static void xcgtxx(anything **params, anything **surf_list)
{
  float act_h, act_w;   /* actual char height and width (VDC) */
  float x1, y1, x2;     /* parrallelogram parts */
  float save, temp_h;   /* temporary variables */
  float temp_array[14]; /* temp array used for SVDI inquires */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* find out actual height and width from SVDI for this size */
  /* ..save the old value */
  vdiqos(&temp_array[1]);
  save = temp_array[6];

  /* ..set current height value */
  temp_h = cur_state->char_height * fabs(cur_state->xscale);
  vdstcs(&temp_h);

  /* ..find out actual height and width */
  vdiqos(&temp_array[1]);
  act_h = temp_array[6] / fabs(cur_state->xscale);
  act_w = temp_array[7] / fabs(cur_state->xscale);

  /* set it back to original */
  vdstcs(&save);

  /* parallelogram */
  /*  x1, y1 = bottom left corner of string */
  x1                  = *(float *)params[1];
  y1                  = *(float *)params[2];
  *(float *)params[8] = x1;
  *(float *)params[9] = y1;

  /*  x2, y2 = bottom right corner of string */
  /* params[4] is an added parameter that contains string size */
  x2                   = x1 + (*(int *)params[16] * act_w);
  *(float *)params[10] = x2;
  *(float *)params[11] = y1;

  /*  x3, y3 = top right corner of string ( x3 = x2 )*/
  *(float *)params[12] = x2;
  *(float *)params[13] = y1 + act_h;

  /*  x4, y4 = top left corner of string ( y4 = y3 )*/
  *(float *)params[14] = x1;
  *(float *)params[15] = y1 + act_h;

  /* concatenation point is bottom right corner ( x2, y2 )*/
  *(float *)params[6] = x2;
  *(float *)params[7] = y1;

  /* validity flags */
  *(int *)params[4] = CVAL;
  *(int *)params[5] = CVAL;

} /* end xcgtxx */

/* INQUIRE PRIMITIVE SUPPORT LEVELS */
static void xcqprl(anything **params, anything **surf_list)
{
  int        qdc_index = 25; /* index for inquiries to vdiqdc  */
  float      value;          /* value returned by vdiqdc */
  static int set = FALSE;    /* flag whether values have been set */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* max points for a polyline */
    dev_descrip.max_pts_polyline = MAX_POLYLINE;

    /* max points for a disjoint polyline */
    dev_descrip.max_pts_disj_polyline = MAX_DJ_POLYLINE;

    /* max points for a polygon */
    vdiqdc(&qdc_index, &value);
    dev_descrip.max_pts_polygon = (value == 99999.) ? MAX_POLYGON : min((int)value, MAX_POLYGON);

    /* max points for a polygon set */
    /* --not supported */
    dev_descrip.max_pts_polygon_set = 0;

    /* max points for a polymarker */
    dev_descrip.max_pts_polymarker = MAX_POLYMARKER;

    /* max number of points for a closed figure */
    /* --not supported */
    dev_descrip.max_pts_closed_fig = 0;

    /* max characters of text  */
    dev_descrip.max_chars_text = MAX_TEXT_LENGTH;

    /* max number of cell color specifiers for cell array */
    dev_descrip.max_pts_cellarray = -1;

    /* cell array fill capability */
    dev_descrip.cellarray_fill_cap = CFILLD;

    /* cell array alignment */
    dev_descrip.cellarray_align_cap = CAXIS;

    /* compound text capability */
    dev_descrip.compound_text_cap = CCNONE;

    /* closed figure capability */
    dev_descrip.compound_fig_cap = CCNONE;

  } /* end if not set */

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* stuff all the values into the return parameters */
  *(int *)params[2]  = dev_descrip.max_pts_polyline;
  *(int *)params[3]  = dev_descrip.max_pts_disj_polyline;
  *(int *)params[4]  = dev_descrip.max_pts_polygon;
  *(int *)params[5]  = dev_descrip.max_pts_polygon_set;
  *(int *)params[6]  = dev_descrip.max_pts_polymarker;
  *(int *)params[7]  = dev_descrip.max_pts_closed_fig;
  *(int *)params[8]  = dev_descrip.max_chars_text;
  *(int *)params[9]  = dev_descrip.max_pts_cellarray;
  *(int *)params[10] = dev_descrip.cellarray_fill_cap;
  *(int *)params[11] = dev_descrip.cellarray_align_cap;
  *(int *)params[12] = dev_descrip.compound_text_cap;
  *(int *)params[13] = dev_descrip.compound_fig_cap;

} /* end xcqprl */

/* INQUIRE LINE CAPABILITY */
static void xcqln(anything **params, anything **surf_list)
{

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* set the validity flag */
  *(int *)params[1] = CVAL;

  /* number of predefined line bundles - bundles not supported */
  *(int *)params[2] = 0;

  /* number of settable line bundles - bundles not supported */
  *(int *)params[3] = 0;

  /* max line bundle index - bundles not supported */
  *(int *)params[4] = 0;

  /* dynamic modification for line bundles- bundles not supported */
  *(int *)params[5] = 0;

  /* nom, min, max, get set in set_dev_descrip() at initialization */
  *(int *)params[6] = dev_descrip.linewidth_nominal;
  *(int *)params[7] = dev_descrip.linewidth_min;
  *(int *)params[8] = dev_descrip.linewidth_max;

} /* end xcqln */

/* INQUIRE LIST OF AVAILABLE LINE TYPES */
static void xcqlnt(anything **params, anything **surf_list)
{
  int        i, j;          /* loop indices */
  int        qdc_index = 6; /* index for inquiries to vdiqdc  */
  float      value;         /* value returned by vdiqdc */
  static int set = FALSE;   /* flag whether values have been set */
  static int ntotal;        /* save total nbr of linestyles */

  static unsigned mask = ~(~0u << 1);

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* inquire supported linestyles from SVDI */
    vdiqdc(&qdc_index, &value);

    /* SVDI always supports linestyle 0 */
    ntotal                           = 0;
    dev_descrip.line_types[ntotal++] = 0;

    /* find all SVDI linestyles supported by looking at bits */
    for (i = 0; i < 5; i++) {
      if (((int)value & mask << i) >> i == 1) {
        dev_descrip.line_types[ntotal++] = i + 1;
      }
    }
    /* map to CGI linestyles */
    /* cgi linestyles: 1 - solid 2 - dash 3 - dot,
     *                 4 - dashdot 5 - dash dot dot
     * vdi linestyles: 0 - solid 1 - dotted 2 - dot dash
     *                 3 - short dash 4 - long dash 5 - medium dash
     */
    for (i = 0; i < ntotal; i++) {
      switch (dev_descrip.line_types[i]) {
      case 0: dev_descrip.line_types[i] = 1; break;
      case 1: dev_descrip.line_types[i] = 3; break;
      case 2: dev_descrip.line_types[i] = 4; break;
      /* ...all SVDI dashes map to cgi dash */
      case 3: dev_descrip.line_types[i] = 2; break;
      case 4: dev_descrip.line_types[i] = 2; break;
      case 5: dev_descrip.line_types[i] = 2; break;
      default: break;
      } /* end switch */
    }

    /* take out duplicates */
    for (i = 0; i < ntotal; i++) {
      if (dev_descrip.line_types[i] == 2) {
        break;
      }
      ntotal = i + 1;
    }
  } /* end if not set */

  /* set the total */
  *(int *)params[4] = ntotal;

  /* do some error checking */
  if (*(int *)params[2] < 1) {
    /* error - set validity flag to invalid */
    *(int *)params[3] = CINVAL;
  }
  else { /* send back only the ones requested -  zero out the rest */

    i                 = *(int *)params[2] - 1;
    j                 = 0;
    *(int *)params[5] = 0; /* actual number of list elements sent back */
    while (j < ntotal) {
      if (j < *(int *)params[1] && i < ntotal) {
        ((int *)params[6])[j++] = dev_descrip.line_types[i++];
        (*(int *)params[5])++;
      }
      else {
        ((int *)params[6])[j++] = 0;
      }
    }
    /* set validity flag */
    *(int *)params[3] = CVAL;
  } /* end else send back */

} /* end xcqlnt */

/* INQUIRE LIST OF AVAILABLE SCALED LINE WIDTHS */
static void xcqslw(anything **params, anything **surf_list)
{
  int        qdc_index = 5; /* index for inquiries to vdiqdc  */
  float      value;         /* value returned by vdiqdc */
  static int set = FALSE;   /* flag whether values have been set */
  static int ntotal;        /* save total nbr of linewidths */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* number of linewidths available */
    vdiqdc(&qdc_index, &value);
    ntotal = (int)value;

  } /* end if not set */

  /* validity flag */
  *(int *)params[3] = CVAL;
  *(int *)params[4] = ntotal;

  /* nlist = 0 means empty list ( continuous range ) */
  /* --leave list (params[6]) untouched */
  *(int *)params[5] = 0;

} /* end xcqslw */

/* INQUIRE MARKER CAPABILITY */
static void xcqmk(anything **params, anything **surf_list)
{
  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* set the validity flag */
  *(int *)params[1] = CVAL;

  /* number of predefined marker bundles - bundles not supported */
  *(int *)params[2] = 0;

  /* number of settable marker bundles - bundles not supported */
  *(int *)params[3] = 0;

  /* max marker bundle index - bundles not supported */
  *(int *)params[4] = 0;

  /* dynamic modification for marker bundles- bundles not supported */
  *(int *)params[5] = 0;

  /* nom, min, max, get set in set_dev_descrip() at initialization */
  *(int *)params[6] = dev_descrip.mark_nominal;
  *(int *)params[7] = dev_descrip.mark_min;
  *(int *)params[8] = dev_descrip.mark_max;

} /* end xcqmk */

/* INQUIRE LIST OF AVAILABLE MARKER TYPES */
static void xcqmkt(anything **params, anything **surf_list)
{
  static int set = FALSE; /* flag whether values have been set */
  static int ntotal;      /* save total nbr of linestyles */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* SVDI only supports 1 marker type - "dot" */
    ntotal                    = 1;
    dev_descrip.mark_types[0] = 1;

  } /* end if not set */

  /* set validity flag */
  *(int *)params[3] = CVAL;

  *(int *)params[4] = ntotal;

  if (*(int *)params[1] > 0) {
    *(int *)params[5]     = 1;
    ((int *)params[6])[0] = dev_descrip.mark_types[0];
  }
  else {
    *(int *)params[5] = 0;
  }

} /* end xcqmkt */

/* INQUIRE LIST OF AVAILABLE SCALED MARKER SIZES */
static void xcqsms(anything **params, anything **surf_list)
{
  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* validity flag */
  *(int *)params[3] = CVAL;

  /* there is only 1 size */
  *(int *)params[4]     = 1;
  *(int *)params[5]     = 1;
  ((int *)params[6])[0] = dev_descrip.mark_nominal;

} /* end xcqsms */

/* INQUIRE LIST OF AVAILABLE CHARACTER HEIGHTS */
static void xcqchh(anything **params, anything **surf_list)
{
  int        qdc_index = 7; /* index for inquiries to vdiqdc  */
  float      value;         /* value returned by vdiqdc */
  static int set = FALSE;   /* flag whether values have been set */
  static int ntotal;        /* save total nbr of char heights */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* number of char sizes */
    vdiqdc(&qdc_index, &value);
    ntotal = (int)value;

  } /* end if not set */

  /* font - SRCP does not support multiple fonts; use null string */
  /* I don't check for font, anything goes */

  /* only string precision is supported */
  if (*(int *)params[2] == CSTRNG) {

    /* total number of character heights */
    *(int *)params[6] = ntotal;

    /* set nlist = 0 means empty list ( continuous range ) */
    /* --leave list (params[8]) untouched */
    *(int *)params[7] = 0;

    /* validity flag */
    *(int *)params[5] = CVAL;
  } /* end if */

  else {
    *(int *)params[5] = CINVAL;
  }

} /* end xcqchh */

/* INQUIRE FILL CAPABILITY */
static void xcqfl(anything **params, anything **surf_list)
{
  static int set = FALSE; /* flag whether values have been set */
  static int ntotal;      /* save total nbr of fill styles */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set */

    set = TRUE;

    ntotal                         = 2;
    dev_descrip.interior_styles[0] = CHOLLO;
    dev_descrip.interior_styles[1] = CSOLID;

  } /* end if not set */

  /* set the validity flag */
  *(int *)params[1] = CVAL;

  /* number of predefined fill bundles - bundles not supported */
  *(int *)params[2] = 0;

  /* number of settable fill bundles - bundles not supported */
  *(int *)params[3] = 0;

  /* max fill bundle index - bundles not supported */
  *(int *)params[4] = 0;

  /* dynamic modification for fill bundles- bundles not supported */
  *(int *)params[5] = 0;

  /* interior styles */
  *(int *)params[6]     = ntotal;
  ((int *)params[7])[0] = dev_descrip.interior_styles[0];
  ((int *)params[7])[1] = dev_descrip.interior_styles[1];

  /* number of predefined patterns - patterns not supported */
  *(int *)params[8] = 0;

  /* number of settable patterns - patterns not supported */
  *(int *)params[9] = 0;

  /* maximum pattern index - patterns not supported */
  *(int *)params[10] = 0;

  /* preferred pattern size divisor - patterns not supported */
  *(int *)params[11] = 0;

  /* maximum pattern size - patterns not supported */
  *(int *)params[12] = 0;
  *(int *)params[13] = 0;

  /* pattern transformation support - patterns not supported */
  *(int *)params[14] = CPTNO;

} /* end xcqfl */

/* INQUIRE COLOR CAPABILITIES */
static void xcqc(anything **params, anything **surf_list)
{
  int        qdc_index;   /* index for inquiries to vdiqdc  */
  float      value;       /* value returned by vdiqdc */
  static int set = FALSE; /* flag whether values have been set */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  if (!set) { /* values have't been set, inquire SVDI */

    set = TRUE;

    /* number of simultaneous colors */
    dev_descrip.num_simul_colors = dev_descrip.num_cols;

    /* number of available colors */
    qdc_index = 27;
    vdiqdc(&qdc_index, &value);
    dev_descrip.num_avail_colors = (int)value;

    /* number of available intensities */
    /* --see 7.1.4 about direct color */
    qdc_index = 3;
    vdiqdc(&qdc_index, &value);
    dev_descrip.num_avail_int = (int)value;

    /* dynamic modification accepted for color table */
    dev_descrip.dynamic_mod_ct = CIRG;

    /* color overwrite capability */
    /* --what is the answer?? */
    /* dev_descrip. color_overwrite = ??; */

    /* monochromatic device */
    qdc_index = 32;
    vdiqdc(&qdc_index, &value);
    dev_descrip.monochrome_device = (value == 0.0) ? CYES : CNO;

  } /* end not set */

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* return the answers */
  *(int *)params[2]     = dev_descrip.num_simul_colors;
  *(int *)params[3]     = dev_descrip.num_avail_colors;
  ((int *)params[4])[0] = dev_descrip.num_avail_int;
  *(int *)params[5]     = dev_descrip.csm_avail; /* set in set_dev_descrip()*/
  *(int *)params[6]     = dev_descrip.dynamic_mod_ct;
  *(int *)params[7]     = dev_descrip.color_overwrite;
  *(int *)params[8]     = dev_descrip.monochrome_device;

} /* end xcqc */

/* INQUIRE LINE ATTRIBUTES */
static void xcqlna(anything **params, anything **surf_list)
{

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* line bundle index - bundles not supported */
  *(int *)params[2] = 0;

  /* line type */
  *(int *)params[3] = cur_state->line_type;

  /* line width specification mode - only support scaled */
  *(int *)params[4] = CSCA;

  /* line width */
  *(float *)params[5] = cur_state->line_width;

  /* color selection mode in which line color was last specified */
  *(int *)params[6] = cur_state->line_csm;

  /* line color in color selection mode last specified */
  ((int *)params[7])[0] = cur_state->line_color[0];
  if (cur_state->line_csm == CDRECT) {
    ((int *)params[7])[1] = cur_state->line_color[1];
    ((int *)params[7])[2] = cur_state->line_color[2];
  }

  /* line clipping mode */
  *(int *)params[8] = CLOCUS;

} /* end xcqlna */

/* INQUIRE MARKER ATTRIBUTES */
static void xcqmka(anything **params, anything **surf_list)
{

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* marker bundle index - bundles not supported */
  *(int *)params[2] = 0;

  /* marker type */
  *(int *)params[3] = cur_state->mark_type;

  /* marker size specification mode - only support scaled */
  *(int *)params[4] = CSCA;

  /* marker size */
  *(float *)params[5] = cur_state->mark_size;

  /* color selection mode in which marker color was last specified */
  *(int *)params[6] = cur_state->mark_csm;

  /* marker color in color selection mode last specified */
  ((int *)params[7])[0] = cur_state->mark_color[0];
  if (cur_state->mark_csm == CDRECT) {
    ((int *)params[7])[1] = cur_state->mark_color[1];
    ((int *)params[7])[2] = cur_state->mark_color[2];
  }

  /* marker clipping mode */
  *(int *)params[8] = CLOCUS;

} /* end xcqmka */

/* INQUIRE TEXT ATTRIBUTES */
static void xcqtxa(anything **params, anything **surf_list)
{

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* text bundle index - bundles not supported */
  *(int *)params[2] = 0;

  /* font index - only the default font is supported */
  *(int *)params[3] = 1;

  /* font precision */
  *(int *)params[4] = CSTRNG;

  /* character expansion factor */
  *(float *)params[5] = 1.0;

  /* character spacing */
  *(float *)params[6] = 0.0;

  /* color selection mode in which text color was last specified */
  *(int *)params[7] = cur_state->text_csm;

  /* text color in color selection mode last specified */
  ((int *)params[8])[0] = cur_state->text_color[0];
  if (cur_state->text_csm == CDRECT) {
    ((int *)params[8])[1] = cur_state->text_color[1];
    ((int *)params[8])[2] = cur_state->text_color[2];
  }

  /* character height  */
  *(float *)params[9] = cur_state->char_height;

  /* character orientation */
  ((float *)params[10])[0] = 0.;
  ((float *)params[10])[1] = 1.;
  ((float *)params[10])[2] = 1.;
  ((float *)params[10])[3] = 0.;

  /* text path */
  *(int *)params[11] = CTPRIT;

  /* horizontal alignment */
  *(int *)params[12] = CTANH;

  /* continuous horizontal alignment */
  *(float *)params[13] = 1.0;

  /* vertical alignment */
  *(int *)params[14] = CTANV;

  /* continuous horizontal alignment */
  *(float *)params[15] = 1.0;

  /* character set index ?? */
  *(int *)params[16] = 1;

  /* alternate character set index ?? */
  *(int *)params[17] = 1;

} /* end xcqtxa */

/* INQUIRE FILL ATTRIBUTES */
static void xcqfla(anything **params, anything **surf_list)
{

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* validity flag */
  *(int *)params[1] = CVAL;

  /* fill bundle index - bundles not supported */
  *(int *)params[2] = 0;

  /* interior style */
  *(int *)params[3] = cur_state->interior_style;

  /* color selection mode in which fill color was last specified */
  *(int *)params[4] = cur_state->fill_csm;

  /* fill color in color selection mode last specified */
  ((int *)params[5])[0] = cur_state->fill_color[0];
  if (cur_state->fill_csm == CDRECT) {
    ((int *)params[5])[1] = cur_state->fill_color[1];
    ((int *)params[5])[2] = cur_state->fill_color[2];
  }

  /* hatch index - not supported */
  *(int *)params[6] = 0;

  /* pattern index - not supported */
  *(int *)params[7] = 0;

  /* --i don't know what any of the following means */
  /* fill bitmap identifier */
  *(int *)params[8] = 0;

  /* fill bitmap region */
  *(float *)params[9]  = 0.;
  *(float *)params[10] = 0.;
  *(float *)params[11] = 0.;
  *(float *)params[12] = 0.;

  /* fill reference point */
  *(float *)params[13] = 0.;
  *(float *)params[14] = 0.;

  /* pattern orientation and size */
  *(float *)params[15] = 0.;
  *(float *)params[16] = 0.;
  *(float *)params[17] = 0.;
  *(float *)params[18] = 0.;

} /* end xcqfla */

/* INQUIRE LIST OF COLOR TABLE ENTRIES */
static void xcqcte(anything **params, anything **surf_list)
{
  int i;     /* loop index */
  int first; /* index of first element to return */
  int nreq;  /* number of list elements requested */
  int j = 0; /* array index */

  /* there is only one surface for inquiries */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* number of requested elements, and the first index */
  nreq  = *(int *)params[1];
  first = *(int *)params[2];

  /* error checking */
  if (first < 1) {
    first = 1;
  }

  if ((first + nreq - 1) > dev_descrip.num_cols) {
    nreq = dev_descrip.num_cols - first + 1;
  }

  /* actual number of colors returned */
  *(int *)params[5] = nreq;

  for (i = first; i < (first + *(int *)params[5]); i++) {

    ((int *)params[6])[j++] = cur_state->color_table[i - 1].r;
    ((int *)params[6])[j++] = cur_state->color_table[i - 1].g;
    ((int *)params[6])[j++] = cur_state->color_table[i - 1].b;
  } /* end for i */

  /* return total number of colors */
  *(int *)params[4] = dev_descrip.num_cols;

  /* validity flag */
  *(int *)params[3] = CVAL;

} /* end xcqcte */

/* INITIALIZE LOGICAL INPUT DEVICE */
static void xcili(anything **params, anything **surf_list)
{
  /* there is only one surface for input */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* only support locator  */
  if (*(int *)params[1] != CLOCAT) {
    /* error 2:001 - Request for an unsupported feature. Function ignored */
    report_error(cur_state, 2, 001, *(short *)params[0]);
  }
  else if (*(int *)params[2] != 1) {
    /* error 3:501 - LID index exceeds number of this class of LID.
       Function ignored */
    report_error(cur_state, 3, 501, *(short *)params[0]);
  }
  else if (dev_descrip.dev_class == CINPUT ||
           dev_descrip.dev_class == COUTIN) { /* device does input */

    cur_state->input_dev_class = CLOCAT;
    cur_state->input_dev_index = 1;
    cur_state->input_dev_state = CREADY;

  } /* end if device does input */

} /* end xcili */

/* REQUEST LOCATOR */
static void xcrqlc(anything **params, anything **surf_list)
{
  float x, y; /* x,y values of locator */

  /* there is only one surface for input */
  cur_state = (surf_statelist *)surf_list[0];
  if (cur_state->cgi_inited != CYES) {
    return;
  }

  /* check if device has been initialized */
  /* ...only support input device index 1 */
  if (*(int *)params[1] != 1 || cur_state->input_dev_state != CREADY) {
    /* error 5:501 - Function illegal in LID state RELEASED.
       Function ignored  */
    report_error(cur_state, 5, 501, *(short *)params[0]);
    *(int *)params[3] = CINVAL;
  }

  else if (dev_descrip.dev_class == CINPUT ||
           dev_descrip.dev_class == COUTIN) { /* device does input */

    /* for now only allow limited timeout */
    /* positive and negative values = wait forever, 0 = no wait */
    if (*(int *)params[2] == 0) {
      *(int *)params[3] = CINVAL;
      return;
    }

    /* get location and trigger */
    vdakgl((int *)params[6], &x, &y);

    /* convert location from NDC to VDC */
    *(float *)params[7] =
        (x - (cur_state->vp1.x * dev_descrip.xndc_max)) / cur_state->xscale + cur_state->vdc1.x;

    *(float *)params[8] =
        (y - (cur_state->vp1.y * dev_descrip.yndc_max)) / cur_state->yscale + cur_state->vdc1.y;

    /* response validity */
    *(int *)params[3] = CVAL;

    /* request status, SVDI only has trigger */
    *(int *)params[4] = CTRIGR;

    /* measure validity */
    /* --i'm not sure when measure validity is invalid */
    *(int *)params[5] = CVAL;

  } /* end else device does input */

  else { /* device doesn't do input */
    /* response validity */
    *(int *)params[3] = CVAL;
  }

} /* end xcrqlc */

/*-------------------------------------------------------------*/
/* >> UTILITY ROUTINES                                         */
/*-------------------------------------------------------------*/

/* init_state */
/* Initialize the state list */
static void init_state(surf_statelist *surf_state)
/* surface statelist to init */
{
  int   i;               /* loop index */
  float tmp_array[1][3]; /* temp array to set SVDI color table */
  int   one = 1;

  /* set default error control: class 1 : DETECTION OFF
     class 2 - MAX_ERROR_CLASS : REPORTING OFF */
  surf_state->err_flag[0] = CEHDOF;
  for (i = 1; i < MAX_ERROR_CLASS; i++) {
    surf_state->err_flag[i] = CEHROF;
  }

  /* initialize to 0 errors */
  surf_state->err_head_ptr = 0;
  surf_state->err_count    = 0;

  /* default vdc extent */
  surf_state->vdc1.x = 0.0;
  surf_state->vdc1.y = 0.0;
  surf_state->vdc2.x = 32767.;
  surf_state->vdc2.y = 32767.;

  /* viewport  - default viewport specification mode is FRACTION OF
   * DRAWING SURFACE.  Thus, the viewport (0,0),(1,1) corresponds to
   * the lower and upper corners, respectively, of the default
   * device viewport, which is the largest unrotated rectangular area
   * visible on the drawing surface.
   */
  /* default viewport is (0,0),(1,1) */
  surf_state->vp1.x = 0.0;
  surf_state->vp1.y = 0.0;
  surf_state->vp2.x = 1.0;
  surf_state->vp2.y = 1.0;

  /* VDC to viewport in NDC mapping */
  /* --device viewport mapping defaults are: isotropy forced,
     horizontal alignment left, vertical alignment bottom */
  /* -- this routine sets effective viewport in VDC and
     NDC (forcing isotropic mapping) and the mapping scale */
  set_mapping(surf_state);

  /* default clip indicator is off (VDC clipping) */
  surf_state->clip_indicator = COFF;

  /* default view surface clip indicator is viewport */
  surf_state->ds_clip_indicator = CVPORT;

  /* default clip rectangle is (0,0) (32767,32767) */
  /* ...can't be changed in SS-CGI */
  surf_state->clip_rect1.x = 0.0;
  surf_state->clip_rect1.y = 0.0;
  surf_state->clip_rect2.x = 32767.;
  surf_state->clip_rect2.y = 32767.;

  /* effective clip rectangle is intersection of VDC extent and
     clip rectangle */
  surf_state->eff_clip_rect1.x = 0.0;
  surf_state->eff_clip_rect1.y = 0.0;
  surf_state->eff_clip_rect2.x = 32767.;
  surf_state->eff_clip_rect2.y = 32767.;

  /* set the clip region */
  set_clipping(surf_state);

  /* default color mode is indexed */
  surf_state->csm = CINDEX;

  /* default line, fill, text, and marker color indices */
  surf_state->line_color[0] = 1;
  surf_state->fill_color[0] = 1;
  surf_state->text_color[0] = 1;
  surf_state->mark_color[0] = 1;

  /* default line type is solid */
  surf_state->line_type = 1;

  /* default interior style is hollow  */
  surf_state->interior_style = CHOLLO;

  /* char size is in VDC - default is 1/100 of default VDC extent */
  surf_state->char_height = 327.;

  /* default line width specification mode is scaled, default
     linewidth is 1. */
  surf_state->line_width = 1.;

  /* default marker specification mode is scaled, default
     marker size is 1. */
  surf_state->mark_size = 1.;

  /* the only default colors defined in CGI are index 0 = background,
   * and index 1 = foreground. Use SVDI's bg and fg until the user
   * sets the color table.
   */
  if (dev_descrip.svdi_type == 0) { /* vector SVDI */

    /* set up the fg and bg mappings */
    surf_state->bg_index = dev_descrip.att_array[2];
    surf_state->fg_index = dev_descrip.att_array[1];

    /* set the current surface color table */
    vdiqco(&one, &dev_descrip.index_array[surf_state->bg_index], tmp_array, &dev_descrip.col_mode);
    surf_state->color_table[0].r = (int)(tmp_array[0][0] * 255.0f);
    surf_state->color_table[0].g = (int)(tmp_array[0][1] * 255.0f);
    surf_state->color_table[0].b = (int)(tmp_array[0][2] * 255.0f);
    vdiqco(&one, &dev_descrip.index_array[surf_state->fg_index], tmp_array, &dev_descrip.col_mode);
    surf_state->color_table[1].r = (int)(tmp_array[0][0] * 255.0f);
    surf_state->color_table[1].g = (int)(tmp_array[0][1] * 255.0f);
    surf_state->color_table[1].b = (int)(tmp_array[0][2] * 255.0f);
  } /* end if vector SVDI */

  else { /* raster SVDI */

    vdiqci(&dev_descrip.att_array[11], &dev_descrip.att_array[12], &dev_descrip.att_array[13],
           &surf_state->bg_index);
    vdiqci(&dev_descrip.att_array[8], &dev_descrip.att_array[9], &dev_descrip.att_array[10],
           &surf_state->fg_index);

    /* set the current surface color table */
    surf_state->color_table[0].r = (int)(dev_descrip.att_array[8] * 255.0f);
    surf_state->color_table[0].g = (int)(dev_descrip.att_array[9] * 255.0f);
    surf_state->color_table[0].b = (int)(dev_descrip.att_array[10] * 255.0f);
    surf_state->color_table[1].r = (int)(dev_descrip.att_array[11] * 255.0f);
    surf_state->color_table[1].g = (int)(dev_descrip.att_array[12] * 255.0f);
    surf_state->color_table[1].b = (int)(dev_descrip.att_array[13] * 255.0f);
  } /* end else raster SVDI */

  /* set/reset the SVDI attribute array */
  /* --set everything to -1 to force resetting of attributes */
  surf_state->vdi_attrib.fg_color   = -1;
  surf_state->vdi_attrib.bg_color   = -1;
  surf_state->vdi_attrib.intensity  = -1.0;
  surf_state->vdi_attrib.line_type  = -1;
  surf_state->vdi_attrib.line_width = -1.0;
  surf_state->vdi_attrib.char_x     = -1.0;
  surf_state->vdi_attrib.char_y     = -1.0;
  surf_state->vdi_attrib.fg_rgb[0]  = -1.0;
  surf_state->vdi_attrib.fg_rgb[1]  = -1.0;
  surf_state->vdi_attrib.fg_rgb[2]  = -1.0;
  surf_state->vdi_attrib.bg_rgb[0]  = -1.0;
  surf_state->vdi_attrib.bg_rgb[1]  = -1.0;
  surf_state->vdi_attrib.bg_rgb[2]  = -1.0;

  /* set color_set flag to false = color table not set by user */
  surf_state->color_set = FALSE;

  /* set picture state to clean */
  surf_state->pic_dirty = CCLEAN;

} /* end init_state */

/* set_dev_descrip */
/* Set the device description table */
static void set_dev_descrip(void)
{
  /* dev_descrip is global */
  int   j;         /* loop index */
  int   qdc_index; /* used for inquiring SVDI */
  float value;     /* value returned from inquire SVDI */

  /* set copy_class */
  qdc_index = 1;
  vdiqdc(&qdc_index, &value);
  if (value == 0.) {
    dev_descrip.copy_class = CHARD;
  }
  else {
    dev_descrip.copy_class = CSOFT;
  }

  /* device corners in device coords */
  dev_descrip.dc1.x = 0;
  dev_descrip.dc1.y = 0;
  qdc_index         = 15;
  vdiqdc(&qdc_index, &value);
  dev_descrip.dc2.x = (int)value;
  qdc_index         = 16;
  vdiqdc(&qdc_index, &value);
  dev_descrip.dc2.y = (int)value;

  /* vector or raster SVDI (0=vector, 1=raster) */
  qdc_index = 29;
  vdiqdc(&qdc_index, &value);
  dev_descrip.svdi_type = (int)value;

  /* vector SVDI only does indexed color, raster does both */
  if (dev_descrip.svdi_type == 0) {
    dev_descrip.csm_avail = CCLRI;
  }
  else {
    dev_descrip.csm_avail = CCLRID;
  }

  /* device class */
  qdc_index = 13;
  vdiqdc(&qdc_index, &value);
  dev_descrip.dev_class = (value == 0.) ? COUTPT : COUTIN;

  /* NDC space */
  vdiqnd(&dev_descrip.xndc_max, &dev_descrip.yndc_max);

  /* output status */
  vdiqos(&dev_descrip.att_array[1]);

  /* min line width  - SVDI returns value in DC */
  qdc_index = 19;
  vdiqdc(&qdc_index, &value);
  dev_descrip.linewidth_min = (int)value;

  /* max line width  - SVDI returns value in DC */
  qdc_index = 31;
  vdiqdc(&qdc_index, &value);
  dev_descrip.linewidth_max = (int)value;

  /* use the default SVDI linewidth as nominal */
  /* --convert to device coordinates */
  value = (dev_descrip.att_array[5] / dev_descrip.xndc_max) * .01 *
          (dev_descrip.dc2.x - dev_descrip.dc1.x);
  dev_descrip.linewidth_nominal = (int)max(value, 1.);

  /* --make sure that for some reason it isn't out of range */
  dev_descrip.linewidth_nominal =
      max(min(dev_descrip.linewidth_nominal, dev_descrip.linewidth_max), dev_descrip.linewidth_min);

  /* nominal,min,max marker size - only one size */
  qdc_index = 20;
  vdiqdc(&qdc_index, &value);
  dev_descrip.mark_nominal = (int)value;
  dev_descrip.mark_min     = (int)value;
  dev_descrip.mark_max     = (int)value;

  /* inquire the number of colors */
  qdc_index = 4;
  vdiqdc(&qdc_index, &value);

  /* check for max colors */
  dev_descrip.num_cols = ((int)value > COLOR_TABLE_SIZE) ? COLOR_TABLE_SIZE : (int)value;

  /* set up a color index array for use in color table stuff */
  for (j = 0; j < 256; j++) {
    dev_descrip.index_array[j] = j;
  }

  /* set SVDI color mode */
  dev_descrip.col_mode = 0;

  /* inquire default SVDI color table.  */
  vdiqco(&dev_descrip.num_cols, dev_descrip.index_array, dev_descrip.color_array,
         &dev_descrip.col_mode);

} /* end set_dev_descrip */

/* reset_vdi */
/* Reset SVDI values. Used if CI is called twice */
static void reset_vdi(surf_statelist *surf_state)
/* current surface statelist */
{

  /* reset the output status */
  vdstos(&dev_descrip.att_array[1]);

  /* only reset color table is user has changed it */
  if (surf_state->color_set) {
    vdstco(&dev_descrip.num_cols, dev_descrip.index_array, dev_descrip.color_array,
           &dev_descrip.col_mode);
  }
} /* end reset_vdi */

/* set_mapping */
/* Set up the VDC to NDC mappings. Also sets the effective viewports
 * in NDC and VDC.
 */
static void set_mapping(surf_statelist *surf_state)
{
  point ndc_vp1, ndc_vp2;

  /* viewport in NDC space */
  ndc_vp1.x = cur_state->vp1.x * dev_descrip.xndc_max;
  ndc_vp1.y = cur_state->vp1.y * dev_descrip.yndc_max;
  ndc_vp2.x = cur_state->vp2.x * dev_descrip.xndc_max;
  ndc_vp2.y = cur_state->vp2.y * dev_descrip.yndc_max;

  /* x and y scale maps VDC to viewport in NDC (ISOTROPIC Mapping) */
  /* ...smallest scale magnitude, but keep the original sign */
  surf_state->xscale = (ndc_vp2.x - ndc_vp1.x) / (surf_state->vdc2.x - surf_state->vdc1.x);
  surf_state->yscale = (ndc_vp2.y - ndc_vp1.y) / (surf_state->vdc2.y - surf_state->vdc1.y);

  if (fabs(surf_state->xscale) < fabs(surf_state->yscale)) {
    surf_state->yscale =
        (surf_state->yscale < 0) ? -1 * fabs(surf_state->xscale) : fabs(surf_state->xscale);
  }
  else {
    surf_state->xscale =
        (surf_state->xscale < 0) ? -1 * fabs(surf_state->yscale) : fabs(surf_state->yscale);
  }

  /* x and y offset for VDC to NDC mapping */
  surf_state->xoffset = (ndc_vp1.x < ndc_vp2.x)
                            ? ndc_vp1.x - surf_state->vdc1.x * surf_state->xscale
                            : ndc_vp2.x - surf_state->vdc2.x * surf_state->xscale;

  surf_state->yoffset = (ndc_vp1.y < ndc_vp2.y)
                            ? ndc_vp1.y - surf_state->vdc1.y * surf_state->yscale
                            : ndc_vp2.y - surf_state->vdc2.y * surf_state->yscale;

  /* effective viewport */
  /* - mode is FRACTION OF DRAWING SURFACE */
  surf_state->eff_vp1.x = cur_state->vp1.x;
  surf_state->eff_vp1.y = cur_state->vp1.y;
  surf_state->eff_vp2.x =
      (surf_state->vdc2.x * surf_state->xscale + surf_state->xoffset) / dev_descrip.xndc_max;
  surf_state->eff_vp2.y =
      (surf_state->vdc2.y * surf_state->yscale + surf_state->yoffset) / dev_descrip.yndc_max;

} /* end set_mapping */

/* set_clipping */
/* Set up the clip region. */
static void set_clipping(surf_statelist *my_cur_state)
{
  point clip1, clip2; /* temp clip values */
  clip1.x = clip1.y = clip2.x = clip2.y = 0;

  /* The clip region depends on clip indicator and drawing surface
   * clip indicator.
   */
  switch (my_cur_state->clip_indicator) {
  case CON:

    /* VDC clipping is on */
    switch (my_cur_state->ds_clip_indicator) {

    case CDCOFF: /* view surface clipping off */

      /* map the effective clip rectangle to NDC */
      my_cur_state->clip_on = TRUE;
      clip1.x = my_cur_state->eff_clip_rect1.x * my_cur_state->xscale + my_cur_state->xoffset;
      clip1.y = my_cur_state->eff_clip_rect1.y * my_cur_state->yscale + my_cur_state->yoffset;
      clip2.x = my_cur_state->eff_clip_rect2.x * my_cur_state->xscale + my_cur_state->xoffset;
      clip2.y = my_cur_state->eff_clip_rect2.y * my_cur_state->yscale + my_cur_state->yoffset;

      break; /* end case CDCOFF */

    case CVPORT: /* clip at viewport */

      /* map the effective clip rectangle to NDC */
      my_cur_state->clip_on = TRUE;

      clip1.x = my_cur_state->eff_clip_rect1.x * my_cur_state->xscale + my_cur_state->xoffset;
      clip1.y = my_cur_state->eff_clip_rect1.y * my_cur_state->yscale + my_cur_state->yoffset;
      clip2.x = my_cur_state->eff_clip_rect2.x * my_cur_state->xscale + my_cur_state->xoffset;
      clip2.y = my_cur_state->eff_clip_rect2.y * my_cur_state->yscale + my_cur_state->yoffset;

      break; /* end case CVPORT */

    case CDCREC: /* clip at display surface */

      /* map the effective clip rectangle to NDC and intersect it
         with the max NDC space */
      my_cur_state->clip_on = TRUE;

      clip1.x = my_cur_state->eff_clip_rect1.x * my_cur_state->xscale + my_cur_state->xoffset;
      clip1.y = my_cur_state->eff_clip_rect1.y * my_cur_state->yscale + my_cur_state->yoffset;
      clip2.x = my_cur_state->eff_clip_rect2.x * my_cur_state->xscale + my_cur_state->xoffset;
      clip2.y = my_cur_state->eff_clip_rect2.y * my_cur_state->yscale + my_cur_state->yoffset;

      clip1.x = max(min(clip1.x, clip2.x), 0.0);
      clip1.y = max(min(clip1.y, clip2.y), 0.0);
      clip2.x = min(max(clip1.x, clip2.x), dev_descrip.xndc_max);
      clip2.y = min(max(clip1.y, clip2.y), dev_descrip.yndc_max);

      break; /* end case CDCREC */

    default:
      /* illegal enumeration type - ignore function */
      break;

    } /* end switch ds_clip_indicator  */

    break; /* end case CON */

  case COFF:

    /* VDC clipping is off */
    switch (my_cur_state->ds_clip_indicator) {

    case CDCOFF: /* display surface clipping off */

      /* no clip */
      my_cur_state->clip_on = FALSE;
      break; /* end case CDCOFF */

    case CVPORT: /* clip at viewport */

      /* clip at NDC effective viewport  */
      my_cur_state->clip_on = TRUE;
      clip1.x               = my_cur_state->eff_vp1.x * dev_descrip.xndc_max;
      clip1.y               = my_cur_state->eff_vp1.y * dev_descrip.yndc_max;
      clip2.x               = my_cur_state->eff_vp2.x * dev_descrip.xndc_max;
      clip2.y               = my_cur_state->eff_vp2.y * dev_descrip.yndc_max;
      break; /* end case CVPORT */

    case CDCREC: /* clip at display surface */

      /* clip at max NDC space */
      my_cur_state->clip_on = TRUE;
      clip1.x               = 0.0;
      clip1.y               = 0.0;
      clip2.x               = dev_descrip.xndc_max;
      clip2.y               = dev_descrip.yndc_max;
      break; /* end case CDCREC */

    default:
      /* illegal enumeration type - ignore function */
      break;

    } /* end switch ds_clip_indicator  */

    break; /* end case COFF */

  default:
    /* illegal enumeration type - ignore function */
    break;

  } /* end switch clip_indicator  */

  /* set clipmin and clipmax such that clipmin is the minimum corner of
   * the clipping window, and clipmax is the maximum. This is done
   * so that the clipping algorithms work correctly with mirroring.
   */
  if (my_cur_state->clip_on) {
    my_cur_state->clipmin.x = min(clip1.x, clip2.x);
    my_cur_state->clipmax.x = max(clip1.x, clip2.x);
    my_cur_state->clipmin.y = min(clip1.y, clip2.y);
    my_cur_state->clipmax.y = max(clip1.y, clip2.y);
  }
} /* end set_clipping */

/* set_foreground_color */
/* Set the SVDI foreground color */
static void set_foreground_color(surf_statelist *surf_state, int *colors)
{
  int zero = 0;
  int one  = 1;

  /* if color selection mode is in error, ignore function */
  if (surf_state->csm == -1) {
    return;
  }

  /* indexed or direct color? */
  /* --vector SVDI only does indexed color */
  if (surf_state->csm == CINDEX || dev_descrip.svdi_type == 0) {

    /* does foreground need to be updated? */
    if (colors[0] != cur_state->vdi_attrib.fg_color) {

      /* color index mappings:
       *   if set index 0  --> set bg_index
       *   if set index 1  --> set fg_index
       *   if set bg_index --> set 0
       *   if set fg_index --> set 1
       */
      if (colors[0] == 0) {
        vdstfc(&cur_state->bg_index);
      }
      else if (colors[0] == 1) {
        vdstfc(&cur_state->fg_index);
      }
      else if (colors[0] == cur_state->bg_index) {
        vdstfc(&zero);
      }
      else if (colors[0] == cur_state->fg_index) {
        vdstfc(&one);
      }
      else {
        vdstfc(&colors[0]);
      }

      /* update att_array */
      cur_state->vdi_attrib.fg_color = colors[0];

      /* flag (direct) that foreground color has changed */
      /* do i need to do this? */
      cur_state->vdi_attrib.fg_rgb[0] = -1.0;
      cur_state->vdi_attrib.fg_rgb[1] = -1.0;
      cur_state->vdi_attrib.fg_rgb[2] = -1.0;
    } /* end does foreground... */
  } /* end indexed color */

  else /* direct color */

    /* does foreground need to be updated? */
    /* -- i need to check this out - might need to store as int */
    if (colors[0] != (int)(cur_state->vdi_attrib.fg_rgb[0] * 255.0f) ||
        colors[1] != (int)(cur_state->vdi_attrib.fg_rgb[1] * 255.0f) ||
        colors[2] != (int)(cur_state->vdi_attrib.fg_rgb[2] * 255.0f)) {

      /* update att_array */
      cur_state->vdi_attrib.fg_rgb[0] = (float)colors[0] / 255.0f;
      cur_state->vdi_attrib.fg_rgb[1] = (float)colors[1] / 255.0f;
      cur_state->vdi_attrib.fg_rgb[2] = (float)colors[2] / 255.0f;

      /* set new foreground color */
      vdfrgb(&cur_state->vdi_attrib.fg_rgb[0], &cur_state->vdi_attrib.fg_rgb[1],
             &cur_state->vdi_attrib.fg_rgb[2]);

      /* flag (indexed) that foreground color has changed */
      /* do i need to do this? */
      cur_state->vdi_attrib.fg_color = -1;
    } /* end does foreground... */

} /* end set_foreground_color */

/* set_background_color */
/* Set the SVDI background color */
static void set_background_color(surf_statelist *surf_state, int *colors)
{
  int   i;                      /* loop index */
  int   index = 0;              /* color index */
  float dr, dg, db, dmin, dist; /* for finding the closet index */
  int   one     = 1;
  float epsilon = .001;

  /* does background need to be updated */
  /* --background color is saved in att_array, even for vector */
  if (colors[0] != (int)(cur_state->vdi_attrib.bg_rgb[0] * 255.0f) ||
      colors[1] != (int)(cur_state->vdi_attrib.bg_rgb[1] * 255.0f) ||
      colors[2] != (int)(cur_state->vdi_attrib.bg_rgb[2] * 255.0f)) {

    /* store new values in att_array */
    cur_state->vdi_attrib.bg_rgb[0] = (float)colors[0] / 255.0f;
    cur_state->vdi_attrib.bg_rgb[1] = (float)colors[1] / 255.0f;
    cur_state->vdi_attrib.bg_rgb[2] = (float)colors[2] / 255.0f;

    /* set the cgi state color table - index 0 */
    cur_state->color_table[0].r = colors[0];
    cur_state->color_table[0].g = colors[1];
    cur_state->color_table[0].b = colors[2];

    /* has the SVDI color table ever been set */
    if (!surf_state->color_set) { /* has never been set */

      /* find out if default SVDI table contains this color */
      if (dev_descrip.svdi_type == 1) { /* raster SVDI */

        vdiqci(&cur_state->vdi_attrib.bg_rgb[0], &cur_state->vdi_attrib.bg_rgb[1],
               &cur_state->vdi_attrib.bg_rgb[2], &i);

        /* does it match */
        dr   = cur_state->vdi_attrib.bg_rgb[0] - dev_descrip.color_array[i][0];
        dg   = cur_state->vdi_attrib.bg_rgb[1] - dev_descrip.color_array[i][1];
        db   = cur_state->vdi_attrib.bg_rgb[2] - dev_descrip.color_array[i][2];
        dmin = dr * dr + dg * dg + db * db;
      } /* end raster SVDI */

      else { /* vector SVDI */

        for (i = 0; i < 8; i++) {
          dr   = cur_state->vdi_attrib.bg_rgb[0] - dev_descrip.color_array[i][0];
          dg   = cur_state->vdi_attrib.bg_rgb[1] - dev_descrip.color_array[i][1];
          db   = cur_state->vdi_attrib.bg_rgb[2] - dev_descrip.color_array[i][2];
          dist = dr * dr + dg * dg + db * db;

          if (i == 0 || dist < dmin) {
            dmin  = dist;
            index = i;
          }
        } /* end for i */
      } /* end vector SVDI */

      /* is it close enough? */
      if (dmin <= epsilon) {
        cur_state->bg_index = index;
      }
      else { /* not close enough, set SVDI color table */

        cur_state->bg_index = 0;
        vdstco(&one, &dev_descrip.index_array[cur_state->bg_index], &cur_state->vdi_attrib.bg_rgb,
               &dev_descrip.col_mode);
        cur_state->color_set = TRUE; /* flag that CT has been set */
      } /* end else not close enough */

    } /* end has never been set */

    else { /* color table has been set */

      /* set SVDI color table */
      vdstco(&one, &dev_descrip.index_array[cur_state->bg_index], &cur_state->vdi_attrib.bg_rgb,
             &dev_descrip.col_mode);

    } /* end else has been set */

    /* set the SVDI background color */
    if (dev_descrip.svdi_type == 1) { /* raster SVDI */
      vdbrgb(&cur_state->vdi_attrib.bg_rgb[0], &cur_state->vdi_attrib.bg_rgb[1],
             &cur_state->vdi_attrib.bg_rgb[2]);
    }
    else /* vector SVDI */

      /* once the ct has been set, bg_index never changes */
      if (!cur_state->color_set) {
        vdstbc(&cur_state->bg_index);
      }

  } /* end does background need to be updated */
} /* end set_background_color */

/* report_error */
/* Stuff the error info into the error queue */
static void report_error(surf_statelist *surf_state, int e_class, int e_num, int f_id)
{
  int err_slot; /* slot in the queue to put the error */

  /* is error reporting on for this class? */
  if (surf_state->err_flag[e_class - 1] == CEHON) {

    /* figure out where the next slot is */
    err_slot = surf_state->err_head_ptr + surf_state->err_count;
    if (err_slot >= ERROR_LIST_SIZE) {
      err_slot = err_slot - ERROR_LIST_SIZE;
    }

    /* store errors if there is room in the queue */
    if (surf_state->err_count < ERROR_LIST_SIZE) {
      surf_state->err_queue[err_slot].err_class = e_class;
      surf_state->err_queue[err_slot].err_num   = e_num;
      surf_state->err_queue[err_slot].func_id   = f_id;
      surf_state->err_count++;
    } /* end if room in the queue */

    else { /* not enuff room in the queue */

      /* if the most recent error is ERROR REPORTS LOST, increase count
         of errors lost by one (count is stored in func_id) */
      if (surf_state->err_queue[err_slot - 1].err_class == -1 &&
          surf_state->err_queue[err_slot - 1].err_num == 0) {
        surf_state->err_queue[err_slot - 1].func_id++;
      }
      else { /* create an ERROR REPORTS LOST error */
        /* The format for an ERROR REPORTS LOST has not yet
         * been defined.  For now, use error class -1, error
         * number 0, and the count in function id.
         * (as defined in SRCP document)
         */
        surf_state->err_queue[err_slot - 1].err_class = -1;
        surf_state->err_queue[err_slot - 1].err_num   = 0;
        surf_state->err_queue[err_slot - 1].func_id   = 2;
      } /* end else create */
    } /* end not enuff room in the queue */

  } /* end if error reporting is on */
} /* end report_error */

/* gettoken */
/* Return one token from data record....written by Pat MaGee, LANL */
static void gettoken(int *index_p, int *numrecs_p, char *data_p, int max_chars, char *outtoken_p)
/* in: 1st char of whitespace before token */
/* out: 1st char past token */
/* in: number of data records */
/* in: 80 character records */
/* in: max # of chars to put into outtoken */
/* out: where to put output token */
{
  char c;         /* char from data record */
  int  num_chars; /* # chars put into outtoken so far */
  int  workstate; /* 0=whitespace,1=nonwhite,2=endofstring */

  num_chars = 0;

  workstate = 0;
  if (*index_p >= (*numrecs_p * 80)) {
    workstate = 2;
  }

  while (workstate == 0) { /* while scanning whitespace */
    c = data_p[*index_p];
    if ((c == ' ') || (c == '\t') || (c == '\r')) {
      /* c is whitespace, keep scanning */
      ++(*index_p);
      if (*index_p >= *numrecs_p * 80) {
        workstate = 2;
      } /* end if off end of data */
    }
    else {
      /* c is not whitespace, is start of token */
      workstate = 1;
    } /* end if c is whitespace */
  } /* end while skipping whitespace */

  while (workstate == 1) { /* while in token */
    c = data_p[*index_p];
    if (num_chars >= max_chars) {
      workstate = 3;
    }
    if ((c == ' ') || (c == '\t') || (c == '\r')) {
      workstate = 0;
    }
    else { /* c is not whitespace, output and keep scanning */
      outtoken_p[num_chars] = c;
      ++num_chars;
      ++(*index_p);
      if (*index_p >= *numrecs_p * 80) {
        workstate = 2;
      } /* end if off end of data */
    } /* end if c is whitespace */

  } /* end while putting chars into outtoken */

  outtoken_p[num_chars] = 0;

} /* end gettoken */

/* poly_clip */
/* Clip a polygon. Based on Sutherland - Hodgman algorithm
 * Returns the new polygon vertices, and the new vertex count.
 * One day, I need to make this more efficient by moving everything
 * inline.
 */
static int poly_clip(point *cmin, point *cmax, float *vx, float *vy, int vlen, float *xout,
                     float *yout, int *lenout)
{
  point s, p, t;
  int   i, j, bnd, curlen;
  float xtemp[MAX_POLYGON];
  float ytemp[MAX_POLYGON];

  /* make sure there is something to clip */
  if (vlen <= 0) {
    return (FALSE);
  }

  for (bnd = 0; bnd < 4; bnd++) { /* loop through all boundaries */

    /* start with last vertex */
    if (bnd == 0) { /* 1st time through use original vertex list */
      curlen  = vlen;
      *lenout = 0;
      s.x     = vx[curlen - 1];
      s.y     = vy[curlen - 1];
    }
    else { /* after 1st time through, use new vertex list */
      curlen  = *lenout;
      *lenout = 0;
      if (curlen > 0) {
        s.x = xout[curlen - 1];
        s.y = yout[curlen - 1];
      }
    }

    for (i = 0; i < curlen; i++) { /* loop through all vertices */

      if (bnd == 0) {
        p.x = vx[i];
        p.y = vy[i];
      }
      else {
        p.x = xout[i];
        p.y = yout[i];
      }

      /* is vertex p 'inside' current boundary */
      if (inside_bnd(&p, cmin, cmax, bnd)) {

        /* is vertex s 'inside' current boundary */
        if (inside_bnd(&s, cmin, cmax, bnd)) {

          /* both p and s are inside, add p to vertex list */
          xtemp[*lenout]     = p.x;
          ytemp[(*lenout)++] = p.y;
        } /* end s is inside */

        else { /* p is inside, s is outside, intersect */

          t = intersect_bnd(&s, &p, cmin, cmax, bnd);

          /* add t and p to the vertex list */
          xtemp[*lenout]     = t.x;
          ytemp[(*lenout)++] = t.y;
          xtemp[*lenout]     = p.x;
          ytemp[(*lenout)++] = p.y;
        } /* end else intersect */

      } /* end if p is inside */

      else { /* p is outside, is s inside */

        if (inside_bnd(&s, cmin, cmax, bnd)) {

          /* p is outside, s is inside, intersect */
          t = intersect_bnd(&s, &p, cmin, cmax, bnd);

          /* add t to the vertex list */
          xtemp[*lenout]     = t.x;
          ytemp[(*lenout)++] = t.y;
        } /* end intersect */
      } /* end else p outside, s inside */

      s.x = p.x;
      s.y = p.y;
    } /* end for i */

    for (j = 0; j < *lenout; j++) {
      xout[j] = xtemp[j];
      yout[j] = ytemp[j];
    } /* end for j */

  } /* end for bnd */

  if (*lenout <= 0) {
    return (FALSE);
  }
  return (TRUE);
} /* end poly_clip */

/* inside_bnd */
/* Used by poly_clip and returns true if the vertex is within the
 * boundary. ( where the boundaries are numbered 1-4,corresponding to
 * top, right, bottom, left )
 */
static int inside_bnd(point *v, point *bmin, point *bmax, int bound_num)
{
  /* this routine assumes a rectangular boundary */
  switch (bound_num) {
  case 0: return (v->y <= bmax->y); break; /* top */
  case 1: return (v->x <= bmax->x); break; /* right */
  case 2: return (v->y >= bmin->y); break; /* bottom */
  case 3: return (v->x >= bmin->x); break; /* left */
  default: return (FALSE); break;
  } /* end switch */
} /* end inside_bnd */

/* intersect_bnd */
/* Used by poly_clip. Returns the the intersection point of a polygon
 * edge with a clip boundary.  The clip boundary is coded 1-4, corresponding
 * to top, right, bottom, and left boundary edges, respectively.
 */
static point intersect_bnd(point *p1, point *p2, point *bmin, point *bmax, int bound_num)
{
  point temp;
  temp.x = 0;
  temp.y = 0;

  switch (bound_num) {
  case 0: /* top */
    temp.x = p1->x + (p2->x - p1->x) * (bmax->y - p1->y) / (p2->y - p1->y);
    temp.y = bmax->y;
    break;

  case 1: /* right */
    temp.y = p1->y + (p2->y - p1->y) * (bmax->x - p1->x) / (p2->x - p1->x);
    temp.x = bmax->x;
    break;

  case 2: /* bottom */
    temp.x = p1->x + (p2->x - p1->x) * (bmin->y - p1->y) / (p2->y - p1->y);
    temp.y = bmin->y;
    break;

  case 3: /* left */
    temp.y = p1->y + (p2->y - p1->y) * (bmin->y - p1->x) / (p2->x - p1->x);
    temp.x = bmin->x;
    break;

  default:
    /* error check ? */
    break;

  } /* end switch */
  return temp;
} /* end intersect_bnd */

/*-------------------------------------------------------------*/
/* >> CDR ROUTINES                                             */
/*-------------------------------------------------------------*/

/* cur_state is global - points to current statelist */
/*
  Open file for sequential access.  FORTRAN unit number is ignored.
  File descriptor which is returned from open is stored in current
  statelist
*/
void cdrofs(int *ifilcd)
{
  int        errnum, errsev;
  char       symbol[1024];
  int        qdc_index;
  float      value;
  char      *devid;
  char      *env;
  static int file_cnt = 1;

  /* if user hasn't named the file, build a default */
  if (!cur_state->filename[0]) { /* filename is blank */
    /* get the SVDI device code */
    qdc_index = 23;
    vdiqdc(&qdc_index, &value);
    devid = get_devid_char(value);
    if (devid != NULL) {
      snprintf(cur_state->filename, 100, "cgi%s%d", devid, file_cnt++);
    }
    else {
      snprintf(cur_state->filename, 100, "cgiout%d", file_cnt++);
    }
  }

  /* copy filename to symbol */
  copy_string(symbol, cur_state->filename, 1024);

  /* check the environment to see if a file name has been assigned */
  env = getenv(symbol);
  if (env != 0 && strlen(env) < 1024) {
    snprintf(symbol, 1024, "%s", env);
  }

  /* open the file  - if it doesn't exist, create it with mode 664 */
  /* -- open returns a file descriptor which is stored in the statelist */
  /* -- O_TRUNC ??? read/writing??? */
  if ((cur_state->file_d = open(symbol, (O_CREAT | O_TRUNC | O_RDWR), 0664)) == -1) {
    errnum = 722;
    errsev = 10;
    char err[50];
    snprintf(err, 50, "SVDI ERROR NUMBER %d SEVERITY CODE %d", errnum, errsev);
    perror(err);
  }
}

/* this routine is here only for compatibility with SVDI drivers */
void cdrof3(int *ifilcd, int *idum) { cdrofs(ifilcd); }

/* this routine is here only for compatibility with SVDI drivers */
void cdroff(int *ifilcd, int *idum1, int *idum2, int *idum3) { cdrofs(ifilcd); }

/*
  Try to read LENGTH words from the file. Return the actual number
  of words read. FORTRAN unit number is ignored, file descriptor is
  gotten from the current statelist
*/
void cdrrfs(int *ifilcd, int *length, char *buffer, int *istat)
{

  /* if the file hasn't been opened, open it */
  if (cur_state->file_d == -1) {
    cdrofs(ifilcd);
  }

  /* read the data as a byte stream */
  if ((*istat = read(cur_state->file_d, buffer, (*length * cdrcom.KCPW))) != -1) {
    *istat = *istat / cdrcom.KCPW;
  }
}

/*
  Write LENGTH words from BUFFER to the file.  FORTRAN unit number is
  ignored, file descriptor is gotten from the current state list
*/
void cdrwfs(int *ifilcd, int *length, char *buffer, int *istat)
{

  /* if the file hasn't been opened, open it */
  if (cur_state->file_d == 0) {
    cdrofs(ifilcd);
  }

  /* write out the data as a byte stream. */
  if ((*istat = write(cur_state->file_d, buffer, (*length * cdrcom.KCPW))) == -1) {
    *istat = 1;
  }
  else {
    *istat = 0;
  }
}

/*
  Close file sequential. FORTRAN unit number is ignore and file
  descriptor is gotten from current state
*/
void cdrcfs(int *ifilcd, int *eof)
{
  int  istat;
  char buf[4];

  /* if eof = 1 then write eof on file */
  if (*eof == 1) {
    *buf = EOF;
    if ((istat = write(cur_state->file_d, buf, 4)) == -1) {
      char err[50];
      snprintf(err, 50, "%s", "CDRCFS error");
      perror(err);
    }
  }

  /* close the file */
  close(cur_state->file_d);
}

/* open routine for the Abekas */
void cdroab(int *ifilcd, int *frame)
{
  char ic[5];
  int  i, inbr;

  /* form the file name from the frame number */
  inbr = *frame;
  for (i = 3; i >= 0; i--) {
    ic[i] = inbr % 10 + '0';
    inbr  = inbr / 10;
  }
  ic[4] = '\0';

  /* set the file name in the state list */
  snprintf(cur_state->filename, 100, "%s.RGB", ic);

  cdrofs(ifilcd);
}

/*
  Metafile buffering routine. (device dependent )
  Metafile driver passes OUTARY as an integer array, always with
  16 bits/word, right justified.  This buffering routine assumes
  16 bits/word and assumes 8 bit bytes.
*/
void nmtbuf(int *numwds, unsigned outary[])
{
  static unsigned mask1 = ~(~0u << 8) << 8; /* mask off higher 8 bits */
  static unsigned mask2 = ~(~0u << 8);      /* mask off lower 8 bits */

  int i;     /* loop variable */
  int istat; /* error reporting */

  /* cur_state is global and points to the current state */

  /* check for buffer flush and if there is something to flush */
  if (*numwds <= 0 && cur_state->buff_ptr > 0) {

    /* write out the data as a byte stream. */
    if ((istat = write(cur_state->file_d, cur_state->buffer, cur_state->buff_ptr)) == -1) {
      char err[50]; /* for error reporting */
      snprintf(err, 50, "%s", "NMTBUF: write error");
      perror(err);
    } /* end write buffer */

    /* reset buff pointer */
    cur_state->buff_ptr = 0;

  } /* end if buffer flush */

  else { /* pack outary into buffer */

    for (i = 0; i < *numwds; i++) {
      cur_state->buffer[cur_state->buff_ptr++] = (outary[i] & mask1) >> 8;
      cur_state->buffer[cur_state->buff_ptr++] = (outary[i] & mask2);

      if (cur_state->buff_ptr > BUFFER_SIZE - 2) {

        /* write out the data as a byte stream. */
        if ((istat = write(cur_state->file_d, cur_state->buffer, cur_state->buff_ptr)) == -1) {
          char err[50]; /* for error reporting */
          snprintf(err, 50, "%s", "NMTBUF: write error");
          perror(err);

        } /* end write buffer */

        /* reset buff pointer */
        cur_state->buff_ptr = 0;

      } /* end if > BUFFER_SIZE */
    } /* end for i=... */
  }
} /* end nmtbuf */
