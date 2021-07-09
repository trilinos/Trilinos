/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef SDCGI_H
#define SDCGI_H

#include "fortyp.h"

void cesc2_(f_integer *funcid, f_integer *n, f_real *args);
void cqchh2_(f_integer *txp, f_integer *nreq, f_integer *first, f_integer *vstat, f_integer *ntotal,
             f_integer *nlist, f_integer *chhit);
void ctx2_(f_real *x, f_real *y, f_integer *text1, f_integer *length);
void cgtxx2_(f_real *x, f_real *y, f_integer *vstat, f_integer *vconc, f_real *f_xconc,
             f_real *yconc, f_real *x1, f_real *y1, f_real *x2, f_real *y2, f_real *x3, f_real *y3,
             f_real *x4, f_real *y4);

void ci_(f_integer *pds);
void ct_(void);
void cxdfac_(void);
void cpds_(f_integer *clear);
void cendpg_(void);
void cbc_(f_integer *red, f_integer *green, f_integer *blue);
void cvdcx_(f_real *x1, f_real *y1, f_real *x2, f_real *y2);
void cv_(f_real *x1, f_real *y1, f_real *x2, f_real *y2);
void ccl_(f_integer *clipi);
void cdscl_(f_integer *clipi);
void cdqerr_(f_integer *nreq, f_integer *vstat, f_integer *nrem, f_integer *nret, f_integer *errcl,
             f_integer *errnm, f_integer *funcid);
void cerhct_(f_integer *n, f_integer *erclas, f_integer *hflag);
void ccixp_(f_integer *cip);
void cesc1_(f_integer *funcid, f_integer *ldr, char *data, f_integer *drec_size);
void cqid_(f_integer *maxchr, f_integer *vstat, f_integer *dclass, char *devid);
void cqd_(f_integer *vstat, f_integer *hscopy, f_integer *disp, f_integer *bcolor, f_integer *dynbc,
          f_integer *dynvdm, f_integer *dx1, f_integer *dy1, f_integer *dx2, f_integer *dy2,
          f_real *width, f_real *height, f_integer *pixloc);
void clf_(f_integer *n, f_integer *funccl, f_integer *funcid, f_integer *vstat, f_integer *supprt);
void clpr_(f_integer *n, char *profid, f_integer *profid_size, f_integer *vstat, f_integer *supprt);
void cqsp_(f_integer *vstat, f_integer *nvip, f_integer *vip, f_integer *nvrp, f_integer *vrfmt,
           f_integer *vrexp, f_integer *vrfrac, f_integer *nip, f_integer *ip, f_integer *nrp,
           f_integer *rfmt, f_integer *rexp, f_integer *rfrac, f_integer *nixp, f_integer *ixp,
           f_integer *ncp, f_integer *cp, f_integer *ncixp, f_integer *cixp);
void clesc_(f_integer *n, f_integer *escid, f_integer *vstat, f_integer *supprt);
void cqp_(f_integer *vstat, f_integer *vip, f_integer *vrfmt, f_integer *vrexp, f_integer *vrfrac,
          f_integer *ip, f_integer *rfmt, f_integer *rexp, f_integer *rfrac, f_integer *ixp,
          f_integer *cp, f_integer *cixp);
void cqcl_(f_integer *vstat, f_integer *clip1, f_integer *clipr, f_integer *sclip1,
           f_integer *sclipr);
void cpl_(f_integer *np, f_real *px, f_real *py);
void cdjpl_(f_integer *np, f_real *px, f_real *py);
void cpm_(f_integer *np, f_real *px, f_real *py);
void ctx1_(f_real *x, f_real *y, f_integer *flag, char *text, f_integer *text_size);
void cpg_(f_integer *np, f_real *px, f_real *py);
void cca_(f_real *x1, f_real *y1, f_real *x2, f_real *y2, f_real *x3, f_real *y3, f_integer *nx,
          f_integer *ny, f_integer *lcp, f_integer *cells);
void clnt_(f_integer *lntyp);
void clnw_(f_real *lnwid);
void clnc_(f_integer *lnclr);
void cmkt_(f_integer *mktyp);
void cmkc_(f_integer *mkclr);
void ctxp_(f_integer *txp);
void ctxc_(f_integer *txclr);
void cchh_(f_real *chhit);
void ccho_(f_integer *xup, f_integer *yup, f_integer *xbase, f_integer *ybase);
void cis_(f_integer *instyl);
void cflc_(f_integer *fclr);
void ccsm_(f_integer *csmode);
void cct_(f_integer *starti, f_integer *nclrs, f_integer *clrs);
void cgtxx1_(f_real *x, f_real *y, char *string, f_integer *vstat, f_integer *vconc, f_real *xconc,
             f_real *yconc, f_real *x1, f_real *y1, f_real *x2, f_real *y2, f_real *x3, f_real *y3,
             f_real *x4, f_real *y4, f_integer *string_size);
void cqprl_(f_integer *vstat, f_integer *maxpl, f_integer *maxdpl, f_integer *maxpg,
            f_integer *maxpgs, f_integer *maxpm, f_integer *maxcf, f_integer *maxchr,
            f_integer *maxcel, f_integer *celfil, f_integer *celaln, f_integer *comptx,
            f_integer *clofig);
void cqln_(f_integer *vstat, f_integer *npdefb, f_integer *nsetb, f_integer *maxbi,
           f_integer *dynmod, f_integer *nomwid, f_integer *minwid, f_integer *maxwid);
void cqlnt_(f_integer *nreq, f_integer *first, f_integer *vstat, f_integer *ntotal,
            f_integer *nlist, f_integer *lntyp);
void cqchh1_(char *font, f_integer *txp, f_integer *nreq, f_integer *first, f_integer *vstat,
             f_integer *ntotal, f_integer *nlist, f_integer *chhit, f_integer *font_size);
void cqfl_(f_integer *vstat, f_integer *npdefb, f_integer *nsetb, f_integer *maxbi,
           f_integer *dynmod, f_integer *ninsty, f_integer *instyl, f_integer *npdefp,
           f_integer *nsetp, f_integer *maxpi, f_integer *pdiv, f_integer *maxpx, f_integer *maxpy,
           f_integer *ptrans);
void cqc_(f_integer *vstat, f_integer *nsimul, f_integer *navail, f_integer *nint, f_integer *cmode,
          f_integer *dynmod, f_integer *overit, f_integer *monoc);
void cqlna_(f_integer *vstat, f_integer *lnbi, f_integer *lntyp, f_integer *lwmode,
            f_integer *lnwid, f_integer *csmode, f_integer *lnclr, f_integer *lcmode);
void cqtxa_(f_integer *vstat, f_integer *txbi, f_integer *fonti, f_integer *fontp, f_real *chexp,
            f_real *chspac, f_integer *csmode, f_integer *txclr[3], f_real *chhit,
            f_real *orient[4], f_integer *txpath, f_integer *horal, f_real *contha,
            f_integer *veral, f_real *contva, f_integer *chsi, f_integer *achsi);
void cqcte_(f_integer *nreq, f_integer *first, f_integer *vstat, f_integer *ntotal,
            f_integer *nlist, f_integer *colors);
void crqlc_(f_integer *idev, f_real *timeout, f_integer *vstat, f_integer *rstat, f_integer *mvalid,
            f_integer *triggr, f_real *x, f_real *y);
void cili_(f_integer *iclass, f_integer *idev);
void cpxa_(f_real *x, f_real *y, f_integer *nx, f_integer *ny, f_integer *xscal, f_integer *yscal,
           f_integer *xdir, f_integer *ydir, f_integer *pxclrs);
#endif
