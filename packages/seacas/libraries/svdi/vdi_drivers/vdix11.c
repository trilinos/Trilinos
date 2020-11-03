/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
** /

/* The following ifdef redefines entry points for use of this X driver
   as an SVDI driver on systems which require an underscore "_" appended
   to the routine name for FORTRAN - C interfacing.
*/
#ifdef ADD_
#define vdinit vdinit_
#define vdterm vdterm_
#define vdfram vdfram_
#define vdiqdc vdiqdc_
#define vdnwpg vdnwpg_
#define vdbell vdbell_
#define vdwait vdwait_
#define vdbufl vdbufl_
#define vdstco vdstco_
#define vdiqco vdiqco_
#define vdescp vdescp_
#define vdiqes vdiqes_
#define vdiqnd vdiqnd_
#define vdmova vdmova_
#define vdlina vdlina_
#define vdpnta vdpnta_
#define vdtext vdtext_
#define vdpoly vdpoly_
#define vdiqcp vdiqcp_
#define vdiqos vdiqos_
#define vdstos vdstos_
#define vdstfc vdstfc_
#define vdstbc vdstbc_
#define vdstin vdstin_
#define vdstls vdstls_
#define vdstlw vdstlw_
#define vdstcs vdstcs_
#define vdaabu vdaabu_
#define vdaloc vdaloc_
#define vdabgl vdabgl_
#define vdakgl vdakgl_
#define vdstla vdstla_
#define vdloge vdloge_
#define vberrh vberrh_
#define vdmoni vdmoni_
#define vbpkg vbpkg_
#define vbiqpk vbiqpk_
#define vbiqdv vbiqdv_
#define vbdev vbdev_
#endif

/* The following ifdefs redefine entry points for use of this X driver
   as a CGI driver via the CGI-SVDI shell .
*/
#if defined(HP)
#define vdmova wx11mv
#define vdmoni wx11mo
#define vdgnam wx11gn
#define vbiqdv wx11iv
#define vbiqpk wx11qp
#define vdlina wx11ln
#define vdtext wx11tx
#define vdpnta wx11pt
#define vdpoly wx11py
#define vdiqcp wx11cp
#define vdstos wx11os
#define vdiqos wx11io
#define vdstfc wx11fc
#define vdstbc wx11bc
#define vdstin wx11in
#define vdstls wx11ls
#define vdstlw wx11lw
#define vdstcs wx11cs
#define vdaabu wx11bu
#define vdaloc wx11lo
#define vdabgl wx11bl
#define vdakgl wx11kl
#define vdstla wx11la
#define vdinit wx11nt
#define vdfram wx11fr
#define vdterm wx11tr
#define vdiqdc wx11dc
#define vdnwpg wx11pg
#define vdbell wx11be
#define vdwait wx11wt
#define vdbufl wx11fl
#define vdstco wx11co
#define vdiqco wx11ic
#define vdescp wx11es
#define vdiqes wx11ie
#define vdiqnd wx11id
#define vimova wx11im
#define vilina wx11il
#define vipnta wx11ip
#define vitext wx11ix
#define viinit wx11ii
#define viterm wx11it
#define vinwpg wx11ig
#define vcjob vcjob
#define vberrh wx11er
#define vdloge wx11le
#define cdrwfs wx11wf
#define cdrrfs wx11rf
#define cdrofs wx11of
#define cdrof3 wx11o3
#define cdrcfs wx11cf
#define cdroff wx11ff
#define cdroab wx11ab
#define bgpbuf wx11bf
#define qmsbuf wx11qm
#define qmsbu1 wx11bf
#define ddcbuf wx11bf
#define h75buf wx11bf
#define btkbuf wx11bf
#define nmtbuf wx11bf
#define vbimbf wx11ib
#define vbpkg wx11pk
#define vbdev wx11dv
#define vdiqrs wx11qr
#define vdstmp wx11mp
#define vdstrs wx11rs
#define vdstrv wx11rv
#define vdbrgb wx11bg
#define vdfrgb wx11fg
#define vdpixl wx11px
#define vdpixi wx11pi
#define vdrpix wx11rp
#define vdrpxi wx11ri
#define vdrscl wx11rl
#define vdiqci wx11ci
#define vbstmp wx1101
#define vifram wx1102
#define vcndcm wx1103
#define vcattr wx1104
#define vbini1 wx1105
#define vb2hls wx1106
#define vb2rgb wx1107
#define vccolt wx1108
#define vccrps wx1109
#define vcscal wx1110
#define vcddim wx1111
#define vipoly wx1112
#define vbout wx1113
#define cgixxx cgix11
#endif
#if defined(CRAY)
#define vdmova WX11MV
#define vdmoni WX11MO
#define vdgnam WX11GN
#define vbiqdv WX11IV
#define vbiqpk WX11QP
#define vdlina WX11LN
#define vdtext WX11TX
#define vdpnta WX11PT
#define vdpoly WX11PY
#define vdiqcp WX11CP
#define vdstos WX11OS
#define vdiqos WX11IO
#define vdstfc WX11FC
#define vdstbc WX11BC
#define vdstin WX11IN
#define vdstls WX11LS
#define vdstlw WX11LW
#define vdstcs WX11CS
#define vdaabu WX11BU
#define vdaloc WX11LO
#define vdabgl WX11BL
#define vdakgl WX11KL
#define vdstla WX11LA
#define vdinit WX11NT
#define vdfram WX11FR
#define vdterm WX11TR
#define vdiqdc WX11DC
#define vdnwpg WX11PG
#define vdbell WX11BE
#define vdwait WX11WT
#define vdbufl WX11FL
#define vdstco WX11CO
#define vdiqco WX11IC
#define vdescp WX11ES
#define vdiqes WX11IE
#define vdiqnd WX11ID
#define vimova WX11IM
#define vilina WX11IL
#define vipnta WX11IP
#define vitext WX11IX
#define viinit WX11II
#define viterm WX11IT
#define vinwpg WX11IG
#define cdrcom CDRCOM
#define vcjob VCJOB
#define vconod VCONOD
#define vberrh WX11ER
#define vdloge WX11LE
#define cdrwfs WX11WF
#define cdrrfs WX11RF
#define cdrofs WX11OF
#define cdrof3 WX11O3
#define cdrcfs WX11CF
#define cdroff WX11FF
#define cdroab WX11AB
#define bgpbuf WX11BF
#define qmsbuf WX11QM
#define qmsbu1 WX11BF
#define ddcbuf WX11BF
#define h75buf WX11BF
#define btkbuf WX11BF
#define nmtbuf WX11BF
#define vbimbf WX11IB
#define vbpkg WX11PK
#define vbdev WX11DV
#define vdiqrs WX11QR
#define vdstmp WX11MP
#define vdstrs WX11RS
#define vdstrv WX11RV
#define vdbrgb WX11BG
#define vdfrgb WX11FG
#define vdpixl WX11PX
#define vdpixi WX11PI
#define vdrpix WX11RP
#define vdrpxi WX11RI
#define vdrscl WX11RL
#define vdiqci WX11CI
#define vbstmp WX1101
#define vifram WX1102
#define vcndcm WX1103
#define vcattr WX1104
#define vbini1 WX1105
#define vb2hls WX1106
#define vb2rgb WX1107
#define vccolt WX1108
#define vccrps WX1109
#define vcscal WX1110
#define vcddim WX1111
#define vipoly WX1112
#define vbout WX1113
#define wx11zz WX11ZZ
#define cgixxx CGIX11
#endif
#if defined(SUN) || defined(DEC) || defined(ALLIANT)
#define vdmova wx11mv_
#define vdmoni wx11mo_
#define vdgnam wx11gn_
#define vbiqdv wx11iv_
#define vbiqpk wx11qp_
#define vdlina wx11ln_
#define vdtext wx11tx_
#define vdpnta wx11pt_
#define vdpoly wx11py_
#define vdiqcp wx11cp_
#define vdstos wx11os_
#define vdiqos wx11io_
#define vdstfc wx11fc_
#define vdstbc wx11bc_
#define vdstin wx11in_
#define vdstls wx11ls_
#define vdstlw wx11lw_
#define vdstcs wx11cs_
#define vdaabu wx11bu_
#define vdaloc wx11lo_
#define vdabgl wx11bl_
#define vdakgl wx11kl_
#define vdstla wx11la_
#define vdinit wx11nt_
#define vdfram wx11fr_
#define vdterm wx11tr_
#define vdiqdc wx11dc_
#define vdnwpg wx11pg_
#define vdbell wx11be_
#define vdwait wx11wt_
#define vdbufl wx11fl_
#define vdstco wx11co_
#define vdiqco wx11ic_
#define vdescp wx11es_
#define vdiqes wx11ie_
#define vdiqnd wx11id_
#define vimova wx11im_
#define vilina wx11il_
#define vipnta wx11ip_
#define vitext wx11ix_
#define viinit wx11ii_
#define viterm wx11it_
#define vinwpg wx11ig_
#define cdrcom cdrcom_
#define vcjob vcjob_
#define vconod vconod_
#define vberrh wx11er_
#define vdloge wx11le_
#define cdrwfs wx11wf_
#define cdrrfs wx11rf_
#define cdrofs wx11of_
#define cdrof3 wx11o3_
#define cdrcfs wx11cf_
#define cdroff wx11ff_
#define cdroab wx11ab_
#define bgpbuf wx11bf_
#define qmsbuf wx11qm_
#define qmsbu1 wx11bf_
#define ddcbuf wx11bf_
#define h75buf wx11bf_
#define btkbuf wx11bf_
#define nmtbuf wx11bf_
#define vbimbf wx11ib_
#define vbpkg wx11pk_
#define vbdev wx11dv_
#define vdiqrs wx11qr_
#define vdstmp wx11mp_
#define vdstrs wx11rs_
#define vdstrv wx11rv_
#define vdbrgb wx11bg_
#define vdfrgb wx11fg_
#define vdpixl wx11px_
#define vdpixi wx11pi_
#define vdrpix wx11rp_
#define vdrpxi wx11ri_
#define vdrscl wx11rl_
#define vdiqci wx11ci_
#define vbstmp wx1101_
#define vifram wx1102_
#define vcndcm wx1103_
#define vcattr wx1104_
#define vbini1 wx1105_
#define vb2hls wx1106_
#define vb2rgb wx1107_
#define vccolt wx1108_
#define vccrps wx1109_
#define vcscal wx1110_
#define vcddim wx1111_
#define vipoly wx1112_
#define vbout wx1113_
#define wx11bf wx11bf_
#define wx11zz wx11zz_
#define cgixxx cgix11_
#endif
/* ---------------------------- */
/*
 SVDI X Windows driver - X11 version
    Implemented in C, callable by Fortran
    (Note: Implements SVDI at the "VI" level -- some routines are "VI"
     routines, e.g. viinit, viterm, vilina, etc., while others are "VD")
    Dino Pavlakos   Jun 1989
*/

#include <stdio.h>

#include <X11/Xatom.h>
#include <X11/Xlib.h>
#include <X11/Xresource.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>

/* svdi state */
/* attributes
 *      vector[0] = foreground color
 *      vector[1] = background color
 *      vector[2] = intensity
 *      vector[3] = line style
 *      vector[4] = line width
 *      vector[5] = character box y
 *      vector[6] = character box x      */
#define MAX_VECTOR 7
    static float vector[MAX_VECTOR] = {0., 7., 1., 0., 0., 0., 0.};

/* current position */
static float xcp = 0.;
static float ycp = 0.;

/* device capabilities
 *   1. Erasability
 *   2. Scan Type
 *   3. Intensities (1-N)
 *   4. Colors (1-N)
 *   5. Line Widths (1-N)
 *   6. Line Styles (0-N)
 *   7. Character Sizes (0-N)
 *   8. Number of Locator Devices
 *   9. Number of Valuator Devices
 *  10. Number of Button Devices
 *  11. Number of Keyboard Devices
 *  12. Number of Stroke Devices
 *  13. Input
 *  14. Input Timing
 *  15. X Dimension of View Surface in Device Coordinates
 *  16. Y Dimension of View Surface in Device Coordinates
 *  17. X Dimension of View Surface in Physical Units (mm)
 *  18. Y Dimension of View Surface in Physical Units (mm)
 *  19. Smallest Line Width (DC) at default intensity
 *  20. Smallest Point (DC) at default intensity
 *  21. Smallest Character Size (DC)
 *  22. Header and Trailer Frames Required (0=no,1=yes)
 *  23. Device Identifier
 *  24. Polygon support level
 *  25. Maximum number of points in a polygon
 *  26. Setable color table
 *  27. Device color palette size (1-N)
 *  28. Direct color space size (0-N)
 *  29. Vector verses Raster VDI
 *  30. Maximum character height (DC)
 *  31. Maximum line width (DC)
 *  32. Color verses monochrome (greyscale) device
 *  33. Device pixel aspect                                      */

#define MAX_DEV_CAP 33
#define VBUF_SIZE 1024
static float dev_cap[MAX_DEV_CAP] = {4., 1., 1.,  256., 1.,
                                     3., 1., 1.,  1.,   1.,
                                     1., 1., 2.,  0.,   0.,
                                     0., 0., 0.,  1.,   1.,
                                     0., 0., 35., 3.,   (float)(VBUF_SIZE - 1),
                                     1., 0., 0.,  0.,   0.,
                                     1., 1., 1.};

/* SVDI internal parameters */
/*  Maximum NDC Values */
static float ndc_xmax = 1.;
static float ndc_ymax = 1.;

/* Logical Color Table */
#define MAX_COLORS 256
static float color_table[MAX_COLORS][3];
static float default_color_table[8][3] = {{0., 0., 0.}, {1., 0., 0.}, {0., 1., 0.}, {1., 1., 0.},
                                          {0., 0., 1.}, {1., 0., 1.}, {0., 1., 1.}, {1., 1., 1.}};
static int   x_colors[MAX_COLORS];
static enum CT { MONO, PSEUDO, FULL } color_type;
static int def_bc_index, def_fc_index;
static int svdicolors_base;

/* aspect, justification */
static float asp;
static int   just;

/* transformation parameters*/
static int   xpad, ypad;
static float scale;

/* Other globals for X */
static Window window_id;       /* window identifier */
static int    window_stat = 0; /* status indicating whether window
                                  has changed (e.g. hidden/exposed,
                                  resized, etc.) since new page or
                                  since the last time the
                                  application checked (vdescp 3501)
                               */
static Display *    display;
static int          screen;
static Screen *     screen_pnt;
static GC           gc; /* graphics context */
static GContext     gcontext;
static XPoint       vlist[VBUF_SIZE]; /* vertex buffer (polyline, polygon) */
static int          nvert = 0;
static unsigned int x_width, x_height; /* window size */
static unsigned int line_width;        /* line width in device coord. */
static Colormap     cmap;              /* color map */
static Visual *     visual;
static int          ncolors;   /* number of colors supported */
static int          f_color;   /* device foreground color */
static int          b_color;   /* device background color */
static int          line_type; /* device line type */
static XFontStruct *font_info; /* text stuff */
/* static char *font_name = "9x15"; */
static char *font_name = "fixed";
static Font  font_id;
static int   font_height, font_width; /* char size in device coord. */

/* macros which map ndc into X Windows device coords. */
#define map_x(xin) ((int)(xpad + scale * (xin)))
#define map_y(yin) ((int)(x_height - (ypad + scale * (yin))))

/* macros which map X Windows coords. into ndc */
#define ndc_map_x(xin) ((float)(((xin)-xpad) / scale))
#define ndc_map_y(yin) ((float)(((x_height - (yin)) - ypad) / scale))

/* macro to convert measure in X window units into ndc units */
#define ndc_units(in) ((float)((in) / scale))
/* macro to convert measure in ndc into X window measure */
#define x_units(in) ((int)((in)*scale))

/* macro to convert ascii(integer) to char (note: machine dependent) */
#define a_to_c(ain) ((char)(ain)) /* for ascii machine */

/* misc macros */
#define min(p1, p2) ((p1) < (p2) ? p1 : p2)
#define max(p1, p2) ((p1) > (p2) ? p1 : p2)

/* flush polyline buffer */
/* implemented as macro to save the overhead of subroutine call */
#define vflush()                                                                                   \
  {                                                                                                \
    if (nvert > 1) {                                                                               \
      XDrawLines(display, window_id, gc, vlist, nvert, CoordModeOrigin);                           \
    }                                                                                              \
    nvert = 0;                                                                                     \
  }

viinit(aspect, justif) float *aspect;
int *justif;
{
  int                  i, j, xpos, ypos;
  int                  index[8];
  unsigned long        valuemask;
  unsigned int         d_width, d_height;
  char *               geometry = NULL;
  XSetWindowAttributes setwinattr;
  XWMHints             wmhints;
  XSizeHints           hints;
  XEvent               x_event;
  XColor               def_cmap[MAX_COLORS];
  int                  nreserve, backing_store = 0;

  asp  = *aspect;
  just = *justif;

  if (asp < 0.) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 721, 5);
    asp = 0.;
  }

  if (just < 0 || just > 9) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 720, 5);
    just = 0;
  }

  /* open the display specified by the env variable DISPLAY */
  if ((display = XOpenDisplay(NULL)) == NULL) {
    fprintf(stderr, "svdi_x: cannot create window on %s\n", XDisplayName(NULL));
    exit(1);
  }
  screen     = DefaultScreen(display);
  screen_pnt = DefaultScreenOfDisplay(display);

  /* get display dimensions */
  d_width  = DisplayWidth(display, screen);
  d_height = DisplayHeight(display, screen);

  /* create a window */
  /* color attributes */
  setwinattr.background_pixel = BlackPixel(display, screen);
  setwinattr.border_pixel     = WhitePixel(display, screen);
  if (DoesBackingStore(screen_pnt)) { /* backing store */
    setwinattr.backing_store = Always;
    backing_store            = 1;
    /* set attribute mask */
    valuemask = CWBackPixel | CWBorderPixel | CWBackingStore;
  }
  else
    valuemask = CWBackPixel | CWBorderPixel;
  x_width  = .65 * d_width; /* default size */
  x_height = .65 * d_height;
  xpos     = .3 * d_width;
  ypos     = .3 * d_height;
  /* use user-supplied default geometry if available */
  geometry = XGetDefault(display, "svdi", "Geometry");
  if (geometry) {
    XParseGeometry(geometry, &xpos, &ypos, &x_width, &x_height);
  }
  window_id = XCreateWindow(display, RootWindow(display, screen), xpos, ypos, x_width, x_height, 3,
                            0, InputOutput, CopyFromParent, valuemask, &setwinattr);

  /* setup window properties */
  hints.flags      = PSize | PMinSize;
  hints.width      = x_width;
  hints.height     = x_height;
  hints.min_width  = 100;
  hints.min_height = 100;
  XSetStandardProperties(display, window_id, "svdi", "svdi", None, 0, NULL, &hints);
  wmhints.flags = InputHint;
  wmhints.input = True;
  XSetWMHints(display, window_id, &wmhints);

  /* setup a graphics context */
  gc = XCreateGC(display, window_id, 0, NULL);

  /* setup defaults */
  line_width = 0;         /* line width */
  line_type  = LineSolid; /* line type */
  XSetLineAttributes(display, gc, line_width, line_type, CapButt, JoinMiter);

  /* setup colors */
  visual  = DefaultVisual(display, screen);   /* get visual info */
  cmap    = DefaultColormap(display, screen); /* get default colormap */
  ncolors = DisplayCells(display, screen);    /* get size of color map */
  if (ncolors > MAX_COLORS)
    ncolors = MAX_COLORS;

  if (visual->class == TrueColor || visual->class == DirectColor ||
      (visual->class == PseudoColor && ncolors >= 8)) {
    /* setup color */
    if (visual->class == TrueColor || visual->class == DirectColor) {
      color_type = FULL;
      /* printf("** full color **\n"); */
    }
    else {
      color_type = PSEUDO;
      /* printf("** pseudo color, %d colors **\n", ncolors); */
      /* For pseudo-color systems, a separate virtual color map is created
         for use by SVDI.  In order to minimize the effect of SVDI color
         table changes on the rest of the screen colors, the original
         colors are copied into the new virtual map, and certain entries
         are preserved.
      */
      /* query default colors */
      for (i = 0; i < ncolors; i++)
        def_cmap[i].pixel = (unsigned long)i;
      XQueryColors(display, cmap, def_cmap, ncolors);
      /* setup a virtual colormap */
      cmap = XCreateColormap(display, window_id, visual, AllocAll);
      XSetWindowColormap(display, window_id, cmap);
      /* copy default colors into new color map */
      XStoreColors(display, cmap, def_cmap, ncolors);
      if (ncolors < 16) {
        /* number of colors is between 8 and 16 -- let driver use
           8 of them */
        svdicolors_base = ncolors - 8;
        ncolors         = 8;
      }
      else {
        /* number of colors is >= 16 -- preserve a few at each end
           of table */
        nreserve        = max(4, ncolors / 16);
        svdicolors_base = nreserve;
        ncolors         = ncolors - 2 * nreserve;
      }
    }
    dev_cap[3]   = ncolors; /* update svdi device capabilities */
    def_bc_index = 0;       /* remember defaults */
    def_fc_index = 7;
    /* setup default SVDI colors */
    for (i = 0; i < 8; i++)
      index[i] = i;
    i = 8;
    j = 0;
    vdstco(&i, index, default_color_table, &j); /* color table */
    i = 0;
    vdstbc(&i); /* background color */
    i = 7;
    vdstfc(&i); /* foreground color */
  }
  else {
    /* default to monochrome */
    color_type = MONO; /* remember color type */
    /* printf("** mono color **\n"); */
    ncolors      = 2;
    dev_cap[3]   = 2; /* update svdi device capabilities */
    def_bc_index = 0; /* remember defaults */
    def_fc_index = 1;
    /* setup default colors */
    color_table[0][0] = 0; /* background */
    color_table[0][1] = 0;
    color_table[0][2] = 0;
    b_color           = BlackPixel(display, screen);
    x_colors[0]       = b_color;
    color_table[1][0] = 1; /* foreground */
    color_table[1][1] = 1;
    color_table[1][2] = 1;
    f_color           = WhitePixel(display, screen);
    x_colors[1]       = f_color;
  }

  /* text stuff */
  if ((font_info = XLoadQueryFont(display, font_name)) == NULL) {
    fprintf(stderr, "svdi_x: font request failed, using default\n");
    gcontext  = XGContextFromGC(gc);
    font_info = XQueryFont(display, gcontext);
  }
  else {
    font_id = (*font_info).fid;
    XSetFont(display, gc, font_id);
  }
  font_height = font_info->max_bounds.ascent + font_info->max_bounds.descent;
  font_width  = font_info->max_bounds.rbearing - font_info->min_bounds.lbearing;
  /* this gives bounding box width, which seems like what
     we want ... but for whatever reason, this computation
     doesn't seem to give correct results on VMS */

  /* setup input event stuff */
  /* specify what kinds of input events to accept */
  XSelectInput(display, window_id,
               ExposureMask | StructureNotifyMask | KeyPressMask | ButtonPressMask);

  /* display the window */
  XMapWindow(display, window_id);

  /* if no backing store, wait for first expose event before drawing */
  if (!backing_store)
    XWindowEvent(display, window_id, ExposureMask, &x_event);

  /* setup parameters which are affected by dynamics */
  x_dynamics();

  /* get synced and flush event queue -- start with a clean slate */
  XSync(display, True);
}

x_dynamics()
{
  XWindowAttributes win_info;
  float             asp1;
  int               just1, xused, yused;

  /* asp and just remember the requested aspect and justification
     at initialization */
  asp1  = asp;
  just1 = just;

  /* query window information needed */
  XGetWindowAttributes(display, window_id, &win_info);
  x_width  = win_info.width;
  x_height = win_info.height;

  /* update SVDI info */
  dev_cap[14] = (float)x_width;
  dev_cap[15] = (float)x_height;

  /* handle aspect */
  if (asp1 == 0.) {
    asp1 = (float)x_width / (float)x_height;
    /* note: What to do if the window changes??? If the application
       requested the entire view surface originally (i.e. aspect = 0.),
       then, when the window changes, do we recompute aspect, which
       in turn may result in new ndc maximums, which in turn may result
       in the application drawing to different limits than currently
       exist (which means the application may end up trying to draw
       outside the window)?  Or, do we just remember the original
       resulting aspect */
    /* In the default case, just use the original aspect so the
       application is guaranteed to see everything which would draw
       based on that aspect, regardless of window */
    asp = asp1;
  }
  if (asp1 > 1.) {
    ndc_xmax = 1.;
    ndc_ymax = 1. / asp1;
  }
  else {
    ndc_xmax = asp1;
    ndc_ymax = 1.;
  }
  scale = min((x_width - 1) / ndc_xmax, (x_height - 1) / ndc_ymax);
  xused = (int)(scale * ndc_xmax);
  yused = (int)(scale * ndc_ymax);

  /* handle justification */
  if (just1 == 0)
    just1 = 5;
  /* y offset */
  if (just1 < 4)
    ypad = 0;
  else if (just1 < 7)
    ypad = (x_height - yused) / 2;
  else
    ypad = x_height - yused;
  /* x offset */
  if (just1 == 1 || just1 == 4 || just1 == 7)
    xpad = 0;
  else if (just1 == 2 || just1 == 5 || just1 == 8)
    xpad = (x_width - xused) / 2;
  else
    xpad = x_width - xused;

  /* spatial attributes */
  vector[4] = (float)(ndc_units((line_width == 0) ? 1 : line_width) * 100);
  vector[5] = (float)(ndc_units(font_height));
  vector[6] = (float)(ndc_units(font_width));
}

vifram(type) int *type;
{
}

viterm()
{
  vflush(); /* flush polyline buffer */

  /* close window, disconnect display */
  XDestroyWindow(display, window_id);
  XCloseDisplay(display);
}

vdiqdc(index, value) int *index;
float *value;
{
  if (*index < 1 || *index > MAX_DEV_CAP) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 726, 5);
    return (-1);
  }
  else
    *value = dev_cap[*index - 1];
}

vinwpg()
{
  vflush();         /* flush polyline buffer */
  x_check_window(); /* check window change */
  window_stat = 0;  /* reset window change indicator */
  XClearWindow(display, window_id);
}

vdbell()
{
  vflush(); /* flush polyline buffer */
  XBell(display, 100);
}

vdwait()
{
  XEvent x_event;
  vdbufl();

  XBell(display, 100); /* ring bell */
  /* wait for keyboard key or mouse button */
  XWindowEvent(display, window_id, KeyPressMask | ButtonPressMask, &x_event);

  x_check_window(); /* check window change */
}

vdbufl()
{
  vflush(); /* flush polyline buffer */
  /* flush output buffer, wait for server to do its thing ...
     we use Sync instead of Flush to wait for server to finish
     processing (versus just sending buffer to server) */
  XSync(display, False);
}

vdstco(num, index_array, color_array, color_mod) int *num, *color_mod;
int   index_array[];
float color_array[][3];
{
  int    i;
  XColor tcolor;

  if (color_type == MONO)
    return;

  /* check valid number of colors */
  if (*num < 1) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 723, 5);
    return (-1);
  }

  /* check valid color mod */
  if (*color_mod != 0 && *color_mod != 1) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 725, 5);
    return (-1);
  }

  /* loop over number of colors passed in */
  for (i = 0; i < *num; i++) {
    /* check index value */
    if (index_array[i] < 0 || index_array[i] > MAX_COLORS - 1) {
      fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 724, 5);
    }

    /* check index against color table size */
    else if (index_array[i] > ncolors - 1) {
      /* ignore */
    }

    /* if rgb mode */
    else if (*color_mod == 0) {
      /* check rgb values */
      if (color_array[i][0] < 0. || color_array[i][0] > 1. || color_array[i][1] < 0. ||
          color_array[i][1] > 1. || color_array[i][2] < 0. || color_array[i][2] > 1.)
        fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 727, 5);
      else /* set colors */
      {
        color_table[index_array[i]][0] = color_array[i][0];
        color_table[index_array[i]][1] = color_array[i][1];
        color_table[index_array[i]][2] = color_array[i][2];
        tcolor.red                     = (unsigned int)(65535 * color_array[i][0]);
        tcolor.green                   = (unsigned int)(65535 * color_array[i][1]);
        tcolor.blue                    = (unsigned int)(65535 * color_array[i][2]);
        if (color_type == PSEUDO) {
          tcolor.pixel = (unsigned int)(svdicolors_base + index_array[i]);
          tcolor.flags = DoRed | DoGreen | DoBlue;
          XStoreColor(display, cmap, &tcolor);
          x_colors[index_array[i]] = tcolor.pixel;
        }
        else { /* FULL */
          if (XAllocColor(display, cmap, &tcolor) == 0)
            fprintf(stderr, "svdi_x: problem doing XAllocColor\n");
          x_colors[index_array[i]] = tcolor.pixel;
        }
      }
    }

    /* hls mode not supported */
    else {
      fprintf(stderr, " HLS option being phased out - not doing anything\n");
      fprintf(stderr, " Contact Computer Graphics Group - Div. 2644\n");
    }
  }
  /* end for loop */
}

vdiqco(num, index_array, color_array, color_mod) int *num, *color_mod;
int   index_array[];
float color_array[][3];
{
  int i;

  /* check valid number of colors */
  if (*num < 1) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 723, 5);
    return (-1);
  }

  /* check valid color mod */
  if (*color_mod != 0 && *color_mod != 1) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 725, 5);
    return (-1);
  }

  /* loop over number of colors passed in */
  for (i = 0; i < *num; i++) {
    /* check valid index value */
    if (index_array[i] < 0 || index_array[i] > ncolors - 1) {
      fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 724, 5);
    }

    /* check index against color table size */
    else if (index_array[i] > ncolors - 1) {
      /* unsupported index value */
      color_array[i][0] = -1.;
    }

    /* if rgb mode */
    else if (*color_mod == 0) {
      color_array[i][0] = color_table[index_array[i]][0];
      color_array[i][1] = color_table[index_array[i]][1];
      color_array[i][2] = color_table[index_array[i]][2];
    }

    /* hls mode */
    else {
      fprintf(stderr, " HLS option being phased out - returning RGB values\n");
      fprintf(stderr, " Contact Computer Graphics Group - Div. 2644\n");
      color_array[i][0] = color_table[index_array[i]][0];
      color_array[i][1] = color_table[index_array[i]][1];
      color_array[i][2] = color_table[index_array[i]][2];
    }
  }
  /* end for loop */
}

vdescp(escape_code, n, args) int *escape_code, *n;
float args[];
{
  vflush();

  if (*n < 0) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 802, 5);
    return (-1);
  }

  switch (*escape_code) {
  case 3500: /* sync, flush event queue */ XSync(display, True); break;
  case 3501: /* check window change (expose, resize, etc.) --
                return change indicator */
    x_check_window();
    args[0]     = (float)window_stat;
    window_stat = 0; /* reset window change indicator */
    break;
  default: break;
  }
}

vdiqes(escape_code, support) int *escape_code, *support;
{
  switch (*escape_code) {
  case 3500: *support = 1; break;
  case 3501: *support = 1; break;
  default: *support = 0; break;
  }
}

vdiqnd(x_ndc, y_ndc) float *x_ndc, *y_ndc;
{
  *x_ndc = ndc_xmax;
  *y_ndc = ndc_ymax;
}

vimova(x, y) float *x, *y;
{
  vflush(); /* flush polyline buffer */

  /* update current position */
  xcp = *x;
  ycp = *y;
}

vilina(x, y) float *x, *y;
{
  /* if polyline buffer full, flush it */
  if (nvert >= VBUF_SIZE)
    vflush();

  /* if vertex buffer empty, then start a polyline at the previous
     current position */
  if (nvert == 0) {
    vlist[0].x = (short)map_x(xcp);
    vlist[0].y = (short)map_y(ycp);
    nvert      = 1;
  }

  /* add point to polyline buffer */
  vlist[nvert].x = (short)map_x(*x);
  vlist[nvert].y = (short)map_y(*y);
  /* add point only if it's not the same as the previous point */
  if ((vlist[nvert].x != vlist[nvert - 1].x) || (vlist[nvert].y != vlist[nvert - 1].y))
    nvert++;

  /* update current position */
  xcp = *x;
  ycp = *y;
}

vipnta(x, y) float *x, *y;
{
  int xx, yy;

  vflush(); /* flush polyline buffer */

  /* draw point */
  xx = map_x(*x);
  yy = map_y(*y);
  /* XDrawPoint doesn't seem to work on the HP --
     nor does XDrawLine of just a single point */
  XDrawLine(display, window_id, gc, xx, yy, xx + 1, yy);

  /* update current position */
  xcp = *x;
  ycp = *y;
}

vitext(length, char_array) int *length, char_array[];
{
  int  len, lenout, i;
  char strout[137];
  len = *length;

  vflush(); /* flush polyline buffer */

  if (len < 1) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 212, 5);
    return (-1);
  }

  if (len > 136) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 213, 5);
    len = 136;
  }

  lenout = 0; /* count characters in string output buffer "strout" */

  for (i = 0; i < len; i++) /* for each character */
  {
    if (char_array[i] < 32 || char_array[i] > 126) /* if special char */
    {
      if (lenout != 0) /* flush string output buffer */
      {
        strout[lenout] = '\0';
        XDrawString(display, window_id, gc, map_x(xcp), map_y(ycp), strout, lenout);
        xcp    = xcp + ndc_units(XTextWidth(font_info, strout, lenout));
        lenout = 0;
      }
      switch (char_array[i]) /* process special character */
      {
      case 8: /* backspace */ xcp = xcp - vector[6]; break;
      case 10: /* line feed */ ycp = ycp - vector[5]; break;
      case 13: /* carriage return */ xcp = 0.; break;
      default: /* other */
        fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 208, 5);
        break;
      }
    }
    else /* add char to string output buffer */
    {
      strout[lenout++] = a_to_c(char_array[i]);
    }
  }
  /* end for */

  /* All done, flush string output buffer */
  if (lenout != 0) /* flush string output buffer */
  {
    strout[lenout] = '\0';
    XDrawString(display, window_id, gc, map_x(xcp), map_y(ycp), strout, lenout);
    xcp = xcp + ndc_units(XTextWidth(font_info, strout, lenout));
  }

  /* make sure we end up where we want to be */
  vimova(&xcp, &ycp);
}

vipoly(x_array, y_array, npts) float x_array[], y_array[];
int *npts;
{
  int np, i;
  np = *npts;

  vflush(); /* flush polyline buffer */

  /* check size of polygon ??? */
  if (np > VBUF_SIZE)
    np = VBUF_SIZE;

  /* set up vertex list in X format */
  for (i = 0; i < np; i++) {
    vlist[i].x = (short)map_x(x_array[i]);
    vlist[i].y = (short)map_y(y_array[i]);
  }

  /* draw polygon */
  XFillPolygon(display, window_id, gc, vlist, np, Complex, CoordModeOrigin);

  /* update current position */
  /* is this where it really is ??? */
  xcp = x_array[0];
  ycp = y_array[0];
}

vdiqcp(x, y) float *x, *y;
{
  *x = xcp;
  *y = ycp;
}

vdiqos(attr_array) float attr_array[];
{
  int i;
  for (i = 0; i < MAX_VECTOR; i++)
    attr_array[i] = vector[i];
}

vdstfc(color_index) int *color_index;
{
  int index;
  index = *color_index;

  vflush(); /* flush polyline buffer */

  /* check valid color index */
  if (index < 0 || index > 255) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 724, 5);
    return (-1);
  }

  if (index > ncolors - 1)
    index = def_fc_index;

  f_color = x_colors[index];
  XSetForeground(display, gc, f_color);
  vector[0] = index;
}

vdstbc(color_index) int *color_index;
{
  int index;
  index = *color_index;

  if (index < 0 || index > 255) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 724, 5);
    return (-1);
  }

  if (index > ncolors - 1)
    index = def_bc_index;

  b_color = x_colors[index];
  XSetWindowBackground(display, window_id, b_color);
  vector[1] = index;
}

vdstin(intensity) float *intensity;
{
  if (*intensity < 0. || *intensity > 1.) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 401, 5);
    return (-1);
  }
}

vdstls(line_style) int *line_style;
{
  char dlist[4]; /* dash pattern list */
  int  nd;

  vflush(); /* flush polyline buffer */

  if (*line_style < 0 || *line_style > 5) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 401, 5);
    return (-1);
  }

  /* set line type accordingly */
  switch (*line_style) {
  case 0: line_type = LineSolid; break;
  case 1: /* dotted */
    line_type = LineOnOffDash;
    dlist[0]  = 3;
    dlist[1]  = 1;
    nd        = 2;
    break;
  case 5:
    line_type = LineOnOffDash;
    dlist[0]  = 4;
    dlist[1]  = 4;
    nd        = 2;
    break;
  default:
    line_type = LineOnOffDash;
    dlist[0]  = 4;
    dlist[1]  = 4;
    nd        = 2;
  }
  XSetLineAttributes(display, gc, line_width, line_type, CapButt, JoinMiter);
  if (line_type != LineSolid)
    XSetDashes(display, gc, 0, dlist, nd);

  vector[3] = *line_style; /* update state list */
}

vdstlw(line_wid) float *line_wid;
{
  float slw;

  vflush(); /* flush polyline buffer */

  if (*line_wid < 0. || *line_wid > 1.) {
    fprintf(stderr, " SVDI Error Number %d, Severity Code %d\n", 401, 5);
    return (-1);
  }

  slw        = (*line_wid) / 100; /* svdi wants line widths scaled by 100 */
  line_width = x_units(slw);
  if (line_width < 1)
    line_width = 1;
  XSetLineAttributes(display, gc, line_width, line_type, CapButt, JoinMiter);
  vector[4] = *line_wid;
}

vdstcs(y_size) float *y_size;
{
  /* doesn't do anything */
  /* character size is fixed for each font */
}

vdaabu(button) int *button;
{
  vdbufl();
  /* not implemented */
}

vdaloc(x, y) float *x, *y;
{
  XEvent        x_event;
  XButtonEvent *x_butp;

  x_butp = (XButtonEvent *)&x_event;
  vdbufl();

  /* wait for button pressed, pull coordinates out of event report */
  XWindowEvent(display, window_id, ButtonPressMask, x_butp);
  *x = ndc_map_x(x_butp->x);
  *y = ndc_map_y(x_butp->y);
}

vdabgl(button, x, y) int *button;
float *x, *y;
{
  vdbufl();
  /* not implemented */
}

vdakgl(charac, x, y) int *charac;
float *x, *y;
{
  XEvent     x_event;
  XKeyEvent *x_keyp;
  char       key;

  x_keyp = (XKeyEvent *)&x_event;
  vdbufl();

  /* wait for key pressed, pull coordinates, key out of event report */
  XWindowEvent(display, window_id, KeyPressMask, x_keyp);
  *x = ndc_map_x(x_keyp->x);
  *y = ndc_map_y(x_keyp->y);
  XLookupString(x_keyp, &key, 1, NULL, NULL); /* lookup which key */
  *charac = (int)key;                         /* convert to integer ascii */
  /* note: XLookupString presumably returns key in ASCII, so we
     just need to repackage as integer */
}

vdstla(x, y) float *x, *y;
{
  XWarpPointer(display, None, window_id, 0, 0, 0, 0, map_x(*x), map_y(*y));
}

vdstos(attr_array) float attr_array[];
{
  /* attr_array[0] is equivalent to *attr_array,
     attr_array[1] is equivalent to *(attr_array+1),
     etc. */
  int i;

  i = (int)attr_array[0];
  vdstfc(&i);
  i = (int)attr_array[1];
  vdstbc(&i);
  vdstin(attr_array + 2);
  i = (int)attr_array[3];
  vdstls(&i);
  vdstlw(attr_array + 4);
  vdstcs(attr_array + 5);
}

x_check_window()
{
  int    change;
  XEvent x_event;

  change = 0;
  /* pick up any window change events */
  while (XCheckWindowEvent(display, window_id, ExposureMask | StructureNotifyMask, &x_event)) {
    change = 1;
  };

  /* if we found any of the ones we care about, then recompute
     dynamic parameters, record change */
  if (change == 1) {
    x_dynamics();
    window_stat = 1; /* set window change indicator */
  }
}

vdpnta(x, y) float *x, *y;
{
  vipnta(x, y);
}
vdnwpg() { vinwpg(); }
vdmova(x, y) float *x, *y;
{
  vimova(x, y);
}
vdmoni(state) int *state;
{
}
vdloge(errnum, errsev) int *errnum, *errsev;
{
}
vberrh(errnum, errsev) int *errnum, *errsev;
{
}
vdlina(x, y) float *x, *y;
{
  vilina(x, y);
}
vdpoly(xarray, yarray, npts) float xarray[], yarray[];
int *npts;
{
  vipoly(xarray, yarray, npts);
}
vdterm() { viterm(); }
vdinit(aspect, justif) float *aspect;
int *justif;
{
  viinit(aspect, justif);
}
vdtext(length, char_array) int *length, char_array[];
{
  vitext(length, char_array);
}
vdfram(type) int *type;
{
  vifram(type);
}
vbpkg(pkg) char pkg[];
{
}
vbdev(dev) char dev[];
{
}
vbiqpk(pkg) char pkg[];
{
}
vbiqdv(dev) char dev[];
{
}

/* null entry points for raster */

vdiqrs() {}

vdstmp() {}

vdstrs() {}

vdstrv() {}

vdbrgb() {}

vdfrgb() {}

vdpixl() {}

vdpixi() {}

vdrpix() {}

vdrpxi() {}

vdrscl() {}

vdiqci() {}

vbstmp() {}
