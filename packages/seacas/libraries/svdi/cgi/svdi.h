/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef SVDI_H
#define SVDI_H
void vbstmp(int *);
void vdakgl(int *charac, float *x, float *y);
void vdbell(void);
void vdbrgb(float *, float *, float *);
void vdbufl(void);
void vdescp(int *escape_code, int *n, float args[]);
void vdfrgb(float *, float *, float *);
void vdinit(float *aspect, int *justif);
void vdiqci(float *, float *, float *, int *);
void vdiqco(int *num, int *index_array, float color_array[][3], int *color_mod);
void vdiqdc(int *index, float *value);
void vdiqes(int *escape_code, int *support);
void vdiqnd(float *x_ndc, float *y_ndc);
void vdiqos(float *attr_array);
void vdiqrs(int *, float *);
void vdlina(float *x, float *y);
void vdmova(float *x, float *y);
void vdnwpg(void);
void vdpixi(int *, int *, int *, int *);
void vdpixl(int *, int *, float *, float *, float *, int *);
void vdpnta(float *x, float *y);
void vdpoly(float *xarray, float *yarray, int *npts);
void vdstbc(int *color_index);
void vdstco(int *num, int *index_array, float color_array[][3], int *color_mod);
void vdstcs(float *y_size);
void vdstfc(int *color_index);
void vdstls(int *line_style);
void vdstlw(float *line_wid);
void vdstos(float *attr_array);
void vdstrs(int *nx1, int *nx2);
void vdstrv(float *xmin, float *xmax, float *ymin, float *ymax);
void vdterm(void);
void vdtext(int *length, int *char_array);
#endif
