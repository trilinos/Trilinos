/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
