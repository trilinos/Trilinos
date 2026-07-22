// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef INLINEGEOMETRIES_H
#define INLINEGEOMETRIES_H  
  enum InlineGeometryType {
      INLINE_CARTESIAN,
      INLINE_SPHERICAL,
      INLINE_CYLINDRICAL,
      RADIAL,
      RADIAL_TRISECTION,
      UNKNOWN
    };

  enum InlineDecompositionType {
      BISECTION,
      PROCESSOR_LAYOUT,
      SEQUENTIAL,
      RANDOM,
      UNKNOWN_DECOMPOSITION
    };
enum inlineParameters{
  INLINE_PARAM_NX,
  INLINE_PARAM_NY,
  INLINE_PARAM_NZ,
  INLINE_PARAM_BX,
  INLINE_PARAM_BY,
  INLINE_PARAM_BZ,
  INLINE_PARAM_RI,
  INLINE_PARAM_RO,
  INLINE_PARAM_THETA,
  INLINE_PARAM_PHI,
  INLINE_PARAM_XMIN,
  INLINE_PARAM_YMIN,
  INLINE_PARAM_ZMIN,
  INLINE_PARAM_ZMAX,
  INLINE_PARAM_NTHETA,
  INLINE_PARAM_NPHI,
  INLINE_PARAM_NR,
  INLINE_PARAM_BR,
  INLINE_PARAM_BTHETA,
  INLINE_PARAM_BPHI,
  INLINE_PARAM_GMINX,
  INLINE_PARAM_GMINY,
  INLINE_PARAM_GMINZ,
  INLINE_PARAM_GMAXX,
  INLINE_PARAM_GMAXY,
  INLINE_PARAM_GMAXZ
};
enum BlockType    {RECTILINEAR, CURVILINEAR};

#endif
