// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef INTREPID2_ETI_MACROS_HPP
#define INTREPID2_ETI_MACROS_HPP

/*
The two macros below allow to ETI with respect to one or more Kokkos devices.
Assume that we define a macro like so:

#define HGRAD_WEDGE_DEG2_FEM_INSTANT(DEVICE, OUTPUT_TYPE, POINT_TYPE, EXTERN)                           \
  EXTERN template class Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM<true, DEVICE,    \
                                                              double, double>; \
  EXTERN template class Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM<false, DEVICE,   \
                                                              double, double>;

Then we can call

INTREPID2_ETI_DEVICE_DEF(HGRAD_WEDGE_DEG2_FEM_INSTANT);

in a header and

INTREPID2_ETI_DEVICE(HGRAD_WEDGE_DEG2_FEM_INSTANT);

in a cpp to perform ETI. Un-ETI'd template parameter combinations will not cause
undefined references, but might get compiled more than once.

 */

#define INTREPID2_ETI_DEVICE_DEF(INSTMACRO)                                    \
  INSTMACRO(Kokkos::DefaultExecutionSpace::device_type, double, double,        \
            extern)

#define INTREPID2_ETI_DEVICE(INSTMACRO)                                        \
  INSTMACRO(Kokkos::DefaultExecutionSpace::device_type, double, double, )



#endif
