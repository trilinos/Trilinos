// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MUELU_MP_VECTOR_HPP
#define STOKHOS_MUELU_MP_VECTOR_HPP

// This header file should be included whenever compiling any MueLu
// code with Stokhos scalar types

// MP includes and specializations
#include "Stokhos_Tpetra_MP_Vector.hpp"
#include "MueLu_config.hpp"
#if defined(HAVE_MUELU_IFPACK2)
#include "Stokhos_Ifpack2_MP_Vector.hpp"
#endif
#if defined(HAVE_MUELU_AMESOS2)
#include "Stokhos_Amesos2_MP_Vector.hpp"
#endif
#if defined(HAVE_MUELU_BELOS)
#include "BelosXpetraAdapterMultiVector_MP_Vector.hpp"
#endif

#endif // STOKHOS_MUELU_MP_VECTOR_HPP
