// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_MP_VECTOR_HPP
#define SACADO_FAD_MP_VECTOR_HPP

#include "Sacado_ConfigDefs.h"

#ifdef SACADO_NEW_FAD_DESIGN_IS_DEFAULT
#include "Sacado_Fad_Exp_MP_Vector.hpp"
#else
#include "Sacado_Fad_DFad_MP_Vector.hpp"
#include "Sacado_Fad_SFad_MP_Vector.hpp"
#include "Sacado_Fad_SLFad_MP_Vector.hpp"
#include "Sacado_Fad_ViewFad_MP_Vector.hpp"
#endif

#endif // SACADO_FAD_MP_VECTOR_HPP
