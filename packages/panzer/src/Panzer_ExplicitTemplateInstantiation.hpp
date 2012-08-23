// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP
#define PANZER_EXPLICIT_TEMPLATE_INSTANTIATION_HPP

#include "Panzer_Traits.hpp"

// ONE template argument 
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  template class name<panzer::Traits::Residual>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  template class name<panzer::Traits::Jacobian>; 

// stochastic galerkin objects
#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_ONE_T(name) \
     template class name<panzer::Traits::SGResidual>; 

   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_ONE_T(name) \
     template class name<panzer::Traits::SGJacobian>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_ONE_T(name)
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_ONE_T(name)
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_ONE_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_ONE_T(name)

// TWO template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  template class name<panzer::Traits::Residual, panzer::Traits>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  template class name<panzer::Traits::Jacobian, panzer::Traits>; 

// stochastic galerkin objects
#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_TWO_T(name) \
     template class name<panzer::Traits::SGResidual, panzer::Traits>; 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_TWO_T(name) \
     template class name<panzer::Traits::SGJacobian, panzer::Traits>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_TWO_T(name) 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_TWO_T(name) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_TWO_T(name) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_TWO_T(name)

// THREE (one user defined) template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Residual, panzer::Traits,ExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_T(name,ExtraT) \
  template class name<panzer::Traits::Jacobian, panzer::Traits,ExtraT>; 

// stochastic galerkin objects
#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_THREE_T(name,ExtraT) \
     template class name<panzer::Traits::SGResidual, panzer::Traits,ExtraT>; 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_THREE_T(name,ExtraT) \
     template class name<panzer::Traits::SGJacobian, panzer::Traits,ExtraT>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_THREE_T(name,ExtraT) 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_THREE_T(name,ExtraT) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_THREE_T(name,ExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_THREE_T(name,ExtraT)

// THREE (two user defined) template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Residual,FirstExtraT,SecondExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Jacobian,FirstExtraT,SecondExtraT>; 

#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
     template class name<panzer::Traits::SGResidual,FirstExtraT,SecondExtraT>; 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
     template class name<panzer::Traits::SGJacobian,FirstExtraT,SecondExtraT>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_THREE_2U_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_THREE_2U_T(name,FirstExtraT,SecondExtraT)

// FOUR (two user defined) template arguments
#define PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Residual, panzer::Traits,FirstExtraT,SecondExtraT>; 

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
  template class name<panzer::Traits::Jacobian, panzer::Traits,FirstExtraT,SecondExtraT>; 

#ifdef HAVE_STOKHOS
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
     template class name<panzer::Traits::SGResidual, panzer::Traits,FirstExtraT,SecondExtraT>; 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
     template class name<panzer::Traits::SGJacobian, panzer::Traits,FirstExtraT,SecondExtraT>; 
#else
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) 
   #define PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) 
#endif

#define PANZER_INSTANTIATE_TEMPLATE_CLASS_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_RESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_JACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGRESIDUAL_FOUR_T(name,FirstExtraT,SecondExtraT) \
  PANZER_INSTANTIATE_TEMPLATE_CLASS_SGJACOBIAN_FOUR_T(name,FirstExtraT,SecondExtraT)

#endif
