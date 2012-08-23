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

#ifndef __Panzer_SGUtilities_hpp__
#define __Panzer_SGUtilities_hpp__

#include "Panzer_config.hpp"

#ifdef HAVE_STOKHOS

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Panzer_Traits.hpp"

#include <vector>

namespace panzer {
namespace sg_utils {

   //! Converts a vector value into a scalar depending on the scalar type
   template <typename ScalarT>
   inline void vectorToValue(const std::vector<double> & in_vector,
                             const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                             ScalarT & value)
   { value = in_vector[0]; }


   //! Converts a vector value into a scalar depending on the scalar type (specific for SGResidual)
   inline void vectorToValue(const std::vector<double> & in_vector,
                             const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                             panzer::Traits::SGResidual::ScalarT & value)
   { 
      value = 0.0; // initialize value

      // if expansion is null, do something sensible
      if(Teuchos::is_null(expansion)) {
         if(in_vector.size()>0)
            value = in_vector[0];
         return;
      }

      TEUCHOS_ASSERT(!Teuchos::is_null(expansion)); // fail if expansion is null: should never happen

      value.reset(expansion); // jamb in expansion here
      value.copyForWrite();

      TEUCHOS_ASSERT((int) in_vector.size()<=value.size()); // fail if too many coefficients

      // now blindly use fastAccessCoeff
      for(std::size_t i=0;i<in_vector.size();i++)
         value.fastAccessCoeff(i) = in_vector[i];
      for(int i=in_vector.size();i<value.size();i++) // make sure to zero out unused values
         value.fastAccessCoeff(i) = 0.0;
   }

   //! Converts a vector value into a scalar depending on the scalar type (specific for SGJacobian)
   inline void vectorToValue(const std::vector<double> & in_vector,
                             const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion,
                             panzer::Traits::SGJacobian::ScalarT & value)
   { 
      value = 0.0; // initialize value

      // if expansion is null, do something sensible
      if(Teuchos::is_null(expansion)) {
         if(in_vector.size()>0)
            value = in_vector[0];
         return;
      }

      TEUCHOS_ASSERT(!Teuchos::is_null(expansion)); // fail if expansion is null: should never happen

      value = 0.0;
      value.val().reset(expansion); // jamb in expansion here
      value.val().copyForWrite();

      TEUCHOS_ASSERT((int) in_vector.size()<=value.val().size()); // fail if too many coefficients

      // now blindly use fastAccessCoeff
      for(std::size_t i=0;i<in_vector.size();i++)
         value.val().fastAccessCoeff(i) = in_vector[i];
      for(int i=in_vector.size();i<value.val().size();i++)
         value.val().fastAccessCoeff(i) = 0.0;
   }
}
}

#endif
#endif 
