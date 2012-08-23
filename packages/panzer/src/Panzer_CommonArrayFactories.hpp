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

#ifndef PANZER_COMMON_ARRAY_FACTORIES_HPP
#define PANZER_COMMON_ARRAY_FACTORIES_HPP

#include "Intrepid_FieldContainer.hpp"
#include "Phalanx_MDField.hpp"

#include <string>

//
// This file contains several common array factories
// useful for build arrays through the BasisValues and
// IntegrationValues classes. In particular these objects
// are used in the <code>setupArrays</code> functions.
// Because these class are used as a template argument the
// types and names used are very specific to the BasisValues
// interface.
//

namespace panzer {
  
  /** Implementation for intrepid field container factory. This
    * is intended to be used only with the BasisValues and
    * IntegrationValues objects. Notice in this case the string
    * argument is not used.
    */
  class IntrepidFieldContainerFactory {
  public:
     template <typename Scalar,typename T0>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     Intrepid::FieldContainer<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;
  };

  /** Implementation for MDField array factory. This
    * is intended to be used only with the BasisValues and
    * IntegrationValues objects.
    */
  class MDFieldArrayFactory {
  public:
     /** Build fields with no prefix, will simply use the string
       * passed into <code>buildArray</code> to name the fields.
       */
     MDFieldArrayFactory() : prefix_("") {}

     /** Build fields with a prefix, will use the string
       * passed into <code>buildArray</code> prefixed with the
       * argument to this constructor to name the fields.
       */
     MDFieldArrayFactory(const std::string & prefix) : prefix_(prefix) {}

 
     template <typename Scalar,typename T0>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0) const;
     template <typename Scalar,typename T0,typename T1>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1) const;
     template <typename Scalar,typename T0,typename T1,typename T2>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3) const;
     template <typename Scalar,typename T0,typename T1,typename T2,typename T3,typename T4>
     PHX::MDField<Scalar> buildArray(const std::string & str,int d0,int d1,int d2,int d3,int d4) const;

  private:
     std::string prefix_;     
  };

} // namespace panzer

#include "Panzer_CommonArrayFactories_impl.hpp"

#endif
