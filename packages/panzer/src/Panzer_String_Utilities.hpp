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

#ifndef PANZER_STRING_UTILITIES_HPP
#define PANZER_STRING_UTILITIES_HPP

#include <vector>
#include <string>

#include "Panzer_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Panzer_SGUtilities.hpp"

#ifdef HAVE_STOKHOS
#include "Stokhos_OrthogPolyExpansion.hpp"
#endif

namespace panzer {

  //! Tokenize a string, put tokens in a vector
  void StringTokenizer(std::vector<std::string>& tokens,
		       const std::string& str,
		       const std::string delimiter = ",",bool trim=false);

  //! Turn a vector of tokens into a vector of doubles
  void TokensToDoubles(std::vector<double> & values,const std::vector<std::string> & tokens);

  //! Turn a vector of tokens into a vector of ints
  void TokensToInts(std::vector<int> & values,const std::vector<std::string> & tokens);

  /** Read in a parameter field and return the correct scalar field. This parses
    * both scalar type data and stochastic galerkin type data
    */
  template <typename ScalarT>
  #ifdef HAVE_STOKHOS
  ScalarT getScalarParameter(const std::string & field,const Teuchos::ParameterList & plist,const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > & expansion)
  #else 
  ScalarT getScalarParameter(const std::string & field,const Teuchos::ParameterList & plist)
  #endif
  {
    #ifdef HAVE_STOKHOS
      if(plist.isType<double>(field)) {
        return plist.get<double>(field);
      }
      else if(plist.isType<std::string>(field)) {
        ScalarT value;


        // this is essentially only for stochastic galerkin
        std::vector<std::string> coeffs_tok;
        std::vector<double> coeffs;
        panzer::StringTokenizer(coeffs_tok,plist.get<std::string>(field),",",true);
        panzer::TokensToDoubles(coeffs,coeffs_tok);

        if(coeffs.size()==1)
           return coeffs[0];

        // we only need expansion if more than one term is needed
        TEUCHOS_ASSERT(expansion!=Teuchos::null);
        panzer::sg_utils::vectorToValue(coeffs,expansion,value);

        return value;
      }
      else {
        TEUCHOS_ASSERT(false);
      }
    #else
      return plist.get<double>(field);
    #endif
  }
}

#endif
