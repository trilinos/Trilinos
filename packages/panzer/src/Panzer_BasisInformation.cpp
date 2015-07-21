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

#include "Panzer_BasisInformation.hpp"

#include "Teuchos_Assert.hpp"

#include <sstream>

namespace panzer {

BasisInformation::
BasisInformation(const std::string & in_basis_type,
	         const int in_basis_order,
	         const shards::CellTopology & cell_topo) :
  topology_(cell_topo)
{
  // Support for deprecated basis descriptions
  std::string basis_type = in_basis_type;
  int basis_order = in_basis_order;

  if (basis_type=="Q1" || basis_type=="T1") {
    basis_type = "HGrad";
    basis_order = 1;
  }
  else if (basis_type == "Q2" || basis_type=="T2") {
    basis_type = "HGrad";
    basis_order = 2;
  }
  else if (basis_type == "TEdge1" || basis_type=="QEdge1") {
    basis_type = "HCurl";
    basis_order = 1;
  }
  else if(basis_type == "Const") {
    basis_type = "Const";
    basis_order = 0;
  }
  // End deprecated basis support

  basis_type_ = basis_type;
  basis_order_ = basis_order;

  if(  basis_type_ == "HGrad")
    element_space_ = HGRAD;
  else if(basis_type_=="HCurl")
    element_space_ = HCURL;
  else if(basis_type_=="HDiv")
    element_space_ = HDIV;
  else if(basis_type_=="Const")
    element_space_ = CONST;
  else { TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
				    "BasisInformation::initializeIntrospection - Invalid basis name \"" 
				    << basis_type_ << "\""); }
}

}
