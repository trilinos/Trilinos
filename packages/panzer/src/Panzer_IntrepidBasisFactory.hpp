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

#ifndef PANZER_INTREPID_BASIS_FACTORY_H
#define PANZER_INTREPID_BASIS_FACTORY_H

#include <sstream>
#include <string>
#include <map>
#include "Teuchos_RCP.hpp"
#include "Intrepid_Basis.hpp"

#include "Shards_CellTopology.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_LINE_C1_FEM.hpp"
#include <Intrepid_HCURL_TRI_I1_FEM.hpp>
#include <Intrepid_HCURL_QUAD_I1_FEM.hpp>
#include <Intrepid_HCURL_HEX_I1_FEM.hpp>

namespace panzer {

  template <typename ScalarT, typename ArrayT>
    Teuchos::RCP<Intrepid::Basis<ScalarT,ArrayT> >
    createIntrepidBasis(const std::string type, int cell_dimension,const Teuchos::RCP<const shards::CellTopology> & cell_topo) {
    
    Teuchos::RCP<Intrepid::Basis<ScalarT,ArrayT> > basis;
    
    if (cell_dimension == 3) {
      if (type == "Q1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C1_FEM<ScalarT,ArrayT> );
      else if (type == "Q2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_HEX_C2_FEM<ScalarT,ArrayT> );
      else if (type == "QEdge1") 
	basis = Teuchos::rcp( new Intrepid::Basis_HCURL_HEX_I1_FEM<ScalarT,ArrayT> );
      else if (type == "T1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TET_C1_FEM<ScalarT,ArrayT> );
      else if (type == "T2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TET_C2_FEM<ScalarT,ArrayT> );
    }
    else if (cell_dimension == 2) {
      if (type == "Q1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_QUAD_C1_FEM<ScalarT,ArrayT> );
      else if (type == "Q2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_QUAD_C2_FEM<ScalarT,ArrayT> );
      else if (type == "QEdge1") 
	basis = Teuchos::rcp( new Intrepid::Basis_HCURL_QUAD_I1_FEM<ScalarT,ArrayT> );
      else if (type == "T1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TRI_C1_FEM<ScalarT,ArrayT> );
      else if (type == "T2")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_TRI_C2_FEM<ScalarT,ArrayT> );
      else if (type == "TEdge1") 
	basis = Teuchos::rcp( new Intrepid::Basis_HCURL_TRI_I1_FEM<ScalarT,ArrayT> );
    }
    else if (cell_dimension == 1) {
      if (type == "Q1")
	basis = Teuchos::rcp( new Intrepid::Basis_HGRAD_LINE_C1_FEM<ScalarT,ArrayT> );
    }


    if (Teuchos::is_null(basis)) {
      std::ostringstream s;
      s << "Failed to create basis.  Either the basis type in the input file, \"" << type << "\" is unsupported or we are out of memory!" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
			 s.str());
    }

    if((*cell_topo)!=basis->getBaseCellTopology()) {
      std::ostringstream s;
      s << "Failed to create basis.  Intrepid basis topology does not match mesh cell topology!";
      TEUCHOS_TEST_FOR_EXCEPTION((*cell_topo)!=basis->getBaseCellTopology(), std::runtime_error,
			 s.str());
    }

    return basis;
  }

  inline 
  std::string getBasisName(const std::string type) {
    const std::string name = "Basis " + type;
    return  name;
  }
  
  inline
  std::string getD1BasisName(const std::string type) {
    const std::string name = "Grad_Basis " + type;
    return  name;
  }

  inline
  std::string getD2BasisName(const std::string type) {
    const std::string name = "D2_Basis " + type;
    return  name;
  }
  
}


#endif
