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

#ifndef PANZER_STK_SCATTER_CELL_AVG_VECTOR_DECL_HPP
#define PANZER_STK_SCATTER_CELL_AVG_VECTOR_DECL_HPP

#include "Phalanx_config.hpp"
#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"

#include "Teuchos_ParameterList.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_STK_Interface.hpp"

#include "Panzer_Evaluator_Macros.hpp"

namespace panzer_stk {

/** This class is a scatter operation to the mesh. It
  * takes a set of field names and an integration rule. 
  * Those quantities are components of vector fields. They are averaged over 
  * the cell and written to the mesh.
  *
  * The constructor takes a STK_Interface RCP and parameter list
  * that is required to contain the following fields:
  * "Scatter Name" string specifying the name of this evaulator
  * "Field Names" of type this is a comma seperated list of strings,
  * "IR" of type <code>Teuchos::RCP<panzer::IntegrationRule></code> and
  * "Mesh" of type <code>Teuchos::RCP<const panzer_stk::STK_Interface></code>.
  */
PANZER_EVALUATOR_CLASS(ScatterCellAvgVector)

  // typedef panzer_stk::STK_Interface::SolutionFieldType VariableField; // this is weird, but the correct thing
  typedef panzer_stk::STK_Interface::VectorFieldType VariableField;
  
  std::size_t numValues_;
 
  std::vector< PHX::MDField<const ScalarT,panzer::Cell,panzer::Point,panzer::Dim> > scatterFields_;
  Teuchos::RCP<STK_Interface> mesh_;

  std::vector<VariableField*> stkFields_;
 
PANZER_EVALUATOR_CLASS_END

}

// **************************************************************
#endif
