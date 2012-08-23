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

#ifndef PANZER_EQUATION_SET_FACTORY_HPP
#define PANZER_EQUATION_SET_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"

namespace panzer {

  /** \brief Allocates and initializes an equation set template manager

     \param ies [in] Input equation set
     \param cell_data [in] The cell data
     \param global_data [in] Global data
     \param build_transient_support [in] If true, the transient evaluators will be built, registered, and required in the Phalanx evaluation graph.

     Returns an RCP to a newly allocated EquationSet_TemplateManager.  
  */
  struct EquationSetFactory {
    virtual ~EquationSetFactory() = 0;

    virtual Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const panzer::InputEquationSet& ies,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     bool build_transient_support) const = 0;

  };

  inline EquationSetFactory::~EquationSetFactory() {}
  
}

#endif
