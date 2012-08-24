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

#ifndef PANZER_ASSEMBLY_ENGINE_HPP
#define PANZER_ASSEMBLY_ENGINE_HPP

#include "Panzer_Base.hpp"
#include "Teuchos_RCP.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"

namespace panzer {
  template <typename LO, typename GO> class FieldManagerBuilder;
  struct AssemblyEngineInArgs;
}

namespace panzer {

  //! Class for the matrix and residual fill.
  template <typename EvalT, typename LO, typename GO>
    class AssemblyEngine : public panzer::Base {

  public:    
    
    AssemblyEngine(const Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> >& fmb,
                   const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof);
    
    void evaluate(const panzer::AssemblyEngineInArgs& input_arguments);

    void evaluateVolume(const panzer::AssemblyEngineInArgs& input_arguments);
    void evaluateNeumannBCs(const panzer::AssemblyEngineInArgs& input_arguments);
    void evaluateDirichletBCs(const panzer::AssemblyEngineInArgs& input_arguments);

    Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> > getManagerBuilder()
      { return m_field_manager_builder; }
    
  protected:
      
    /** Evaluate both Dirichlet and Neumann conditions.
      *
      * \param[in] bc_type Type of Dirichlet condition to evaluate
      * \param[in] input_arguments Get solver parameters (alpha,beta, linear object containers)
      * \param[in] preEval_loc Linear object container used by Dirichlet conditions for
      *                        keeping track of rows that have been modified.
      */
    void evaluateBCs(const panzer::BCType bc_type, 
		     const panzer::AssemblyEngineInArgs& input_arguments,
                     const Teuchos::RCP<LinearObjContainer> preEval_loc=Teuchos::null);

  protected:
    
      Teuchos::RCP<panzer::FieldManagerBuilder<LO,GO> > 
      m_field_manager_builder;

      Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > 
      m_lin_obj_factory;
    
  };
  
}

#include "Panzer_AssemblyEngine_impl.hpp"

#endif
