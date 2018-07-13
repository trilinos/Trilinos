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

#include "Teuchos_RCP.hpp"

#include "Panzer_Base.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_LinearObjContainer.hpp"

namespace panzer {
  class FieldManagerBuilder;
  class AssemblyEngineInArgs;
}

namespace panzer {

  //! Class for the matrix and residual fill.
  template <typename EvalT>
    class AssemblyEngine : public panzer::Base {

  public:    
    struct EvaluationFlags {
      EvaluationFlags(int flags) : value_(flags) {
        TEUCHOS_ASSERT(flags>0 && flags <= EvaluationFlags::All);
      }
      static constexpr int Initialize=1;
      static constexpr int VolumetricFill=2;
      static constexpr int BoundaryFill=4;
      static constexpr int Scatter=8;
      static constexpr int All=15;
      int getValue() const {return value_;}
    protected:
      int value_;
    };

    AssemblyEngine(const Teuchos::RCP<panzer::FieldManagerBuilder>& fmb,
                   const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > & lof);
    
    void evaluate(const panzer::AssemblyEngineInArgs& input_arguments, const EvaluationFlags flags=EvaluationFlags(EvaluationFlags::All));

    void evaluateVolume(const panzer::AssemblyEngineInArgs& input_arguments);

    /** This method returns the global counter used to indicate which rows are boundary conditions.
      * Note that this method does all the communication neccessary to evaluate the dirichlet boundary
      * conditions. The dirichlet values are set in the global "F" vector, and the count values are
      * in the return linear obj containers "X" vector.
      */
    Teuchos::RCP<LinearObjContainer> evaluateOnlyDirichletBCs(const panzer::AssemblyEngineInArgs& input_arguments);

    void evaluateNeumannBCs(const panzer::AssemblyEngineInArgs& input_arguments);

    void evaluateInterfaceBCs(const panzer::AssemblyEngineInArgs& input_arguments);

    //! This method returns the global counter used to indicate which rows are boundary conditions
    Teuchos::RCP<LinearObjContainer> evaluateDirichletBCs(const panzer::AssemblyEngineInArgs& input_arguments);

    Teuchos::RCP<panzer::FieldManagerBuilder> getManagerBuilder()
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
    
      Teuchos::RCP<panzer::FieldManagerBuilder> m_field_manager_builder;

      Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > m_lin_obj_factory;

      // These members improve performance by ensuring that "buildPrimitiveGhostedLinearObjContainer"
      // is not called uneccessarily often
      bool countersInitialized_;
      Teuchos::RCP<LinearObjContainer> localCounter_;
      Teuchos::RCP<LinearObjContainer> globalCounter_;
      Teuchos::RCP<LinearObjContainer> summedGhostedCounter_;
  };
  
}

// #include "Panzer_AssemblyEngine_impl.hpp"

#endif
