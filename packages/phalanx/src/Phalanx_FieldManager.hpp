// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_FIELD_MANAGER_HPP
#define PHX_FIELD_MANAGER_HPP

#include <cstddef>
#include <string>
#include <map>
#include <iostream>
#include <vector>
#include <algorithm>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_Field.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_EvaluationContainer_TemplateManager.hpp"

namespace PHX {

  template<typename Traits>
  class FieldManager {
    
  public:

    typedef typename PHX::EvaluationContainer_TemplateManager<Traits>::iterator iterator;

    FieldManager();

    ~FieldManager();
    
    void requireFieldForAllEvaluationTypes(const PHX::FieldTag& t);
    
    template<typename EvalT>
    void requireField(const PHX::FieldTag& t);

    void registerEvaluatorForAllEvaluationTypes(const Teuchos::RCP< PHX::Evaluator<Traits> >& e);
    
    template<typename EvalT>
    void registerEvaluator(const Teuchos::RCP< PHX::Evaluator<Traits> >& e);

    void registerEvaluator(typename PHX::FieldManager<Traits>::iterator it,
			   const Teuchos::RCP< PHX::Evaluator<Traits> >& e);
    
    template<typename DataT, typename EvalT> 
    void getFieldData(PHX::Field<DataT>& f);
    
    template<typename DataT, typename EvalT, 
	     typename Tag0, typename Tag1, typename Tag2, typename Tag3,
	     typename Tag4, typename Tag5, typename Tag6, typename Tag7> 
    void getFieldData(PHX::MDField<DataT,Tag0,Tag1,Tag2,Tag3,Tag4,Tag5,
		      Tag6,Tag7>& f);
    
    template<typename DataT, typename EvalT> 
    void getFieldData(const PHX::FieldTag& t, Teuchos::ArrayRCP<DataT>& d);
    
    //! Allocates memory for a single evaluation type
    template<typename EvalT>
    void postRegistrationSetupForType(typename Traits::SetupData d);

    //! Allocates memory for all evaluation types
    void postRegistrationSetup(typename Traits::SetupData d);

    template<typename EvalT>
    void evaluateFields(typename Traits::EvalData d);

    template<typename EvalT>
    void preEvaluate(typename Traits::PreEvalData d);

    template<typename EvalT>
    void postEvaluate(typename Traits::PostEvalData d);

    //! Return iterator to first EvaluationContainer
    typename FieldManager::iterator begin();

    //! Return iterator to last EvaluationContainer
    typename FieldManager::iterator end();

    //! Writes graphviz dot file for the evaluation type
    template<typename EvalT>
    void writeGraphvizFile(const std::string filename = "graph.dot",
			   bool writeEvaluatedFields = true,
			   bool writeDependentFields = false,
			   bool debugRegisteredEvaluators = false) const;

    //! Writes graphviz dot file for all evaluation types (adds eval type to filename).
    void writeGraphvizFile(const std::string base_filename = "graph",
			   const std::string file_extension = ".dot",
			   bool writeEvaluatedFields = true,
			   bool writeDependentFields = false,
			   bool debugRegisteredEvaluators = false) const;

    void print(std::ostream& os) const;

  private:

    typedef PHX::EvaluationContainer_TemplateManager<Traits> SCTM;

    std::size_t m_num_evaluation_types;

    PHX::EvaluationContainer_TemplateManager<Traits> m_eval_containers;

  };

  template<typename Traits>
  std::ostream& operator<<(std::ostream& os, 
			   const PHX::FieldManager<Traits>& vm);

} 

#include "Phalanx_FieldManager_Def.hpp"

#endif 
