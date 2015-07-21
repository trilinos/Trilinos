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


#include "Phalanx_config.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Evaluator_Manager.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_DimTag.hpp"

// Evaluators
#include "evaluators/Evaluator_Constant.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

PHX_DIM_TAG_DECLARATION(Cell)
PHX_DIM_TAG_IMPLEMENTATION(Cell)

PHX_DIM_TAG_DECLARATION(Node)
PHX_DIM_TAG_IMPLEMENTATION(Node)

PHX_DIM_TAG_DECLARATION(QP)
PHX_DIM_TAG_IMPLEMENTATION(QP)

TEUCHOS_UNIT_TEST(kokkos, MemoryAssignment)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  EvaluatorManager<MyTraits> em;
  
  RCP<DataLayout> nodes = rcp(new MDALayout<Cell,Node>(100,4));
  RCP<DataLayout> qp = rcp(new MDALayout<Cell,QP>(100,4));
  
  RCP<FieldTag> den_n = rcp(new Tag<double>("Density", nodes));
  RCP<FieldTag> den_qp = rcp(new Tag<double>("Density", qp));
  
  em.requireField(*den_n);
  em.requireField(*den_qp);
  
  { 
    ParameterList p;
    p.set<string>("Name", "Density");
    p.set<double>("Value", 2.0);
    p.set< RCP<DataLayout> >("Data Layout", nodes);
    Teuchos::RCP< PHX::Evaluator<MyTraits> > ptr = 
      rcp(new Constant<MyTraits::Residual, MyTraits>(p));
    em.registerEvaluator(ptr);
  }
  { 
    ParameterList p;
    p.set<string>("Name", "Density");
    p.set<double>("Value", 2.0);
    p.set< RCP<DataLayout> >("Data Layout", qp);
    Teuchos::RCP< PHX::Evaluator<MyTraits> > ptr = 
      rcp(new Constant<MyTraits::Residual, MyTraits>(p));
    em.registerEvaluator(ptr);
  }
  
  em.setEvaluationTypeName(PHX::typeAsString<MyTraits::Residual>());

  em.sortAndOrderEvaluators();
  
  TEST_ASSERT(em.sortingCalled());

  const std::vector< Teuchos::RCP<PHX::FieldTag> >& tags = em.getFieldTags();

  TEST_EQUALITY(tags.size(),2);

  em.writeGraphvizFile("graph.dot",true,true,false);
  
  cout << "\n" << em << endl;  
}
