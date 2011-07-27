// @HEADER
// ************************************************************************
// 
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                  Copyright (2008) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_Evaluator_Manager.hpp"
#include "Phalanx_TypeStrings.hpp"

// Evaluators
#include "evaluators/Evaluator_Constant.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// From test/Utilities directory
#include "Traits.hpp"

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Node)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Node)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(QP)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(QP)

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of EvaluatorManager Testing
    // *********************************************************************
    {
      cout << "\nStarting EvaluatorManager Testing\n";

      
      cout << "\nConstructing EvaluatorManager...";
      EvaluatorManager<MyTraits> em;
      cout << "Passed!" << endl;

      RCP<DataLayout> nodes = rcp(new MDALayout<Cell,Node>(100,4));
      RCP<DataLayout> qp = rcp(new MDALayout<Cell,QP>(100,4));

      RCP<FieldTag> den_n = rcp(new Tag<double>("Density", nodes));
      RCP<FieldTag> den_qp = rcp(new Tag<double>("Density", qp));
      
      
      cout << "\nTesting requireField()...";
      em.requireField(*den_n);
      em.requireField(*den_qp);
      cout << "Passed!" << endl;

      cout << "\nTesting registerEvaluator()...";
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
      cout << "Passed!" << endl;
      
      cout << "\nTesting setEvaluationTypeName()...";
      em.setEvaluationTypeName(PHX::typeAsString<MyTraits::Residual>());
      cout << "Passed!" << endl;


      cout << "\nTesting sortAndOrderEvaluators()...";
      em.sortAndOrderEvaluators();
      cout << "Passed!" << endl;

      cout << "\nTesting writeGraphvizFile()...";
      em.writeGraphvizFile("graph.dot",true,true,false);
      cout << "Passed!" << endl;

      cout << "\nPrinting EvaluatorManager:\n";
      cout << em << endl;

    }
    
    // *********************************************************************
    // *********************************************************************
    std::cout << "\nTest passed!\n" << std::endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const std::exception& e) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Exception Caught!" << endl;
    std::cout << "Error message is below\n " << e.what() << endl;
    std::cout << "************************************************" << endl;
  }
  catch (...) {
    std::cout << "************************************************" << endl;
    std::cout << "************************************************" << endl;
    std::cout << "Unknown Exception Caught!" << endl;
    std::cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
