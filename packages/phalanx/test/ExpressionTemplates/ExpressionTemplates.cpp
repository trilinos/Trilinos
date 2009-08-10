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
#include "Phalanx_ExpressionTemplates_Traits.hpp"
#include "Phalanx_ExpressionTemplates_Operands.hpp"
#include "Phalanx_ExpressionTemplates_Operators.hpp"
#include "Phalanx_ExpressionTemplates_Array.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    
    double a_value = 1.0;
    double b_value = 1.0;
    ExprScalar<int,double> a(a_value);
    ExprScalar<int,double> b(b_value);

    ExprAdd<int,double,ExprScalar<int,double>,ExprScalar<int,double> > 
      exa(a,b);

    TEST_FOR_EXCEPTION(exa.size() != 0, std::runtime_error,
		       "Wrong template instatiation selected!");

    cout << "exa.size() = " << exa.size() << endl;

    ExprMult<int,double,ExprScalar<int,double>,ExprScalar<int,double> > 
      exm(a,b);
    
    TEST_FOR_EXCEPTION(exm.size() != 0, std::runtime_error,
		       "Wrong template instatiation selected!");

    cout << "exm.size() = " << exm.size() << endl;

    const int vec_size = 1000000;
    ExprArray<std::size_t,double> ea1(vec_size);
    ExprArray<std::size_t,double> ea2(vec_size);
    ExprArray<std::size_t,double> ea3(vec_size);
    ExprArray<std::size_t,double> ea4(vec_size);
    ExprArray<std::size_t,double> ea5(vec_size);

    for (int i=0; i < ea1.size(); ++i) {
      ea1[i] = 1.0;
      ea2[i] = 2.0;
      ea3[i] = 3.0;
      ea4[i] = 4.0;
      ea5[i] = 5.0;
    }

    ea1 = ea2;

    for (int i=0; i < ea1.size(); ++i)
      TEST_FOR_EXCEPTION(ea1[i] - ea2[i] > 1.0e-12, std::runtime_error,
			 "Error operator= has failed!");
    
    ea1 = ea2 + ea3;

    // for (int i=0; i < ea1.size(); ++i)
//       cout << ea1[i] << endl;


    // Timings

    // Expression templates
    RCP<Time> expr_time = TimeMonitor::getNewTimer("Expr. Templates Time");
    {
      TimeMonitor tm(*expr_time);
      ea1 = ea2 + ea3 + ea4 + ea5;
    }

    // for loops
    RCP<Time> loop_time = TimeMonitor::getNewTimer("Hand Coded Loops Time");
    {
      TimeMonitor tm(*loop_time);
      for (int i=0; i < ea1.size(); ++i)
	ea1[i] = ea2[i] + ea3[i] + ea4[i] + ea5[i];
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
