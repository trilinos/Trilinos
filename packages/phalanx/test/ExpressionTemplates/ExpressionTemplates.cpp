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
#include "Teuchos_GlobalMPISession.hpp"

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

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

    ExprMult<int,double,ExprScalar<int,double>,ExprScalar<int,double> > 
      exm(a,b);
    
    TEST_FOR_EXCEPTION(exm.size() != 0, std::runtime_error,
		       "Wrong template instatiation selected!");

    const int vec_size = 100;
    ExprArray<std::size_t,double> ea1(vec_size);
    ExprArray<std::size_t,double> ea2(vec_size);
    ExprArray<std::size_t,double> ea4(vec_size);
    ExprArray<std::size_t,double> result(vec_size);

    for (int i=0; i < vec_size; ++i) {
      ea1[i] = 1.0;
      ea2[i] = 2.0;
      ea4[i] = 4.0;
    }

    // operator=
    cout << "Testing operator=()...";
    result = ea1;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 1.0 > 1.0e-12, std::runtime_error,
			 "Error operator= has failed!");
   
    cout << "Passed!" << endl;

    // ************
    // Addition
    // ************

    // ExprArray + ExprArray
    cout << "Testing \"ExprArray + ExprArray\"...";
    result = ea1 + ea2;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 3.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray + ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprScalar + ExprArray
    cout << "Testing \"ExprScalar + ExprArray\"...";
    result = 2.0 + ea1;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 3.0 > 1.0e-12, std::runtime_error,
			 "\"ExprScalar + ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprArray + ExprScalar
    cout << "Testing \"ExprArray + ExprScalar\"...";
    result = ea1 + 3.0;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 4.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray + ExprScalar\" has failed!");

    cout << "Passed!" << endl;

    // ************
    // Subtraction
    // ************

    // ExprArray - ExprArray
    cout << "Testing \"ExprArray - ExprArray\"...";
    result = ea2 - ea1;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 1.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray - ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprScalar - ExprArray
    cout << "Testing \"ExprScalar - ExprArray\"...";
    result = 2.0 - ea1;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 1.0 > 1.0e-12, std::runtime_error,
			 "\"ExprScalar - ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprArray - ExprScalar
    cout << "Testing \"ExprArray - ExprScalar\"...";
    result = ea1 - 3.0;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] + 2.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray - ExprScalar\" has failed!");

    cout << "Passed!" << endl;

    // ************
    // Multiplication
    // ************

    // ExprArray * ExprArray
    cout << "Testing \"ExprArray * ExprArray\"...";
    result = ea2 * ea1;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 2.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray * ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprScalar * ExprArray
    cout << "Testing \"ExprScalar * ExprArray\"...";
    result = 2.0 * ea1;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 2.0 > 1.0e-12, std::runtime_error,
			 "\"ExprScalar * ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprArray * ExprScalar
    cout << "Testing \"ExprArray * ExprScalar\"...";
    result = ea1 * 3.0;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 3.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray * ExprScalar\" has failed!");

    cout << "Passed!" << endl;

    // ************
    // Division
    // ************

    // ExprArray / ExprArray
    cout << "Testing \"ExprArray / ExprArray\"...";
    result = ea4 / ea2;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 2.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray / ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprScalar / ExprArray
    cout << "Testing \"ExprScalar / ExprArray\"...";
    result = 2.0 / ea2;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 1.0 > 1.0e-12, std::runtime_error,
			 "\"ExprScalar / ExprArray\" has failed!");

    cout << "Passed!" << endl;

    // ExprArray / ExprScalar
    cout << "Testing \"ExprArray / ExprScalar\"...";
    result = ea4 / 2.0;

    for (int i=0; i < vec_size; ++i)
      TEST_FOR_EXCEPTION(result[i] - 2.0 > 1.0e-12, std::runtime_error,
			 "\"ExprArray / ExprScalar\" has failed!");

    cout << "Passed!" << endl;


    // ************
    // Timings
    // ************
    const bool do_timing_comparison = false;
    if (do_timing_comparison) {
      RCP<Time> raw_time = TimeMonitor::getNewTimer("Raw Pointer Time");
      RCP<Time> expr_time = TimeMonitor::getNewTimer("Expression Time");


      const int my_size = 50000000;
      ExprArray<std::size_t,double> ea1(my_size);
      ExprArray<std::size_t,double> ea2(my_size);
      ExprArray<std::size_t,double> ea4(my_size);
      ExprArray<std::size_t,double> result(my_size);
      for (int i=0; i < my_size; ++i) {
	result[i] = 0.0;
	ea1[i] = 1.0;
	ea2[i] = 1.0;
      }

      ArrayRCP<double> raw_result = arcp<double>(my_size);
      ArrayRCP<double> raw_ea1 = arcp<double>(my_size);
      ArrayRCP<double> raw_ea2 = arcp<double>(my_size);
      for (int i=0; i < my_size; ++i) {
	raw_result[i] = 0.0;
	raw_ea1[i] = 1.0;
	raw_ea2[i] = 1.0;
      }

      double* r_ptr = raw_result.get();
      double* ea1_ptr = raw_ea1.get();
      double* ea2_ptr = raw_ea2.get();

      const int num_loops = 5;
      for (int loops = 0; loops < num_loops; ++loops) {
	
	{
	  TimeMonitor tm(*expr_time);
	  result = 4.0 * ea1 + 3.0 * ea2; 
	}
	
	{
	  TimeMonitor tm(*raw_time);
	  for (int i=0; i < my_size; ++i)
	    r_ptr[i] = 4.0 * ea1_ptr[i] + 3.0 * ea2_ptr[i];
	  	  
	}
      }

      const double t_expr = expr_time->totalElapsedTime();
      const double t_raw = raw_time->totalElapsedTime();
      cout << "\nt_expr / t_raw = "
	   << t_expr / t_raw << "% overhead!" << endl;
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
