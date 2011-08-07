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
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Traits.hpp"

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
    // Start of debug naming testing
    // *********************************************************************
    {
      cout << "\nTesting Debug naming scheme.\n";
      // Evaluation Types
      cout << typeAsString<MyTraits::Residual>() << endl;
      cout << typeAsString<MyTraits::Jacobian>() << endl;
      // Data Types
      cout << typeAsString<double>() << endl;
      cout << typeAsString<MyTraits::FadType>() << endl;
      cout << typeAsString< MyVector<double> >() << endl;
      cout << typeAsString< MyVector<MyTraits::FadType> >() << endl;
      cout << typeAsString< MyTensor<double> >() << endl;
      cout << typeAsString< MyTensor<MyTraits::FadType> >() << endl;
    }
    
    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    cout << "\nRun has completed successfully!\n" << endl; 
    // *********************************************************************
    // *********************************************************************

  }
  catch (const exception& e) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Exception Caught!" << endl;
    cout << "Error message is below\n " << e.what() << endl;
    cout << "************************************************" << endl;
  }
  catch (...) {
    cout << "************************************************" << endl;
    cout << "************************************************" << endl;
    cout << "Unknown Exception Caught!" << endl;
    cout << "************************************************" << endl;
  }

  TimeMonitor::summarize();
    
  return 0;
}
