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


#include "Phalanx_ConfigDefs.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Element_Linear2D.hpp"
#include "MeshBuilder.hpp"
#include "Epetra_SerialComm.h"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  
  Teuchos::GlobalMPISession mpi_session(&argc, &argv);
  
  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);
    
    cout << "\nStarting Element_Linear2D Testing:" << endl;
 
    cout << "Testing area integration..." << endl;
    {
      // A Trapezoid centered at (0,0)
      vector<double> xc(4);
      xc[0] = -1.0;
      xc[1] =  1.0;
      xc[2] =  0.5;
      xc[3] = -0.5;
      vector<double> yc(4);
      yc[0] = -1.0;
      yc[1] = -1.0;
      yc[2] = 1.0;
      yc[3] = 1.0;

      vector<int> gid(4);
      gid[0] = 0;
      gid[1] = 1;
      gid[2] = 2;
      gid[3] = 3;

      Element_Linear2D e(gid, 0, 0, xc,yc);
      
      // Integrate the area to test jacobian transform
      double area = 0.0;
      for (int qp=0; qp < e.numQuadraturePoints(); ++qp) {
	area += e.quadratureWeights()[qp] * e.detJacobian()[qp];
	cout << "|J| = " << e.detJacobian()[qp] << endl;
	
      }
      cout.setf(ios::scientific);
      cout.precision(8);
      cout << "    Area = " << area << " (should be 3.0)" << endl;
      cout << "    Difference = " << fabs(area - 3.0) << " (should be 0.0)" << endl;
      TEST_FOR_EXCEPTION((fabs(area - 3.0) > 1e-8), std::logic_error, 
			 "Area integration failed!");
      cout << "    **Area integration...passed!" << endl;
    }
    

    cout << "\nTesting gradient interpolation..." << endl;
    {
      // A Square centered at (0,0)
      vector<double> xc(4);
      xc[0] = -1.0;
      xc[1] =  1.0;
      xc[2] =  1.0;
      xc[3] = -1.0;
      vector<double> yc(4);
      yc[0] = -1.0;
      yc[1] = -1.0;
      yc[2] =  1.0;
      yc[3] =  1.0;
      
      vector<int> gid(4);
      gid[0] = 0;
      gid[1] = 1;
      gid[2] = 2;
      gid[3] = 3;
      
      Element_Linear2D e(gid, 0, 0, xc,yc);
      
      std::vector<double> u(4);
      // values give a slope of 1.0 in x direction
      // values give a slope of 0.0 in y direction
      u[0] = 0.0;
      u[1] = 2.0;
      u[2] = 2.0;
      u[3] = 0.0;
      
      const shards::Array<double,shards::NaturalOrder,QuadPoint,Node,Dim>& 
	dphi = e.basisFunctionGradientsRealSpace();
      
      for (int qp=0; qp < e.numQuadraturePoints(); ++qp) {
	
	double dudx = 0.0;
	double dudy = 0.0;
	
	for (int node=0; node < e.numNodes(); ++node) {
	  dudx += u[node] * dphi(qp,node,0);
	  dudy += u[node] * dphi(qp,node,1);
	}
	
	cout << "    dudx = " << dudx << endl;
	cout << "    dudy = " << dudy << endl;
	
	TEST_FOR_EXCEPTION( fabs(dudx - 1.0) > 1.0e-8, std::logic_error,
			    "X Derivative is incorrect!" );
	TEST_FOR_EXCEPTION( fabs(dudy) > 1.0e-8, std::logic_error,
			    "Y Derivative is incorrect!" );
      }
      
      // values give a slope of 0.0 in x direction
      // values give a slope of 1.0 in y direction
      u[0] = 0.0;
      u[1] = 0.0;
      u[2] = 2.0;
      u[3] = 2.0;
      for (int qp=0; qp < e.numQuadraturePoints(); ++qp) {
	
	double dudx = 0.0;
	double dudy = 0.0;
	
	for (int node=0; node < e.numNodes(); ++node) {
	  dudx += u[node] * dphi(qp,node,0);
	  dudy += u[node] * dphi(qp,node,1);
	}
	
	cout << "    dudx = " << dudx << endl;
	cout << "    dudy = " << dudy << endl;
	
	TEST_FOR_EXCEPTION( fabs(dudx) > 1.0e-8, std::logic_error,
			    "X Derivative is incorrect!" );
	TEST_FOR_EXCEPTION( fabs(dudy - 1.0) > 1.0e-8, std::logic_error,
			    "Y Derivative is incorrect!" );
	
      }
    }
    cout << "  **Gradient interpolation passed!" << endl;
    
    cout << "\nTesting mesh integration..." << endl;
    {
      RCP<Epetra_Comm> comm = rcp(new Epetra_SerialComm);
      MeshBuilder mb(comm, 100, 200, 3.0, 5.0, 8);
      
      std::vector<Element_Linear2D>& cells = *(mb.myElements());

      double area = 0.0;
      for (std::vector<Element_Linear2D>::iterator cell = cells.begin();
	   cell != cells.end(); ++cell) {
	for (int qp = 0; qp < cell->numQuadraturePoints(); ++qp) {
	  area += 
	    cell->quadratureWeights()[qp] * cell->detJacobian()[qp];
	}
      }
      cout << "    area = " << area << ", should be " << 15.0 << endl;
      TEST_FOR_EXCEPTION( ((area - 15.0) > 1.0e-12), std::logic_error,
			  "Mesh area integration failed!");
    }
    cout << "  **Mesh Integration passed!" << endl;
    

    // *********************************************************************
    // Finished all testing
    // *********************************************************************
    std::cout << "\nRun has completed successfully!\n" << std::endl; 
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

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
