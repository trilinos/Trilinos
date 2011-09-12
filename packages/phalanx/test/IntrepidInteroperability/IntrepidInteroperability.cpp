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
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

// Intrepid
#include "Intrepid_Types.hpp"
// #include "Intrepid_RealSpace.hpp"
// #include "Intrepid_FieldContainer.hpp"
#include "Intrepid_CubatureTensor.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"     // bvbw points to shards
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_Types.hpp"  
#include "Intrepid_Utils.hpp"  

#include "Shards_Array.hpp"

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Dim)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Dim)
  
SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(IP)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(IP)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(NODE)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(NODE)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Point)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Point)

SHARDS_ARRAY_DIM_TAG_SIMPLE_DECLARATION(Cell)
SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION(Cell)

int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  
  GlobalMPISession mpi_session(&argc, &argv);

  typedef MDField<double>::size_type size_type;

  try {
    
    RCP<Time> total_time = TimeMonitor::getNewTimer("Total Run Time");
    TimeMonitor tm(*total_time);

    // *********************************************************************
    // Start of Intrepid testing
    // *********************************************************************
    {
      // Inptrepid integration objects
      PHX::MDField<double,IP,Dim> cub_points;
      PHX::MDField<double,IP> cub_weights;
      PHX::MDField<double,Cell,NODE,Dim> node_coordinates;
      PHX::MDField<double,Cell,IP,Dim,Dim> jacobian;
      PHX::MDField<double,Cell,IP,Dim,Dim> jacobian_inv;
      PHX::MDField<double,Cell,IP> jacobian_det;
      PHX::MDField<double,Cell,IP> weighted_measure;
      // Diffusion terms
      PHX::MDField<double,Cell,IP,Dim> grad_of_basis_at_cub_points;
      PHX::MDField<double,Cell,NODE,IP,Dim> transformed_grad_of_basis_at_cub_points;
      PHX::MDField<double,Cell,NODE,IP,Dim> weighted_transformed_grad_of_basis_at_cub_points;
      // src terms
      PHX::MDField<double,NODE,IP> value_of_basis_at_cub_points;
      PHX::MDField<double,Cell,NODE,IP> transformed_value_of_basis_at_cub_points;
      PHX::MDField<double,Cell,NODE,IP> weighted_transformed_value_of_basis_at_cub_points;
      
      shards::CellTopology cellType = 
	shards::getCellTopologyData< shards::Quadrilateral<> >();
      
      Intrepid::DefaultCubatureFactory<double,MDField<double> > cubFactory;
      
      // Next few lines fail to compile due to ArrayType issues in intrepid....

      const int cubDegree = 2;
      
      Teuchos::RCP<Intrepid::Cubature<double,MDField<double> > > myCub = 
	cubFactory.create(cellType, cubDegree);
      
      //myCub->getCubature(PHX::MDField<double>(cub_points), PHX::MDField<double>(cub_weights));
      



      // Ignore what is below...

      //const int spaceDim = myCub->getDimension();
      //const int numCubPoints = myCub->getNumPoints(); 
      
      //Intrepid::Basis_HGRAD_QUAD_C1_FEM<double, Intrepid::FieldContainer<double> > quadBasis;
      
      //const int numNodes = quadBasis.getCardinality();
      
      // for (THashList::iterator cell = workset.begin; 
// 	   cell !=  workset.end; ++cell) {

//       }
      

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
