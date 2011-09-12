/*
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

*/

#include "Phalanx_ConfigDefs.hpp"
#include "Phalanx.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Dimension.hpp"
#include "DimTagCommon.hpp"
#include "Array.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Spatial : public phdmesh::DimTag {
  const char * name() const ;
  static const Spatial & descriptor();
};

struct Quadrature : public phdmesh::DimTag {
  const char * name() const ;
  static const Quadrature & descriptor();
};

struct Node : public phdmesh::DimTag {
  const char * name() const ;
  static const Node & descriptor();
};

struct Cell : public phdmesh::DimTag {
  const char * name() const ;
  static const Cell & descriptor();
};

const char * Spatial::name() const 
{ static const char n[] = "Spatial" ; return n ; }
const Spatial & Spatial::descriptor() 
{ static const Spatial myself ; return myself ; }

const char * Quadrature::name() const 
{ static const char n[] = "Quadrature" ; return n ; }
const Quadrature & Quadrature::descriptor() 
{ static const Quadrature myself ; return myself ; }

const char * Node::name() const 
{ static const char n[] = "Node" ; return n ; }
const Node & Node::descriptor() 
{ static const Node myself ; return myself ; }

const char * Cell::name() const 
{ static const char n[] = "Cell" ; return n ; }
const Cell & Cell::descriptor() 
{ static const Cell myself ; return myself ; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char *argv[]) 
{
  using namespace std;
  using namespace PHX;
  using namespace phdmesh;
  
  try {
    
    Teuchos::RCP<Teuchos::Time> total_time = 
      Teuchos::TimeMonitor::getNewTimer("Total Run Time");
    Teuchos::TimeMonitor tm(*total_time);
    
    {
      cout << "\nStarting MultiDimensionalArray Example!\n" << endl;
      
      // Simple Examples

      {
	// 1. Create a rank 0 tensor: i.e. a scalar source at each qp in mesh
	// 100 cells, 4 qp per cell
	DimNatural<Cell,Quadrature> dim_s(100, 4);
	double* memory = new double[dim_s.size()];
	Array<DimNatural<Cell,Quadrature>, double* > s(dim_s, memory);
	
	// Usage
	uint num_cells;
	uint num_qp;
	s.dimension().size(num_cells, num_qp);
	
	s(51,0) = 2.1;
	
	for (uint cell = 0; cell < num_cells; ++cell)
	  for (uint qp = 0; qp < num_qp; ++qp)
	    s(cell,qp) = 2.0;

	// This is ugly - need to add bracket operator
	for (uint i = 0; i < s.dimension().size(); ++i)
	  s.data()[i] = 2.0;

	delete [] memory;
      }

      {
	// 2. Create a rank 1 tensor: i.e. a vector at each qp in mesh
	// 100 cells, 4 qp per cell, 3 entries per qp
	DimNatural<Cell,Quadrature,Spatial> dim_v(100, 4, 3);
	double* memory = new double[dim_v.size()];
	Array<DimNatural<Cell,Quadrature,Spatial>, double* > v(dim_v, memory);
	
	// Usage
	uint num_cells;
	uint num_qp;
	uint num_sp_dim;
	v.dimension().size(num_cells, num_qp, num_sp_dim);
	
	v(51,0,1) = 3.1;

	for (uint cell = 0; cell < num_cells; ++cell)
	  for (uint qp = 0; qp < num_qp; ++qp)
	    for (uint dim = 0; dim < num_sp_dim; ++dim)
	      v(cell,qp, dim) = 3.0;
	
	delete [] memory;
      }

      {
	// 3. Create a rank 2 tensor: i.e. a matrix at each qp in a mesh
	// 100 cells, 4 qp per cell, 3X3 matrix = 9 entries per qp
	DimNatural<Cell,Quadrature,Spatial, Spatial> dim_m(100, 4, 3, 3);
	
	double* memory = new double[dim_m.size()];
	
	Array<DimNatural<Cell,Quadrature,Spatial,Spatial>, double* > 
	  A(dim_m, memory);
	
	// Usage
	uint num_cells;
	uint num_qp;
	uint num_col;
	uint num_row;
	A.dimension().size(num_cells, num_qp, num_col, num_row);
	
	A(51,0,1,1) = 4.1;

	for (uint cell = 0; cell < num_cells; ++cell)
	  for (uint qp = 0; qp < num_qp; ++qp)
	    for (uint col = 0; col < num_col; ++col)
	      for (uint row = 0; row < num_row; ++row)
		A(cell,qp, col, row) = 4.0;
	
	

	delete [] memory;
      }

      // 

      typedef DimNatural<Cell,Quadrature,Spatial,Spatial> QP_Matrix;
      typedef DimNatural<Cell,Quadrature,Spatial> QP_Vector;

      typedef DimNatural<Cell,Node,Spatial,Spatial> N_Matrix;
      typedef DimNatural<Cell,Node,Spatial> N_Vector;

      const int num_cells = 100;
      const int num_qp = 8;

      QP_Vector DimV3_QP(num_cells, num_qp, 3);
      QP_Matrix DimM3_QP(num_cells, num_qp, 3, 3);

      cout << "DimV3_QP = " << DimV3_QP << endl;
      cout << "DimM3_QP = " << DimM3_QP << endl;

      int num_vectors = 3;
      int num_matrices = 1;
      int total_memory = num_vectors * num_cells * num_qp * 3 + 
	num_matrices * num_cells * num_qp * 3 * 3;
      double* memory_pool = new double[total_memory];

      Array<QP_Vector,double*> x(DimV3_QP, &memory_pool[0]); 
      Array<QP_Vector,double*> b(DimV3_QP, &memory_pool[DimV3_QP.size()]); 
      Array<QP_Vector,double*> A(DimV3_QP, &memory_pool[2 *DimV3_QP.size()]); 
      


      // Evaluation function
      {
	uint num_cell;
	uint num_qp;
	uint num_dim;
	x.dimension().size(num_cell, num_qp, num_dim);
	
	for (uint cell=0; cell < num_cell; ++cell)
	  for (uint qp=0; qp < num_qp; ++qp)
	    for (uint dim=0; dim < num_dim; ++dim) {
	      //cout << "cell = " << cell << ", qp = " << qp 
	      //   << ", dim = " << dim << endl;
	      x(cell, qp, dim) = 0.0;
	    }
      }

      delete memory_pool;

    }
    
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

  Teuchos::TimeMonitor::summarize();
    
  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
