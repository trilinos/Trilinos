// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
#ifndef TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_DP_HPP
#define TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_DP_HPP

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>

#include "fem_assembly_typedefs.hpp"
#include "fem_assembly_MeshDatabase.hpp"
#include "fem_assembly_Element.hpp"
#include "fem_assembly_utility.hpp"
#include "fem_assembly_commandLineOpts.hpp"

namespace TpetraExamples
{

using comm_ptr_t = Teuchos::RCP<const Teuchos::Comm<int> >;



int executeTotalElementLoopDP(const comm_ptr_t& comm, const struct CmdLineOpts& opts)
{
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  const global_ordinal_t GO_INVALID = Teuchos::OrdinalTraits<global_ordinal_t>::invalid();

  // The output stream 'out' will ignore any output not from Process 0.
  RCP<Teuchos::FancyOStream> pOut = getOutputStream(*comm);
  Teuchos::FancyOStream& out = *pOut;  

  out << "================================================================================" << std::endl
      << "=  Total Element Loop (Dynamic Profile)"    << std::endl
      << "================================================================================" << std::endl
      << std::endl;

  // Processor decomp (only works on perfect squares)
  int numProcs  = comm->getSize();
  int sqrtProcs = sqrt(numProcs); 

  if(sqrtProcs*sqrtProcs != numProcs) 
  {
    if(0 == comm->getRank())
      std::cerr << "Error: Invalid number of processors provided, num processors must be a perfect square." << std::endl;
    return -1;
  }
  int procx = sqrtProcs;
  int procy = sqrtProcs;

  // Generate a simple 3x3 mesh
  int nex = opts.numElementsX;
  int ney = opts.numElementsY;

  MeshDatabase mesh(comm,nex,ney,procx,procy);

  if(opts.verbose) mesh.print(std::cout);

  // Build Tpetra Maps
  // -----------------
  // -- https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1Map.html#a24490b938e94f8d4f31b6c0e4fc0ff77
  RCP<const map_t> row_map = rcp(new map_t(GO_INVALID, mesh.getOwnedNodeGlobalIDs(), 0, comm));

  if(opts.verbose) row_map->describe(out);

  // Graph Construction
  // ------------------
  // - Loop over every element in the mesh.
  //   - Get list of nodes associated with each element.
  //   - Insert node contributions if the node (row) is owned locally.
  //   
  auto domain_map = row_map;
  auto range_map  = row_map;

  auto owned_element_to_node_ids = mesh.getOwnedElementToNode();
  auto ghost_element_to_node_ids = mesh.getGhostElementToNode();

  RCP<TimeMonitor> timerGlobal = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("X) Global")));
  RCP<TimeMonitor> timerElementLoopGraph = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("1) ElementLoop  (Graph)")));

  RCP<graph_t> crs_graph = rcp(new graph_t(row_map, 0));
  
  // Using 4 because we're using quads for this example, so there will be 4 nodes associated with each element.
  Teuchos::Array<global_ordinal_t> global_ids_in_row(4);

  // Insert node contributions for every OWNED element:
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
    // Populate global_ids_in_row:
    // - Copy the global node ids for current owned element into an array.  
    // - Since each element's contribution is a clique, we can re-use this for 
    //   each row associated with this element's contribution.
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
      global_ids_in_row[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
    }

    // Add the contributions from the current row into the graph if the node is owned.
    // - For example, if Element 0 contains nodes [0,1,4,5] and nodes 0 and 4 are owned
    //   by the current processor, then:
    //   - node 0 inserts [0, 1, 4, 5]
    //   - node 1 <skip>
    //   - node 4 inserts [0, 1, 4, 5]
    //   - node 5 <skip>
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
      if(mesh.nodeIsOwned(global_ids_in_row[element_node_idx]))
      {
       crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
      }
    }
  }

  // Insert the node contributions for every GHOST element:
  for(size_t element_gidx=0; element_gidx<mesh.getNumGhostElements(); element_gidx++)
  {
    for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
    {
      global_ids_in_row[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
    }
    for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
    {
      if(mesh.nodeIsOwned(global_ids_in_row[element_node_idx]))
      {
       crs_graph->insertGlobalIndices(global_ids_in_row[element_node_idx], global_ids_in_row());
      }
    }
  }

  timerElementLoopGraph = Teuchos::null;

  // 'finalize' the crs_graph by calling fillComplete().
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("2) FillComplete (Graph)"));
    crs_graph->fillComplete();
  }

  // Print out the crs_graph in detail...
  if(opts.verbose) crs_graph->describe(out, Teuchos::VERB_EXTREME);

  // Matrix Fill
  // -------------------
  // In this example, we're using a simple stencil of values for the stiffness matrix:
  //
  //    +-----+-----+-----+-----+
  //    |  2  | -1  |     | -1  |
  //    +-----+-----+-----+-----+
  //    | -1  |  2  | -1  |     |
  //    +-----+-----+-----+-----+
  //    |     | -1  |  2  | -1  |
  //    +-----+-----+-----+-----+
  //    | -1  |     | -1  |  2  |
  //    +-----+-----+-----+-----+
  //
  // For matrix fill, we create the crs_matrix object and will fill it 
  // in the same manner as we filled in the graph but in this case, nodes 
  // associated with each element will receive contributions according to 
  // the row in this stencil.
  //
  // In this example, the calls to sumIntoGlobalValues() on 1 core will look like:
  //   Element 0
  // - sumIntoGlobalValues( 0,  [  0  1  5  4  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 1,  [  0  1  5  4  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 5,  [  0  1  5  4  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 4,  [  0  1  5  4  ],  [  -1  0  -1  2  ])
  // Element 1
  // - sumIntoGlobalValues( 1,  [  1  2  6  5  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 2,  [  1  2  6  5  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 6,  [  1  2  6  5  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 5,  [  1  2  6  5  ],  [  -1  0  -1  2  ])
  // Element 2
  // - sumIntoGlobalValues( 2,  [  2  3  7  6  ],  [  2  -1  0  -1  ])
  // - sumIntoGlobalValues( 3,  [  2  3  7  6  ],  [  -1  2  -1  0  ])
  // - sumIntoGlobalValues( 7,  [  2  3  7  6  ],  [  0  -1  2  -1  ])
  // - sumIntoGlobalValues( 6,  [  2  3  7  6  ],  [  -1  0  -1  2  ])
  //
  // Similarly to the Graph construction above, we loop over both local and global
  // elements and insert rows for only the locally owned rows.
  //
  RCP<TimeMonitor> timerElementLoopMatrix = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("3) ElementLoop  (Matrix)")));

  RCP<matrix_t> crs_matrix = rcp(new matrix_t(crs_graph));

  scalar_2d_array_t element_matrix;
  Kokkos::resize(element_matrix, 4, 4);

  Teuchos::Array<global_ordinal_t> column_global_ids(4);     // global column ids list
  Teuchos::Array<Scalar> column_scalar_values(4);         // scalar values for each column

  // Loop over owned elements:
  for(size_t element_gidx=0; element_gidx<mesh.getNumOwnedElements(); element_gidx++)
  {
    // Get the stiffness matrix for this element
    ReferenceQuad4(element_matrix);

    // Fill the global column ids array for this element
    for(size_t element_node_idx=0; element_node_idx<owned_element_to_node_ids.extent(1); element_node_idx++)
    {
      column_global_ids[element_node_idx] = owned_element_to_node_ids(element_gidx, element_node_idx);
    }
    
    // For each node (row) on the current element:
    // - populate the values array
    // - add values to crs_matrix if the row is owned.
    //   Note: hardcoded 4 here because we're using quads.
    for(size_t element_node_idx=0; element_node_idx<4; element_node_idx++)
    { 
      global_ordinal_t global_row_id = owned_element_to_node_ids(element_gidx, element_node_idx);
      if(mesh.nodeIsOwned(global_row_id)) 
      {
        for(size_t col_idx=0; col_idx<4; col_idx++)
        {
          column_scalar_values[col_idx] = element_matrix(element_node_idx, col_idx);
        }
        crs_matrix->sumIntoGlobalValues(global_row_id, column_global_ids, column_scalar_values);
      }
    }
  }

  // Loop over ghost elements:
  // - This loop is the same as the element loop for owned elements, but this one
  //   is for ghost elements.
  for(size_t element_gidx=0; element_gidx<mesh.getNumGhostElements(); element_gidx++)
  {
    ReferenceQuad4(element_matrix);

    for(size_t element_node_idx=0; element_node_idx<ghost_element_to_node_ids.extent(1); element_node_idx++)
    {
      column_global_ids[element_node_idx] = ghost_element_to_node_ids(element_gidx, element_node_idx);
    }
    
    for(size_t element_node_idx=0; element_node_idx<4; element_node_idx++)
    { 
      global_ordinal_t global_row_id = ghost_element_to_node_ids(element_gidx, element_node_idx);
      if(mesh.nodeIsOwned(global_row_id)) 
      {
        for(size_t col_idx=0; col_idx<4; col_idx++)
        {
          column_scalar_values[col_idx] = element_matrix(element_node_idx, col_idx);
        }
        crs_matrix->sumIntoGlobalValues(global_row_id, column_global_ids, column_scalar_values);
      }
    }
  }
  timerElementLoopMatrix = Teuchos::null;

  // After the contributions are added, 'finalize' the matrix using fillComplete()
  {
    TimeMonitor timer(*TimeMonitor::getNewTimer("4) FillComplete (Matrix)"));
    crs_matrix->fillComplete();
  }

  timerGlobal = Teuchos::null;

  // Print out crs_matrix details.
  if(opts.verbose) crs_matrix->describe(out, Teuchos::VERB_EXTREME);

  // Save crs_matrix as a MatrixMarket file.
  if(opts.saveMM)
  {
    std::ofstream ofs("crsMatrix_TotalElementLoop_DP.out", std::ofstream::out);
    Tpetra::MatrixMarket::Writer<matrix_t>::writeSparse(ofs, crs_matrix);
    ofs.close();
  }

  return 0;
}


} // namespace TpetraExamples


#endif // TPETRAEXAMPLES_FEM_ASSEMBLY_TOTALELEMENTLOOP_DP_HPP
