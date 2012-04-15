/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_Graph_Impl_hpp_
#define _fei_Graph_Impl_hpp_

#include <fei_iosfwd.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Graph.hpp>
#include <fei_mpi.h>
#include <fei_EqnComm.hpp>

namespace fei {

  /** Basic algebraic matrix-graph representation. */
  class Graph_Impl : public fei::Graph {
  public:
    /** constructor */
    Graph_Impl(MPI_Comm comm, int firstLocalRow, int lastLocalRow);

    /** destructor */
    virtual ~Graph_Impl();

    /** Add indices to a specified row of the table */
    int addIndices(int row,
		   int len,
		   const int* indices);

    /** Add a symmetric block of indices. The array of indices will serve as
	both row-numbers, and as column-numbers in those rows.
    */
    int addSymmetricIndices(int numIndices,
			    int* indices,
			    bool diagonal=false);

    /** gather all remotely-owned table portions to owning processors */
    int gatherFromOverlap();

    /** Retrieve the local portion of the graph. i.e., The rows which
       correspond to locally-owned IDs. */
    table_type* getLocalGraph()
      {
	return( localGraphData_ );
      }

    /** Retrieve the remotely-owned portion of the graph. */
    std::vector<remote_table_type*>& getRemoteGraph()
      {
	return( remoteGraphData_);
      }

    /** Write locally-owned portion of the graph to a specified ostream.*/
    int writeLocalGraph(FEI_OSTREAM& os,
			bool debug=false,
			bool prefixLinesWithPoundSign=true);

    /** Write remotely-owned portion of the graph to a specified ostream.*/
    int writeRemoteGraph(FEI_OSTREAM& os);

    /** Get the number of locally-owned rows. */
    int getNumLocalRows();

    /** Get the number of nonzeros in the locally-owned portion of
       the graph. */
    int getNumLocalNonzeros() const;

    /** Get the length of a specified locally-owned row. */
    int getLocalRowLength(int row);

  private:
    void addDiagonals(int numIndices, int* indices);

    table_type* localGraphData_;
    std::vector<remote_table_type*> remoteGraphData_;
    fei::SharedPtr<fei::EqnComm> eqnComm_;

    int firstLocalRow_, lastLocalRow_;
    int localProc_, numProcs_;
    MPI_Comm comm_;
  };//class Graph_Impl

} //namespace fei

#endif // _fei_Graph_Impl_hpp_

