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


#ifndef _fei_Graph_hpp_
#define _fei_Graph_hpp_

#include <fei_iosfwd.hpp>
#include <snl_fei_RaggedTable_specialize.hpp>

#include <map>

namespace fei {

  /** Basic algebraic matrix-graph representation. */
  class Graph {
  public:

    /** destructor */
    virtual ~Graph(){}

    /** alias for the 'table_type' data container */
    typedef snl_fei::RaggedTable<snl_fei::MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >
      table_type;

    /** alias for table_row_type, which is a row of the table */
    typedef fei::ctg_set<int> table_row_type;

    /** alias for the type of the remotely-owned portion of the table data */
    typedef snl_fei::RaggedTable<std::map<int,fei::ctg_set<int>*>,fei::ctg_set<int> > 
      remote_table_type;

    /** Add indices to a specified row of the table */
    virtual int addIndices(int row,
		   int len,
		   const int* indices) = 0;

    /** Add a symmetric block of indices. The array of indices will serve as
	both row-numbers, and as column-numbers in those rows.
    */
    virtual int addSymmetricIndices(int numIndices,
			    int* indices,
			    bool diagonal=false) = 0;

    /** gather all remotely-owned table portions to owning processors */
    virtual int gatherFromOverlap() = 0;

    /** Retrieve the local portion of the graph. i.e., The rows which
       correspond to locally-owned IDs. */
    virtual table_type* getLocalGraph() = 0;

    /** Retrieve the remotely-owned portion of the graph. */
    virtual std::vector<remote_table_type*>& getRemoteGraph() = 0;

    /** Write locally-owned portion of the graph to a specified ostream.*/
    virtual int writeLocalGraph(FEI_OSTREAM& os,
			bool debug=false,
			bool prefixLinesWithPoundSign=true) = 0;

    /** Write remotely-owned portion of the graph to a specified ostream.*/
    virtual int writeRemoteGraph(FEI_OSTREAM& os) = 0;

  };//class Graph

} //namespace fei

#endif // _fei_Graph_hpp_

