// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef PRINTDATA_HPP
#define PRINTDATA_HPP

#include "Zoltan2_config.h"
#include "Tpetra_CrsGraph.hpp"
#include "Teuchos_ArrayView.hpp"

#include <iostream>
#include <string>

using std::string;
using Teuchos::ArrayView;

template <typename lno_t, typename gno_t>
 void printTpetraGraph(const Tpetra::CrsGraph<lno_t, gno_t> &graph,
   std::ostream &os, size_t maxSize, string info)
{
  size_t nrows = graph.getNodeNumRows();
  if (nrows > maxSize)
    return;

  const RCP<const typename Tpetra::Map<lno_t, gno_t> > &rowMap=
    graph.getRowMap();
  const RCP<const typename Tpetra::Map<lno_t, gno_t> > &colMap=
    graph.getColMap();

  if (info.size() > 0)
    os << info << std::endl;

  if (graph.isGloballyIndexed()){
    ArrayView<const gno_t> indices;
    for (size_t i=0; i < nrows; i++){
      gno_t gid = rowMap->getGlobalElement(i);
      graph.getGlobalRowView(gid, indices);
      os << "Row " << gid << ": ";
      for (typename ArrayView<const gno_t>::size_type j=0; j < indices.size(); j++){
        os << indices[j] << " ";
      }
      os << std::endl;
    }
  }
  else{
    ArrayView<const lno_t> indices;
    for (size_t i=0; i < nrows; i++){
      gno_t gid = rowMap->getGlobalElement(i);
      graph.getLocalRowView(i, indices);
      os << "Row " << gid << ": ";
      for (typename ArrayView<const lno_t>::size_type j=0; j < indices.size(); j++){
        os << colMap->getGlobalElement(indices[j]) << " ";
      }
      os << std::endl;
    }
  }
}

template <typename lno_t, typename gno_t>
  void printTpetraGraph(const RCP<const Comm<int> > &comm,
  const Tpetra::CrsGraph<lno_t, gno_t> &graph, std::ostream &os,
  size_t maxSize, string info)
{
  int rank = comm->getRank();
  std::ostringstream oss;
  oss << "rank " << rank;

  comm->barrier();
  if (rank==0)
    os << info << std::endl;
  comm->barrier();

  for (int p=0; p < comm->getSize(); p++){
    if (p == rank)
      printTpetraGraph<lno_t, gno_t>(graph, os, maxSize, oss.str());
    comm->barrier();
  }
}

#endif
