// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  size_t nrows = graph.getLocalNumRows();
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
