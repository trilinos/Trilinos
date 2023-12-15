// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_LWGRAPH_DEF_HPP
#define MUELU_LWGRAPH_DEF_HPP

#include "MueLu_LWGraph_decl.hpp"
#include <Xpetra_CrsGraph.hpp>
#include <Xpetra_CrsGraphFactory.hpp>

namespace MueLu {

//! Print the object with some verbosity level to an FancyOStream object.
// using MueLu::Describable::describe; // overloading, not hiding
// void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const {
template <class LocalOrdinal, class GlobalOrdinal, class Node>
void LWGraph<LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
  //    MUELU_DESCRIBE;

  if (verbLevel & Parameters0) {
    // out0 << "Prec. type: " << type_ << std::endl;
  }

  if (verbLevel & Parameters1) {
    // out0 << "Linear Algebra: " << toString(lib_) << std::endl;
    // out0 << "PrecType: " << type_ << std::endl;
    // out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
    // out0 << "Overlap: " << overlap_ << std::endl;
  }

  if (verbLevel & Debug) {
    RCP<const Map> col_map = importMap_.is_null() ? domainMap_ : importMap_;

    for (LO i = 0; i < rows_.size() - 1; i++) {
      for (LO j = rows_[i]; j < rows_[i + 1]; j++)
        out << domainMap_->getGlobalElement(i) << " " << col_map->getGlobalElement(columns_[j]) << std::endl;
    }
  }
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > LWGraph<LocalOrdinal, GlobalOrdinal, Node>::GetCrsGraph() const {
  ArrayRCP<size_t> rowPtrs;
  rowPtrs.resize(rows_.size());
  for (size_t i = 0; i < Teuchos::as<size_t>(rows_.size()); i++)
    rowPtrs[i] = rows_[i];
  auto graph = Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(GetDomainMap(), GetImportMap(), rowPtrs, Teuchos::arcp_const_cast<LO>(getEntries()));
  graph->fillComplete();
  return graph;
}

}  // namespace MueLu

#endif  // MUELU_LWGRAPH_DEF_HPP
