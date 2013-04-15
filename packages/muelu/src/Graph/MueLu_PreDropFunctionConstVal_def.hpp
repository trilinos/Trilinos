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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
#define MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP

#include <Xpetra_CrsGraphFactory.hpp>

#include "MueLu_PreDropFunctionConstVal_decl.hpp"
#include "MueLu_Graph.hpp"
#include "Teuchos_ScalarTraits.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PreDropFunctionConstVal(const Scalar threshold)
    : threshold_(threshold) { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Drop(size_t lrow, GlobalOrdinal grow, size_t k, LocalOrdinal lcid, GlobalOrdinal gcid, const Teuchos::ArrayView<const LocalOrdinal> & indices, const Teuchos::ArrayView<const Scalar> & vals) {
    if(Teuchos::ScalarTraits<Scalar>::magnitude(vals[k]) > Teuchos::ScalarTraits<Scalar>::magnitude(threshold_) || grow == gcid ) {
      return false; // keep values
    }
    return true;    // values too small -> drop them
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetThreshold() const {
    return threshold_;
  }


  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << "PreDropFunctionConstVal: threshold = " << threshold_ << std::endl;
    return out.str();
  }

  /*template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void PreDropFunctionConstVal<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;
    if (verbLevel & Parameters0) {
      out0 << "PreDropFunctionConstVal: threshold = " << threshold_ << std::endl;
    }
  }*/

}

#define MUELU_PREDROPFUNCTIONCONSTVAL_SHORT
#endif // MUELU_PREDROPFUNCTIONCONSTVAL_DEF_HPP
