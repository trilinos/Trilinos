// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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

// WARNING: This code is experimental. Backwards compatibility should not be expected.

#ifndef XPETRA_OPERATORFACTORY_HPP
#define XPETRA_OPERATORFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Operator.hpp"
#include "Xpetra_CrsOperator.hpp"
#include "Xpetra_Map.hpp"
#include "Xpetra_Vector.hpp"
#include "Xpetra_Exceptions.hpp"

namespace Xpetra {
  
  template <class Scalar, class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps>
  class OperatorFactory {
#undef XPETRA_OPERATORFACTORY_SHORT
#include "Xpetra_UseShortNames.hpp"

  private:
    //! Private constructor. This is a static class.
    OperatorFactory() {}
    
  public:
    
    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Operator> Build(const RCP<const Map> &rowMap, size_t maxNumEntriesPerRow, Xpetra::ProfileType pftype = Xpetra::DynamicProfile) {
      // if const block size && blocksize == 1

      return rcp( new CrsOperator(rowMap, maxNumEntriesPerRow, pftype) );

      // elseif
      
      // return vbr

      // else

      // TEUCHOS_TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"?");
    }

    //! Constructor specifying (possibly different) number of entries in each row.
    static RCP<Operator> Build(const RCP<const Map> &rowMap, const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc, ProfileType pftype = Xpetra::DynamicProfile) {
      return rcp( new CrsOperator(rowMap, NumEntriesPerRowToAlloc, pftype) );
    }

    //! Constructor for creating a diagonal Xpetra::Operator using the entries of a given vector for the diagonal
    static RCP<Operator> Build(const RCP<const Vector> & diagonal) {
      Teuchos::ArrayRCP<const Scalar> vals = diagonal->getData(0);
      Teuchos::RCP<CrsOperator> mtx = Teuchos::rcp(new CrsOperator(diagonal->getMap(), 1, Xpetra::StaticProfile));
      LocalOrdinal NumMyElements = diagonal->getMap()->getNodeNumElements();
      Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = diagonal->getMap()->getNodeElementList();
      for (LocalOrdinal i = 0; i < NumMyElements; ++i) {
          mtx->insertGlobalValues(MyGlobalElements[i],
                                  Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                                  Teuchos::tuple<Scalar>(vals[i]) );
      }
      mtx->fillComplete();
      return mtx;
    }

  };

}

#define XPETRA_OPERATORFACTORY_SHORT
#endif
