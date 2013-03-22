// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
/*
  Direct translation of parts of Galeri to use Tpetra or Xpetra rather than Epetra. Epetra also supported.
*/

#ifndef GALERI_XPETRAMATRIXTRAITS_HPP
#define GALERI_XPETRAMATRIXTRAITS_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include "Galeri_config.h"
#ifdef HAVE_GALERI_XPETRA
// needed for the specialized traits:
#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"
#include "Xpetra_Matrix.hpp"
#include "Xpetra_MatrixFactory.hpp"
#endif

namespace Galeri {
  
  namespace Xpetra {
    
    /* Default traits */
    /* These traits work for the following couples of (Map,Matrix):
       - Map = Tpetra::Map<...>,       and Matrix = Tpetra::CrsMatrix<...>  
       - Map = Xpetra::TpetraMap<...> and Matrix = Xpetra::TpetraCrsMatrix<...>
       - Map = Xpetra::EpetraMap,     and Matrix = Xpetra::EpetraCrsMatrix
    */
    template <class Map, class Matrix>
    class MatrixTraits 
    {
    public:
      static Teuchos::RCP<Matrix> Build(const Teuchos::RCP<const Map> &rowMap, size_t maxNumEntriesPerRow) // TODO: pftype
      { return rcp( new Matrix(rowMap, maxNumEntriesPerRow) );
      };
    };

#ifdef HAVE_GALERI_XPETRA

    /* Specialized traits for:
       - Map = Xpetra::Map<...>, Matrix = Xpetra::CrsMatrix<...> */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    class MatrixTraits < ::Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, ::Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps> >
    {
    public:
      static Teuchos::RCP< ::Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow)
      // Use the CrsMatrixFactory to decide what kind of matrix to create (Xpetra::TpetraCrsMatrix or Xpetra::EpetraCrsMatrix).
      { return ::Xpetra::CrsMatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps>::Build(rowMap,  maxNumEntriesPerRow); };
    };

    /* Specialized traits for:
       - Map = Xpetra::Map<...>, Matrix = Xpetra::Matrix<...> */
    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    class MatrixTraits < ::Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node>, ::Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps> >
    {
    public:
      static Teuchos::RCP< ::Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps> > Build(const Teuchos::RCP<const ::Xpetra::Map<LocalOrdinal,GlobalOrdinal, Node> > &rowMap, size_t maxNumEntriesPerRow)
      // Use the CrsMatrixFactory to decide what kind of matrix to create (Xpetra::TpetraCrsMatrix or Xpetra::EpetraCrsMatrix).
      { return ::Xpetra::MatrixFactory<Scalar,LocalOrdinal,GlobalOrdinal, Node, LocalMatOps>::Build(rowMap, maxNumEntriesPerRow); };
    };

#endif

  } // namespace Xpetra
} // namespace Galeri


#endif //ifndef GALERI_XPETRAMATRIXTRAITS_HPP
