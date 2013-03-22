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
#ifndef XPETRA_EXPORTFACTORY_HPP
#define XPETRA_EXPORTFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Export.hpp"

#ifdef HAVE_XPETRA_TPETRA
#include "Xpetra_TpetraExport.hpp"
#endif
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraExport.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ExportFactory {

  private:
    //! Private constructor. This is a static class.
    ExportFactory() {}

  public:

    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      XPETRA_MONITOR("ExportFactory::Build");
      TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

      XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
      XPETRA_FACTORY_END;
    }

  };

  template <>
  class ExportFactory<int, int> {

    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef Kokkos::DefaultNode::DefaultNodeType Node;

  private:
    //! Private constructor. This is a static class.
    ExportFactory() {}

  public:

    static RCP<Export<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source, const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target) {
      XPETRA_MONITOR("ExportFactory::Build");

      TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

#ifdef HAVE_XPETRA_TPETRA
      if (source->lib() == UseTpetra)
        return rcp( new TpetraExport<LocalOrdinal, GlobalOrdinal, Node>(source, target));
#endif

#ifdef HAVE_XPETRA_EPETRA
      if (source->lib() == UseEpetra)
        return rcp( new EpetraExport(source, target));
#endif

      XPETRA_FACTORY_END;
    }

  };

}

#define XPETRA_EXPORTFACTORY_SHORT
#endif
