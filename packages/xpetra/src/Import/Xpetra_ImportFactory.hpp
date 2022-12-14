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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_IMPORTFACTORY_HPP
#define XPETRA_IMPORTFACTORY_HPP

#include "Xpetra_ConfigDefs.hpp"

#include "Xpetra_Import.hpp"

#include "Xpetra_TpetraImport.hpp"
#ifdef HAVE_XPETRA_EPETRA
#include "Xpetra_EpetraImport.hpp"
#endif

#include "Xpetra_Exceptions.hpp"

namespace Xpetra {

  template <class LocalOrdinal,
            class GlobalOrdinal,
            class Node = KokkosClassic::DefaultNode::DefaultNodeType>
  class ImportFactory {
  private:
    //! Private constructor. This is a static class.
    ImportFactory() {}

  public:

    //! Constructor specifying the number of non-zeros for all rows.
    static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source,
                                                                 const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target,
                                                                 const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
      XPETRA_MONITOR("ImportFactory::Build");

      TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target, plist));

      XPETRA_FACTORY_ERROR_IF_EPETRA(source->lib());
      XPETRA_FACTORY_END;
    }

  };

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES))

  // Specialization on LO=GO=int with serial node.
  // Used for Epetra and Tpetra
  // For any other node definition the general default implementation is used which allows Tpetra only
  template <>
  class ImportFactory<int, int, EpetraNode> {
    typedef int LocalOrdinal;
    typedef int GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    ImportFactory() {}

  public:

    static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source,
                                                                 const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target,
                                                                 const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
      XPETRA_MONITOR("ImportFactory::Build");
      TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target,plist));

      if (source->lib() == UseEpetra)
        return rcp( new EpetraImportT<int,Node>(source, target));

      XPETRA_FACTORY_END;
    }

  };
#endif

// we need the Epetra specialization only if Epetra is enabled
#if (defined(HAVE_XPETRA_EPETRA) && !defined(XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES))
  template <>
  class ImportFactory<int, long long, EpetraNode> {

    typedef int LocalOrdinal;
    typedef long long GlobalOrdinal;
    typedef EpetraNode Node;

  private:
    //! Private constructor. This is a static class.
    ImportFactory() {}

  public:

    static RCP<Import<LocalOrdinal, GlobalOrdinal, Node> > Build(const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &source,
                                                                 const RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> > &target,
                                                                 const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null) {
      XPETRA_MONITOR("ImportFactory::Build");
      TEUCHOS_TEST_FOR_EXCEPTION(source->lib() != target->lib(), Xpetra::Exceptions::RuntimeError, "");

      if (source->lib() == UseTpetra)
        return rcp( new TpetraImport<LocalOrdinal, GlobalOrdinal, Node>(source, target, plist));

      if (source->lib() == UseEpetra)
        return rcp( new EpetraImportT<long long,Node>(source, target));

      XPETRA_FACTORY_END;
    }

  };
#endif
}

#define XPETRA_IMPORTFACTORY_SHORT
#endif
