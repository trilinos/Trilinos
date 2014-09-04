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
#ifndef MUELU_SCHWARZSMOOTHER_DEF_HPP
#define MUELU_SCHWARZSMOOTHER_DEF_HPP

#include <algorithm>

#include "MueLu_ConfigDefs.hpp"
#if defined (HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2) && defined(HAVE_MUELU_IFPACK2)

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include <Ifpack2_Factory.hpp>

#include "MueLu_SchwarzSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SchwarzSmoother(const std::string& type, const Teuchos::ParameterList& paramList, const LocalOrdinal& overlapLevel)
    : type_(type), paramList_(paramList), overlapLevel_(overlapLevel) {

    this->SetParameterList(paramList);

    if(type_ == "RELAXATION" || type_ == "ILUT"  || type_ == "SCHWARZ" ||
       type_ == "CHEBYSHEV"  || type_ == "RILUK" || type_ == "KRYLOV"  ) {

      this->GetOStream(Warnings0) << "MueLu::SchwarzSmoother: using \"" << type_ << "\"" << std::endl;

    }

    else {

      if (type_ == "" || Amesos2::query(type_) == false) {
        std::string oldtype = type_;
#if defined(HAVE_AMESOS2_SUPERLU)
        type_ = "Superlu";
#elif defined(HAVE_AMESOS2_KLU2)
        type_ = "Klu";
#elif defined(HAVE_AMESOS2_SUPERLUDIST)
        type_ = "Superludist";
#else
        throw Exceptions::RuntimeError("Amesos2 have been compiled without SuperLU_DIST, SuperLU or Klu. By default, MueLu tries"
                                       "to use one of these libraries. Amesos2 must be compiled with one of these solvers or"
                                       "a valid Amesos2 solver have to be specified explicitly.");
#endif
        if (oldtype != "")
          this->GetOStream(Warnings0) << "MueLu::SchwarzSmoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
        else
          this->GetOStream(Warnings0) << "MueLu::SchwarzSmoother: using \"" << type_ << "\"" << std::endl;
      }

      // Check the validity of the solver type parameter
      TEUCHOS_TEST_FOR_EXCEPTION(Amesos2::query(type_) == false, Exceptions::RuntimeError, "The Amesos2 library reported that the solver '" << type_ << "' is not available. "
                                 "Amesos2 have been compiled without the support of this solver or the solver name is misspelled.");

    }

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SchwarzSmoother() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::SchwarzSmoother::Setup(): Setup() has already been called" << std::endl;

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    RCP<Tpetra_CrsMatrix> tA = Utils::Op2NonConstTpetraCrs(A_);

    // make local matrix

    // local communicator
    Teuchos::RCP<const Teuchos::Comm<int> > localComm;
#ifdef HAVE_MPI
    localComm = Teuchos::rcp(new Teuchos::MpiComm<int> (MPI_COMM_SELF));
#else
    localComm = Teuchos::rcp(new Teuchos::SerialComm<int> ());
#endif

    // FIXME (mfh 10 Dec 2013) This looks suspiciously like a
    // duplicate of Ifpack2::AdditiveSchwarz's initialization code.
    // Why not just use Ifpack2::AdditiveSchwarz?

    // first, determine the overlap map
    Array<GO> ExtElements;
    Teuchos::RCP< Tpetra::Map<LO,GO,NO> >       TmpMap;
    Teuchos::RCP< Tpetra::CrsGraph<LO,GO,NO> >  TmpGraph;
    Teuchos::RCP< Tpetra::Import<LO,GO,NO> >    TmpImporter;
    Teuchos::RCP< const Tpetra::Map<LO,GO,NO> > TmpRowMap, TmpColMap;
    Teuchos::RCP< const Tpetra::CrsGraph<LO,GO,NO> > tA_CrsGraph = tA->getCrsGraph();
    const GO global_invalid = Teuchos::OrdinalTraits<GO>::invalid();
    const size_t numMyRowsA = tA->getNodeNumRows();

    // FIXME (mfh 10 Dec 2013) Isn't the 'for' loop below just what
    // Ifpack2::LocalFilter does?  Looks like replicated code to me.
    for (int i = 0; i <= overlapLevel_; ++i) {
      // Get the current maps
      if (i == 0) {
        TmpRowMap = tA->getRowMap ();
        TmpColMap = tA->getColMap ();
      }
      else {
        TmpRowMap = TmpGraph->getRowMap ();
        TmpColMap = TmpGraph->getColMap ();
      }
      const size_t size = TmpColMap->getNodeNumElements () - TmpRowMap->getNodeNumElements ();
      Array<GO> mylist(size);
      size_t count = 0;
      // define the set of rows that are in ColMap but not in RowMap
      for (LO j = 0 ; (size_t) j < TmpColMap->getNodeNumElements(); j++) {
        const GO GID = TmpColMap->getGlobalElement(j);
        if (tA->getRowMap()->getLocalElement(GID) == global_invalid) {
          typedef typename Array<GO>::iterator iter_type;
          const iter_type end = ExtElements.end();
          const iter_type pos = std::find (ExtElements.begin(), end, GID);
          if (pos == end) {
            ExtElements.push_back(GID);
            mylist[count] = GID;
            ++count;
          }
        }
      }
      TmpMap = Teuchos::rcp( new Tpetra::Map<LO,GO,NO>(global_invalid, mylist(0,count),
                                                       Teuchos::OrdinalTraits<GO>::zero(),
                                                       tA->getComm(), tA->getNode()) );
      TmpGraph = Teuchos::rcp( new Tpetra::CrsGraph<LO,GO,NO>(TmpMap, 0));
      TmpImporter = Teuchos::rcp( new Tpetra::Import<LO,GO,NO>(tA->getRowMap(), TmpMap));
      TmpGraph->doImport(*tA_CrsGraph, *TmpImporter, Tpetra::INSERT);
      TmpGraph->fillComplete (tA->getDomainMap(), TmpMap);
    }

    // build the map containing all the nodes (original
    // matrix + extended matrix)
    Array<GO> mylist(numMyRowsA + ExtElements.size());
    for (LO i = 0; (size_t)i < numMyRowsA; ++i) {
      mylist[i] = tA->getRowMap()->getGlobalElement(i);
    }
    for (LO i = 0; i < ExtElements.size(); ++i) {
      mylist[i + numMyRowsA] = ExtElements[i];
    }
    OverlapMap_ = Teuchos::rcp( new Tpetra::Map<LO,GO,NO>(global_invalid, mylist(),
                                                          Teuchos::OrdinalTraits<GO>::zero(),
                                                          tA->getComm(), tA->getNode()) );

    // get maps and setup local matrix
    UniqueMap_  = tA->getRowMap();
    LO numOverlapRows = OverlapMap_->getNodeNumElements();
    LO numUniqueRows = UniqueMap_->getNodeNumElements();
    TEUCHOS_TEST_FOR_EXCEPTION(numUniqueRows>numOverlapRows, Exceptions::RuntimeError, "More unique elements than overlapped elements!");
    Teuchos::ArrayView<const GO> globalRowList = OverlapMap_->getNodeElementList();
    localRowMap_ = Teuchos::rcp(new Tpetra::Map<LO,GO,NO>(numOverlapRows, 0, localComm,
                                                          Tpetra::LocallyReplicated, tA->getNode()));
    Teuchos::RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > OverlapA
      = Teuchos::rcp(new Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>(OverlapMap_,OverlapMap_,100));
    Teuchos::RCP< Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> > LocalA
      = Teuchos::rcp(new Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>(localRowMap_,localRowMap_,100));
    // Import/Export objects
    TpetraExporter_ = Teuchos::rcp (new Tpetra::Export<LO,GO,NO> (OverlapMap_, UniqueMap_));
    TpetraImporter_ = Teuchos::rcp (new Tpetra::Import<LO,GO,NO> (UniqueMap_, OverlapMap_));
    // do import to get overlapped matrix (OverlapA)
    OverlapA->doImport (*tA, *TpetraImporter_, Tpetra::INSERT);
    OverlapA->fillComplete ();

    // extract to local matrix (LocalA)
    for (LO i = 0; i < numOverlapRows; ++i) {
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> values;
      OverlapA->getLocalRowView (i, indices, values);
      LocalA->insertLocalValues (i, indices, values);
    }
    LocalA->fillComplete ();

    if (type_ == "RELAXATION" || type_ == "ILUT"  || type_ == "SCHWARZ" ||
        type_ == "CHEBYSHEV"  || type_ == "RILUK" || type_ == "KRYLOV"  ) {
      Ifpack2::Factory factory;
      ifpack2prec_ = factory.template create<Tpetra_CrsMatrix> (type_, LocalA);
      ifpack2prec_->setParameters (paramList_);
      ifpack2prec_->initialize ();
      ifpack2prec_->compute ();
    }
    else {
      prec_ = Amesos2::create<Tpetra_CrsMatrix, Tpetra_MultiVector> (type_, LocalA);
      TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Amesos2::create returned Teuchos::null");

      prec_->setParameters (Teuchos::rcpFromRef(paramList_));
      prec_->symbolicFactorization ();
      prec_->numericFactorization ();
    }

    SmootherPrototype::IsSetup (true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::SchwarzSmoother::Apply(): Setup() has not been called");

    // Apply
    if (InitialGuessIsZero) {
      Tpetra::MultiVector<SC,LO,GO,NO> &tX = Utils::MV2NonConstTpetraMV(X);
      Tpetra::MultiVector<SC,LO,GO,NO> const &tB = Utils::MV2TpetraMV(B);
      // do import/export of multivector and construct the local vector
      size_t numvecs = X.getNumVectors();
      Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > OverlapB = rcp( new Tpetra::MultiVector<SC,LO,GO,NO>(OverlapMap_,numvecs) );
      OverlapB->doImport(tB,*TpetraImporter_,Tpetra::INSERT);
      OverlapB->replaceMap(localRowMap_);
      Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > OverlapX = Teuchos::rcp( new Tpetra::MultiVector<SC,LO,GO,NO>(localRowMap_,numvecs) );
      if(type_ == "RELAXATION" || type_ == "ILUT"  || type_ == "SCHWARZ" ||
         type_ == "CHEBYSHEV"  || type_ == "RILUK" || type_ == "KRYLOV"  ) {
        ifpack2prec_->apply (*OverlapB, *OverlapX);
      }
      else {
        // solve
        prec_->setX(OverlapX);
        prec_->setB(OverlapB);
        prec_->solve();
      }
      // import to global vector
      OverlapX->replaceMap (OverlapMap_);
      tX.doExport (*OverlapX, *TpetraExporter_, Tpetra::ZERO);
    }

    else {
      RCP<MultiVector> Residual = Utils::Residual(*A_,X,B);
      RCP<MultiVector> Correction = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(A_->getDomainMap(), X.getNumVectors());
      Tpetra::MultiVector<SC,LO,GO,NO> &tX = Utils::MV2NonConstTpetraMV(*Correction);
      Tpetra::MultiVector<SC,LO,GO,NO> const &tB = Utils::MV2TpetraMV(*Residual);
      // do import/export of multivector and construct the local vector
      size_t numvecs = X.getNumVectors();
      Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > OverlapB =
        rcp( new Tpetra::MultiVector<SC,LO,GO,NO>(OverlapMap_,numvecs) );
      OverlapB->doImport(tB,*TpetraImporter_,Tpetra::INSERT);
      OverlapB->replaceMap(localRowMap_);
      Teuchos::RCP< Tpetra::MultiVector<SC,LO,GO,NO> > OverlapX =
        Teuchos::rcp( new Tpetra::MultiVector<SC,LO,GO,NO>(localRowMap_,numvecs) );
      if (type_ == "RELAXATION" || type_ == "ILUT"  || type_ == "SCHWARZ" ||
          type_ == "CHEBYSHEV"  || type_ == "RILUK" || type_ == "KRYLOV"  ) {
        ifpack2prec_->apply (*OverlapB, *OverlapX);
      }
      else {
        // solve
        prec_->setX(OverlapX);
        prec_->setB(OverlapB);
        prec_->solve();
      }
      // extract to global vector
      OverlapX->replaceMap(OverlapMap_);
      tX.doExport(*OverlapX,*TpetraExporter_,Tpetra::ZERO);
      // update
      X.update((Scalar)1.0,*Correction,(Scalar)1.0);
    }

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    return rcp(new SchwarzSmoother(*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    if(type_ == "RELAXATION" || type_ == "ILUT"  || type_ == "SCHWARZ" ||
       type_ == "CHEBYSHEV"  || type_ == "RILUK" || type_ == "KRYLOV"  ) {
      if (SmootherPrototype::IsSetup()) {
        out << ifpack2prec_->description();
      } else {
        out << SmootherPrototype::description();
        out << "{type = " << type_ << "}";
      }
    }
    else {
      if (SmootherPrototype::IsSetup() == true) {
        out << prec_->description();

      } else {
        out << SmootherPrototype::description();
        out << "{type = " << type_ << "}";
      }
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void SchwarzSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
    }

    if ((verbLevel & External) && prec_ != Teuchos::null) {
      Teuchos::OSTab tab2(out);
      out << *prec_ << std::endl;
    }

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
  }

} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_SCHWARZSMOOTHER_DEF_HPP
