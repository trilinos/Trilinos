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
#ifndef MUELU_IFPACK2SMOOTHER_DEF_HPP
#define MUELU_IFPACK2SMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <Teuchos_ParameterList.hpp>

#include <Tpetra_RowMatrix.hpp>

#include <Ifpack2_Chebyshev.hpp>
#include <Ifpack2_Factory.hpp>
#include <Ifpack2_Parameters.hpp>

#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Ifpack2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Ifpack2Smoother(const std::string& type, const Teuchos::ParameterList& paramList, const LO& overlap)
    : type_(type), overlap_(overlap)
  {
    SetParameterList(paramList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
    Factory::SetParameterList(paramList);

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...
      SetPrecParameters();
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SetPrecParameters(const Teuchos::ParameterList& list) const {
    ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
    paramList.setParameters(list);

    RCP<ParameterList> precList = this->RemoveFactoriesFromList(this->GetParameterList());

    prec_->setParameters(*precList);

    paramList.setParameters(*precList);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (this->IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::Ifpack2Smoother::Setup(): Setup() has already been called";

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    typedef Teuchos::ScalarTraits<SC> STS;
    SC negone = -STS::one();

    SC lambdaMax = negone;

    // If we are doing "user" partitioning, we assume that what the user
    // really wants to do is make tiny little subdomains with one row
    // asssigned to each subdomain. The rows used for these little
    // subdomains correspond to those in the 2nd block row.  Then,
    // if we overlap these mini-subdomains, we will do something that
    // looks like Vanka (grabbing all velocities associated with each
    // each pressure unknown). In addition, we put all Dirichlet points
    // as a little mini-domain.

    bool isBlockedMatrix = false;
    RCP<Matrix> merged2Mat;
    if (type_ == "SCHWARZ") {
      ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

      std::string sublistName = "subdomain solver parameters";
      if (paramList.isSublist(sublistName)) {
        ParameterList& subList = paramList.sublist(sublistName);

        std::string partName  = "partitioner: type";
        if (subList.isParameter(partName) && subList.get<std::string>(partName) == "user") {
          isBlockedMatrix = true;

          RCP<BlockedCrsMatrix> bA = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
          TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                                     "Matrix A must be of type BlockedCrsMatrix.");

          size_t numVels = bA->getMatrix(0,0)->getNodeNumRows();
          size_t numPres = bA->getMatrix(1,0)->getNodeNumRows();
          size_t numRows = A_->getNodeNumRows();

          ArrayRCP<LocalOrdinal> blockSeeds(numRows, Teuchos::OrdinalTraits<LocalOrdinal>::invalid());

          for (size_t rowOfB = numVels; rowOfB < numVels+numPres; ++rowOfB)
            blockSeeds[rowOfB] = rowOfB - numVels;

          RCP<BlockedCrsMatrix> bA2 = rcp_dynamic_cast<BlockedCrsMatrix>(A_);
          TEUCHOS_TEST_FOR_EXCEPTION(bA2.is_null(), Exceptions::BadCast,
                                     "Matrix A must be of type BlockedCrsMatrix.");

          RCP<CrsMatrix> mergedMat = bA2->Merge();
          merged2Mat = rcp(new CrsMatrixWrap(mergedMat));

          // Add Dirichlet rows to the list of seeds
          ArrayRCP<const bool> boundaryNodes;
          boundaryNodes = Utils::DetectDirichletRows(*merged2Mat, 0.0);
          for (LO i = 0; i < boundaryNodes.size(); i++)
            if (boundaryNodes[i]) {
              blockSeeds[i] = numPres;
              numPres++;
            }

          subList.set("partitioner: map",         blockSeeds);
          subList.set("partitioner: local parts", as<int>(numPres));
        }
      }
    } // if (type_ == "SCHWARZ")

    if (type_ == "CHEBYSHEV") {
      std::string maxEigString   = "chebyshev: max eigenvalue";
      std::string eigRatioString = "chebyshev: ratio eigenvalue";

      ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());

      // Get/calculate the maximum eigenvalue
      if (paramList.isParameter(maxEigString)) {
        if (paramList.isType<double>(maxEigString))
          lambdaMax = paramList.get<double>(maxEigString);
        else
          lambdaMax = paramList.get<SC>(maxEigString);
        this->GetOStream(Statistics1) << maxEigString << " (cached with smoother parameter list) = " << lambdaMax << std::endl;

      } else {
        lambdaMax = A_->GetMaxEigenvalueEstimate();
        if (lambdaMax != negone) {
          this->GetOStream(Statistics1) << maxEigString << " (cached with matrix) = " << lambdaMax << std::endl;
          paramList.set(maxEigString, lambdaMax);
        }
      }

      // Calculate the eigenvalue ratio
      const SC defaultEigRatio = 20;

      SC ratio = defaultEigRatio;
      if (paramList.isParameter(eigRatioString)) {
        if (paramList.isType<double>(eigRatioString))
          ratio = paramList.get<double>(eigRatioString);
        else
          ratio = paramList.get<SC>(eigRatioString);
      }
      if (currentLevel.GetLevelID()) {
        // Update ratio to be
        //   ratio = max(number of fine DOFs / number of coarse DOFs, defaultValue)
        //
        // NOTE: We don't need to request previous level matrix as we know for sure it was constructed
        RCP<const Matrix> fineA = currentLevel.GetPreviousLevel()->Get<RCP<Matrix> >("A");
        size_t nRowsFine   = fineA->getGlobalNumRows();
        size_t nRowsCoarse = A_->getGlobalNumRows();

        SC levelRatio = as<SC>(as<float>(nRowsFine)/nRowsCoarse);
        if (STS::magnitude(levelRatio) > STS::magnitude(ratio))
          ratio = levelRatio;
      }

      this->GetOStream(Statistics1) << eigRatioString << " (computed) = " << ratio << std::endl;
      paramList.set(eigRatioString, ratio);
    }

    RCP<const Tpetra::RowMatrix<SC, LO, GO, NO> > tpA;
    if (isBlockedMatrix == true) tpA = Utils::Op2NonConstTpetraRow(merged2Mat);
    else                         tpA = Utils::Op2NonConstTpetraRow(A_);
    prec_ = Ifpack2::Factory::create(type_, tpA, overlap_);
    SetPrecParameters();
    prec_->initialize();
    prec_->compute();

    SmootherPrototype::IsSetup(true);

    if (type_ == "CHEBYSHEV" && lambdaMax == negone) {
      typedef Tpetra::RowMatrix<SC, LO, GO, NO> MatrixType;

      Teuchos::RCP<Ifpack2::Chebyshev<MatrixType> > chebyPrec = rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType> >(prec_);
      if (chebyPrec != Teuchos::null) {
        lambdaMax = chebyPrec->getLambdaMaxForApply();
        A_->SetMaxEigenvalueEstimate(lambdaMax);
        this->GetOStream(Statistics1) << "chebyshev: max eigenvalue (calculated by Ifpack2)" << " = " << lambdaMax << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(lambdaMax == negone, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Setup(): no maximum eigenvalue estimate");
    }

    this->GetOStream(Statistics0) << description() << std::endl;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack2
    // TODO:  It might be nice to switch back the internal
    //        "zero starting solution" option of the ifpack2 object prec_ to his
    //        initial value at the end but there is no way right now to get
    //        the current value of the "zero starting solution" in ifpack2.
    //        It's not really an issue, as prec_  can only be used by this method.
    // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
    // we should remove the if/else/elseif and just test if this
    // option is supported by current ifpack2 preconditioner
    Teuchos::ParameterList paramList;
    bool supportInitialGuess = false;
    if (type_ == "CHEBYSHEV") {
      paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);
      SetPrecParameters(paramList);
      supportInitialGuess = true;

    } else if (type_ == "RELAXATION") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
      SetPrecParameters(paramList);
      supportInitialGuess = true;

    } else if (type_ == "KRYLOV") {
      paramList.set("krylov: zero starting solution", InitialGuessIsZero);
      SetPrecParameters(paramList);
      supportInitialGuess = true;

    }
    //TODO JJH 30Apr2014  Calling SetPrecParameters(paramList) when the smoother
    //is Ifpack2::AdditiveSchwarz::setParameterList() will destroy the subdomain
    //(aka inner) solver.  This behavior is documented but a departure from what
    //it previously did, and what other Ifpack2 solvers currently do.  So I have
    //moved SetPrecParameters(paramList) into the if-else block above.

    // Apply
    if (InitialGuessIsZero || supportInitialGuess) {
      Tpetra::MultiVector<SC,LO,GO,NO>&       tpX = Utils::MV2NonConstTpetraMV(X);
      const Tpetra::MultiVector<SC,LO,GO,NO>& tpB = Utils::MV2TpetraMV(B);

      prec_->apply(tpB, tpX);

    } else {
      typedef Teuchos::ScalarTraits<Scalar> TST;
      RCP<MultiVector> Residual   = Utils::Residual(*A_, X, B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

      Tpetra::MultiVector<SC,LO,GO,NO>&       tpX = Utils::MV2NonConstTpetraMV(*Correction);
      const Tpetra::MultiVector<SC,LO,GO,NO>& tpB = Utils::MV2TpetraMV(*Residual);

      prec_->apply(tpB, tpX);

      X.update(TST::one(), *Correction, TST::one());
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    RCP<Ifpack2Smoother> smoother = rcp(new Ifpack2Smoother(*this) );
    smoother->SetParameterList(this->GetParameterList());
    return smoother;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    if (SmootherPrototype::IsSetup()) {
      out << prec_->description();
    } else {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl;
      Teuchos::OSTab tab2(out);
      out << this->GetParameterList();
      out0 << "Overlap: "        << overlap_ << std::endl;
    }

    if (verbLevel & External)
      if (prec_ != Teuchos::null) {
        Teuchos::OSTab tab2(out);
        out << *prec_ << std::endl;
      }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_IFPACK2
#endif // MUELU_IFPACK2SMOOTHER_DEF_HPP
