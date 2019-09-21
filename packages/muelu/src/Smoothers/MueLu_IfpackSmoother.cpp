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
#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_IFPACK)
#include <Ifpack.h>
#include <Ifpack_Chebyshev.h>
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Node>
  IfpackSmoother<Node>::IfpackSmoother(std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap)
    : type_(type), overlap_(overlap)
  {
    SetParameterList(paramList);
  }

  template <class Node>
  void IfpackSmoother<Node>::SetParameterList(const Teuchos::ParameterList& paramList) {
    Factory::SetParameterList(paramList);

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...
      prec_->SetParameters(const_cast<ParameterList&>(this->GetParameterList()));
    }
  }

  template <class Node>
  void IfpackSmoother<Node>::SetPrecParameters(const Teuchos::ParameterList& list) const {
    ParameterList& paramList = const_cast<ParameterList&>(this->GetParameterList());
    paramList.setParameters(list);

    RCP<ParameterList> precList = this->RemoveFactoriesFromList(this->GetParameterList());

    prec_->SetParameters(*precList);

    // We would like to have the following line here:
    //      paramList.setParameters(*precList);
    // For instance, if Ifpack sets somem parameters internally, we would like to have
    // them listed when we call this->GetParameterList()
    // But because of the way Ifpack handles the list, we cannot do that.
    // The bad scenario goes like this:
    //   * SmootherFactory calls Setup
    //   * Setup calls SetPrecParameters
    //   * We call prec_->SetParameters(*precList)
    //     This actually updates the internal parameter list  with default prec_ parameters
    //     This means that we get a parameter ("chebyshev: max eigenvalue", -1) in the list
    //   * Setup calls prec_->Compute()
    //     Here we may compute the max eigenvalue, but we get no indication of this. If we
    //     do compute it, our parameter list becomes outdated
    //   * SmootherFactory calls Apply
    //   * Apply constructs a list with a list with an entry "chebyshev: zero starting solution"
    //   * We call prec_->SetParameters(*precList)
    // The last call is the problem. At this point, we have a list with an outdated entry
    // "chebyshev: max eigenvalue", but prec_ uses this entry and replaces the computed max
    // eigenvalue with the one from the list, resulting in -1.0 eigenvalue.
    //
    // Ifpack2 does not have this problem, as it does not populate the list with new entries
  }

  template <class Node>
  void IfpackSmoother<Node>::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");

    if (type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
        type_ == "LINESMOOTHING_BANDED RELAXATION" ||
        type_ == "LINESMOOTHING_BANDEDRELAXATION"  ||
        type_ == "LINESMOOTHING_BLOCK_RELAXATION"  ||
        type_ == "LINESMOOTHING_BLOCK RELAXATION"  ||
        type_ == "LINESMOOTHING_BLOCKRELAXATION") {
      this->Input(currentLevel, "CoarseNumZLayers");              // necessary for fallback criterion
      this->Input(currentLevel, "LineDetection_VertLineIds"); // necessary to feed block smoother
    } // if (type_ == "LINESMOOTHING_BANDEDRELAXATION")
  }

  template <class Node>
  void IfpackSmoother<Node>::Setup(Level &currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::IfpackSmoother::Setup(): Setup() has already been called" << std::endl;

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    double lambdaMax = -1.0;
    if (type_ == "Chebyshev") {
      std::string maxEigString   = "chebyshev: max eigenvalue";
      std::string eigRatioString = "chebyshev: ratio eigenvalue";

      try {
        lambdaMax = Teuchos::getValue<Scalar>(this->GetParameter(maxEigString));
        this->GetOStream(Statistics1) << maxEigString << " (cached with smoother parameter list) = " << lambdaMax << std::endl;

      } catch (Teuchos::Exceptions::InvalidParameterName&) {
        lambdaMax = A_->GetMaxEigenvalueEstimate();

        if (lambdaMax != -1.0) {
          this->GetOStream(Statistics1) << maxEigString << " (cached with matrix) = " << lambdaMax << std::endl;
          this->SetParameter(maxEigString, ParameterEntry(lambdaMax));
        }
      }

      // Calculate the eigenvalue ratio
      const Scalar defaultEigRatio = 20;

      Scalar ratio = defaultEigRatio;
      try {
        ratio = Teuchos::getValue<Scalar>(this->GetParameter(eigRatioString));

      } catch (Teuchos::Exceptions::InvalidParameterName&) {
        this->SetParameter(eigRatioString, ParameterEntry(ratio));
      }

      if (currentLevel.GetLevelID()) {
        // Update ratio to be
        //   ratio = max(number of fine DOFs / number of coarse DOFs, defaultValue)
        //
        // NOTE: We don't need to request previous level matrix as we know for sure it was constructed
        RCP<const Matrix> fineA = currentLevel.GetPreviousLevel()->Get<RCP<Matrix> >("A");
        size_t nRowsFine   = fineA->getGlobalNumRows();
        size_t nRowsCoarse = A_->getGlobalNumRows();

        ratio = std::max(ratio, as<Scalar>(nRowsFine)/nRowsCoarse);

        this->GetOStream(Statistics1) << eigRatioString << " (computed) = " << ratio << std::endl;
        this->SetParameter(eigRatioString, ParameterEntry(ratio));
      }
    } // if (type_ == "Chebyshev")

    if (type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
        type_ == "LINESMOOTHING_BANDED RELAXATION" ||
        type_ == "LINESMOOTHING_BANDEDRELAXATION"  ||
        type_ == "LINESMOOTHING_BLOCK_RELAXATION"  ||
        type_ == "LINESMOOTHING_BLOCK RELAXATION"  ||
        type_ == "LINESMOOTHING_BLOCKRELAXATION" ) {
      ParameterList& myparamList = const_cast<ParameterList&>(this->GetParameterList());

      LO CoarseNumZLayers = currentLevel.Get<LO>("CoarseNumZLayers",Factory::GetFactory("CoarseNumZLayers").get());
      if (CoarseNumZLayers > 0) {
        Teuchos::ArrayRCP<LO> TVertLineIdSmoo = currentLevel.Get< Teuchos::ArrayRCP<LO> >("LineDetection_VertLineIds", Factory::GetFactory("LineDetection_VertLineIds").get());

        // determine number of local parts
        LO maxPart = 0;
        for(size_t k = 0; k < Teuchos::as<size_t>(TVertLineIdSmoo.size()); k++) {
          if(maxPart < TVertLineIdSmoo[k]) maxPart = TVertLineIdSmoo[k];
        }

        size_t numLocalRows = A_->getNodeNumRows();
        TEUCHOS_TEST_FOR_EXCEPTION(numLocalRows % TVertLineIdSmoo.size() != 0, Exceptions::RuntimeError, "MueLu::Ifpack2Smoother::Setup(): the number of local nodes is incompatible with the TVertLineIdsSmoo.");

        if (numLocalRows == Teuchos::as<size_t>(TVertLineIdSmoo.size())) {
          myparamList.set("partitioner: type","user");
          myparamList.set("partitioner: map",&(TVertLineIdSmoo[0]));
          myparamList.set("partitioner: local parts",maxPart+1);
        } else {
          // we assume a constant number of DOFs per node
          size_t numDofsPerNode = numLocalRows / TVertLineIdSmoo.size();

          // Create a new Teuchos::ArrayRCP<LO> of size numLocalRows and fill it with the corresponding information
          Teuchos::ArrayRCP<LO> partitionerMap(numLocalRows, Teuchos::OrdinalTraits<LocalOrdinal>::invalid());
          for (size_t blockRow = 0; blockRow < Teuchos::as<size_t>(TVertLineIdSmoo.size()); ++blockRow)
            for (size_t dof = 0; dof < numDofsPerNode; dof++)
              partitionerMap[blockRow * numDofsPerNode + dof] = TVertLineIdSmoo[blockRow];
          myparamList.set("partitioner: type","user");
          myparamList.set("partitioner: map",&(partitionerMap[0]));
          myparamList.set("partitioner: local parts",maxPart + 1);
        }

        if (type_ == "LINESMOOTHING_BANDED_RELAXATION" ||
            type_ == "LINESMOOTHING_BANDED RELAXATION" ||
            type_ == "LINESMOOTHING_BANDEDRELAXATION")
          type_ = "block relaxation";
        else
          type_ = "block relaxation";
      } else {
        // line detection failed -> fallback to point-wise relaxation
        this->GetOStream(Runtime0) << "Line detection failed: fall back to point-wise relaxation" << std::endl;
        myparamList.remove("partitioner: type",false);
        myparamList.remove("partitioner: map", false);
        myparamList.remove("partitioner: local parts",false);
        type_ = "point relaxation stand-alone";
      }

    } // if (type_ == "LINESMOOTHING_BANDEDRELAXATION")

    RCP<Epetra_CrsMatrix> epA = Utilities::Op2NonConstEpetraCrs(A_);

    Ifpack factory;
    prec_ = rcp(factory.Create(type_, &(*epA), overlap_));
    TEUCHOS_TEST_FOR_EXCEPTION(prec_.is_null(), Exceptions::RuntimeError, "Could not create an Ifpack preconditioner with type = \"" << type_ << "\"");
    SetPrecParameters();
    prec_->Compute();

    SmootherPrototype::IsSetup(true);

    if (type_ == "Chebyshev" && lambdaMax == -1.0) {
      Teuchos::RCP<Ifpack_Chebyshev> chebyPrec = rcp_dynamic_cast<Ifpack_Chebyshev>(prec_);
      if (chebyPrec != Teuchos::null) {
        lambdaMax = chebyPrec->GetLambdaMax();
        A_->SetMaxEigenvalueEstimate(lambdaMax);
        this->GetOStream(Statistics1) << "chebyshev: max eigenvalue (calculated by Ifpack)" << " = " << lambdaMax << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(lambdaMax == -1.0, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Setup(): no maximum eigenvalue estimate");
    }

    this->GetOStream(Statistics1) << description() << std::endl;
  }

  template <class Node>
  void IfpackSmoother<Node>::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack
    Teuchos::ParameterList paramList;
    bool supportInitialGuess = false;
    if (type_ == "Chebyshev") {
      paramList.set("chebyshev: zero starting solution",  InitialGuessIsZero);
      supportInitialGuess = true;

    } else if (type_ == "point relaxation stand-alone") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
      supportInitialGuess = true;
    }

    SetPrecParameters(paramList);

    // Apply
    if (InitialGuessIsZero || supportInitialGuess) {
      Epetra_MultiVector&       epX = Utilities::MV2NonConstEpetraMV(X);
      const Epetra_MultiVector& epB = Utilities::MV2EpetraMV(B);

      prec_->ApplyInverse(epB, epX);

    } else {
      RCP<MultiVector> Residual   = Utilities::Residual(*A_, X, B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

      Epetra_MultiVector&       epX = Utilities::MV2NonConstEpetraMV(*Correction);
      const Epetra_MultiVector& epB = Utilities::MV2EpetraMV(*Residual);

      prec_->ApplyInverse(epB, epX);

      X.update(1.0, *Correction, 1.0);
    }
  }

  template <class Node>
  RCP<MueLu::SmootherPrototype<double, int, int, Node> > IfpackSmoother<Node>::Copy() const {
    RCP<IfpackSmoother<Node> > smoother = rcp(new IfpackSmoother<Node>(*this) );
    smoother->SetParameterList(this->GetParameterList());
    return Teuchos::rcp_dynamic_cast<MueLu::SmootherPrototype<double, int, int, Node> >(smoother);
  }

  template <class Node>
  std::string IfpackSmoother<Node>::description() const {
    std::ostringstream out;
    // The check "GetVerbLevel() == Test" is to avoid
    // failures in the EasyInterface test.
    if (prec_ == Teuchos::null || this->GetVerbLevel() == Test) {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    } else {
      out << prec_->Label();
    }
    return out.str();
  }

  template <class Node>
  void IfpackSmoother<Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
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
           << "RCP<A_>: " << A_ << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

  template <class Node>
  size_t IfpackSmoother<Node>::getNodeSmootherComplexity() const {
    // FIXME: This is a placeholder
    return Teuchos::OrdinalTraits<size_t>::invalid();
  }
  

} // namespace MueLu

// The IfpackSmoother is only templated on the Node, since it is an Epetra only object
// Therefore we do not need the full ETI instantiations as we do for the other MueLu
// objects which are instantiated on all template parameters.
#if defined(HAVE_MUELU_EPETRA)
template class MueLu::IfpackSmoother<Xpetra::EpetraNode>;
#endif

#endif
