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

  IfpackSmoother::IfpackSmoother(std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap)
    : type_(type), overlap_(overlap)
  {
    SetParameterList(paramList);
  }

  void IfpackSmoother::SetParameterList(const Teuchos::ParameterList& paramList) {
    Factory::SetParameterList(paramList);

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...
      prec_->SetParameters(const_cast<ParameterList&>(this->GetParameterList()));
    }
  }

  void IfpackSmoother::SetPrecParameters(const Teuchos::ParameterList& list) const {
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

  void IfpackSmoother::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");
  }

  void IfpackSmoother::Setup(Level &currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true)
      GetOStream(Warnings0) << "MueLu::IfpackSmoother::Setup(): Setup() has already been called";

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    double lambdaMax = -1.0;
    if (type_ == "Chebyshev") {
      std::string maxEigString   = "chebyshev: max eigenvalue";
      std::string eigRatioString = "chebyshev: ratio eigenvalue";

      try {
        lambdaMax = Teuchos::getValue<Scalar>(this->GetParameter(maxEigString));
        this->GetOStream(Statistics1) << maxEigString << " (cached with smoother parameter list) = " << lambdaMax << std::endl;

      } catch (Teuchos::Exceptions::InvalidParameterName) {
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

      } catch (Teuchos::Exceptions::InvalidParameterName) {
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
    }

    RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);

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

    this->GetOStream(Statistics0) << description() << std::endl;
  }

  void IfpackSmoother::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
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
      Epetra_MultiVector&       epX = Utils::MV2NonConstEpetraMV(X);
      const Epetra_MultiVector& epB = Utils::MV2EpetraMV(B);

      prec_->ApplyInverse(epB, epX);

    } else {
      RCP<MultiVector> Residual   = Utils::Residual(*A_, X, B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());

      Epetra_MultiVector&       epX = Utils::MV2NonConstEpetraMV(*Correction);
      const Epetra_MultiVector& epB = Utils::MV2EpetraMV(*Residual);

      prec_->ApplyInverse(epB, epX);

      X.update(1.0, *Correction, 1.0);
    }
  }

  RCP<MueLu::SmootherPrototype<double, int, int> > IfpackSmoother::Copy() const {
    RCP<IfpackSmoother> smoother = rcp(new IfpackSmoother(*this) );
    smoother->SetParameterList(this->GetParameterList());
    return smoother;
  }

  std::string IfpackSmoother::description() const {
    std::ostringstream out;
    // The check "GetVerbLevel() == Test" is to avoid
    // failures in the EasyInterface test.
    if (prec_ == Teuchos::null || GetVerbLevel() == Test) {
      out << SmootherPrototype::description();
      out << "{type = " << type_ << "}";
    } else {
      out << prec_->Label();
    }
    return out.str();
  }

  void IfpackSmoother::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
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

} // namespace MueLu

#endif
