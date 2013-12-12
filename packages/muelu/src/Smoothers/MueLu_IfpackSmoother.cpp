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

    paramList.setParameters(*precList);
  }

  void IfpackSmoother::DeclareInput(Level &currentLevel) const {
    this->Input(currentLevel, "A");
  }

  void IfpackSmoother::Setup(Level &currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true) GetOStream(Warnings0, 0) << "Warning: MueLu::IfpackSmoother::Setup(): Setup() has already been called";

    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    double lambdaMax = -1.0;
    if (type_ == "Chebyshev")
      try {
        lambdaMax = Teuchos::getValue<Scalar>(this->GetParameter("chebyshev: max eigenvalue"));
        this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue (cached with smoother parameter list)" << " = " << lambdaMax << std::endl;

      } catch (Teuchos::Exceptions::InvalidParameterName) {
        lambdaMax = A_->GetMaxEigenvalueEstimate();
        if (lambdaMax != -1.0) {
          this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue (cached with matrix)" << " = " << lambdaMax << std::endl;
          this->SetParameter("chebyshev: max eigenvalue", ParameterEntry(lambdaMax));
        }
      }

    RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);

    Ifpack factory;
    prec_ = rcp(factory.Create(type_, &(*epA), overlap_));
    SetPrecParameters();
    prec_->Compute();

    SmootherPrototype::IsSetup(true);

    if (type_ == "Chebyshev" && lambdaMax == -1.0) {
      Teuchos::RCP<Ifpack_Chebyshev> chebyPrec = rcp_dynamic_cast<Ifpack_Chebyshev>(prec_);
      if (chebyPrec != Teuchos::null) {
        lambdaMax = chebyPrec->GetLambdaMax();
        A_->SetMaxEigenvalueEstimate(lambdaMax);
        this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue (calculated by Ifpack)" << " = " << lambdaMax << std::endl;
      }
      TEUCHOS_TEST_FOR_EXCEPTION(lambdaMax == -1.0, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Setup(): no maximum eigenvalue estimate");
    }

    this->GetOStream(Statistics0, 0) << description() << std::endl;
  }

  void IfpackSmoother::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack
    Teuchos::ParameterList  paramList;
    if (type_ == "Chebyshev") {
      paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);

    } else if (type_ == "point relaxation stand-alone") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);

    } else if  (type_ == "ILU") {
      // do nothing

    } else {
      // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
      // we should remove the if/else/elseif and just test if this
      // option is supported by current ifpack2 preconditioner
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,"IfpackSmoother::Apply(): Ifpack preconditioner '"+type_+"' not supported");
    }
    SetPrecParameters(paramList);

    // Apply
    if (InitialGuessIsZero) {
      Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
      Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
      prec_->ApplyInverse(epB, epX);
    } else {
      RCP<MultiVector> Residual = Utils::Residual(*A_,X,B);
      RCP<MultiVector> Correction = MultiVectorFactory::Build(A_->getDomainMap(), X.getNumVectors());
      Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(*Correction);
      Epetra_MultiVector const &epB = Utils::MV2EpetraMV(*Residual);
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
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
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
