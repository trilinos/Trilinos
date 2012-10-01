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
#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK
#include <Ifpack.h>
#include "Xpetra_MultiVectorFactory.hpp"

#include "MueLu_IfpackSmoother.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  IfpackSmoother::IfpackSmoother(std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap, RCP<FactoryBase> AFact) //TODO: empty paramList valid for Ifpack??
    : type_(type), paramList_(paramList), overlap_(overlap), AFact_(AFact)
  { }

  IfpackSmoother::~IfpackSmoother() { }

  void IfpackSmoother::SetParameters(Teuchos::ParameterList const & paramList) {
    paramList_ = paramList;

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...

      Teuchos::ParameterList nonConstParamList = paramList; // because Ifpack SetParameters() input argument is not const...
      prec_->SetParameters(nonConstParamList);
    }
  }

  Teuchos::ParameterList const & IfpackSmoother::GetParameters() { return paramList_; }


  void IfpackSmoother::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
  }

  void IfpackSmoother::Setup(Level &currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true) GetOStream(Warnings0, 0) << "Warning: MueLu::IfpackSmoother::Setup(): Setup() has already been called";

    A_ = currentLevel.Get< RCP<Matrix> >("A", AFact_.get());

    if (type_ == "Chebyshev") {
      Scalar maxEigenValue = paramList_.get("chebyshev: max eigenvalue", (Scalar)-1.0);
      if (maxEigenValue == -1.0) {
        maxEigenValue = Utils::PowerMethod(*A_,true,10,1e-4);
        paramList_.set("chebyshev: max eigenvalue",maxEigenValue);
          
        GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue" << " = " << maxEigenValue << std::endl;
      }
    }

    RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
    Ifpack factory;
    prec_ = rcp(factory.Create(type_, &(*epA), overlap_));
    prec_->SetParameters(paramList_);
    prec_->Compute();

    SmootherPrototype::IsSetup(true);
  }

  void IfpackSmoother::Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack
    Teuchos::ParameterList  paramList;
    if (type_ == "Chebyshev") {
      paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);
    } else if (type_ == "point relaxation stand-alone") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
    } else if  (type_ == "ILU") {
      ;
      /*
      if (InitialGuessIsZero == false) {
        if (IsPrint(Warnings0, 0)) {
          static int warning_only_once=0;
          if ((warning_only_once++) == 0)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::Ifpack2Smoother::Apply(): ILUT has no provision for a nonzero initial guess." << std::endl;
        }
      }
      */
    } else {
      // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
      // we should remove the if/else/elseif and just test if this
      // option is supported by current ifpack2 preconditioner
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,"IfpackSmoother::Apply(): Ifpack preconditioner '"+type_+"' not supported");
    }
    prec_->SetParameters(paramList);

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
  } //Apply

  RCP<MueLu::SmootherPrototype<double, int, int> > IfpackSmoother::Copy() const {
    return rcp(new IfpackSmoother(*this) );
  }

  std::string IfpackSmoother::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }
    
  void IfpackSmoother::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }
      
    if (verbLevel & Parameters1) { 
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      out0 << "Overlap: "        << overlap_ << std::endl;
    }
      
    if (verbLevel & External) {
      if (prec_ != Teuchos::null) { Teuchos::OSTab tab2(out); out << *prec_ << std::endl; }
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<A_>: " << A_ << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
      //TODO: add AFact_
    }
  }

} // namespace MueLu

#endif
