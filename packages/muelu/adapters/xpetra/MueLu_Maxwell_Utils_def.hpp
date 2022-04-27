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
#ifndef MUELU_MAXWELL_UTILS_DEF_HPP
#define MUELU_MAXWELL_UTILS_DEF_HPP

#include <sstream>

#include "MueLu_ConfigDefs.hpp"

#include "Xpetra_Map.hpp"
#include "Xpetra_CrsMatrixUtils.hpp"
#include "Xpetra_MatrixUtils.hpp"

#include "MueLu_Maxwell_Utils_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_ThresholdAFilterFactory.hpp"
#include "MueLu_Utilities.hpp"

#ifdef HAVE_MUELU_KOKKOS_REFACTOR
#include "MueLu_Utilities_kokkos.hpp"
#endif

// Stratimikos
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
// Thyra includes
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_MueLuPreconditionerFactory.hpp>
#include "Teuchos_AbstractFactoryStd.hpp"
// Ifpack2 includes
#ifdef HAVE_MUELU_IFPACK2
#include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif
#endif

namespace MueLu {


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell_Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::detectBoundaryConditionsSM(RCP<Matrix> & SM_Matrix_,
                                                                                         RCP<Matrix> & D0_Matrix_,
                                                                                         magnitudeType rowSumTol,
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
                                                                                         bool useKokkos_,
                                                                                         Kokkos::View<bool*, typename Node::device_type> & BCrowsKokkos_, 
                                                                                         Kokkos::View<bool*, typename Node::device_type> & BCcolsKokkos_,   
                                                                                         Kokkos::View<bool*, typename Node::device_type> & BCdomainKokkos_,
#endif
                                                                                         int & BCedges_, 
                                                                                         int & BCnodes_,
                                                                                         Teuchos::ArrayRCP<bool> & BCrows_,
                                                                                         Teuchos::ArrayRCP<bool> & BCcols_,
                                                                                         Teuchos::ArrayRCP<bool> & BCdomain_,
                                                                                         bool & allEdgesBoundary_,
                                                                                         bool & allNodesBoundary_) {
    // clean rows associated with boundary conditions
    // Find rows with only 1 or 2 nonzero entries, record them in BCrows_.
    // BCrows_[i] is true, iff i is a boundary row
    // BCcols_[i] is true, iff i is a boundary column
    int BCedgesLocal = 0;
    int BCnodesLocal = 0;
#ifdef HAVE_MUELU_KOKKOS_REFACTOR
    if (useKokkos_) {
      BCrowsKokkos_ = Utilities_kokkos::DetectDirichletRows(*SM_Matrix_,Teuchos::ScalarTraits<magnitudeType>::eps(),/*count_twos_as_dirichlet=*/true);

      if (rowSumTol > 0.)
        Utilities_kokkos::ApplyRowSumCriterion(*SM_Matrix_, rowSumTol, BCrowsKokkos_);

      BCcolsKokkos_ = Kokkos::View<bool*,typename Node::device_type>(Kokkos::ViewAllocateWithoutInitializing("dirichletCols"), D0_Matrix_->getColMap()->getLocalNumElements());
      BCdomainKokkos_ = Kokkos::View<bool*,typename Node::device_type>(Kokkos::ViewAllocateWithoutInitializing("dirichletDomains"), D0_Matrix_->getDomainMap()->getLocalNumElements());
      Utilities_kokkos::DetectDirichletColsAndDomains(*D0_Matrix_,BCrowsKokkos_,BCcolsKokkos_,BCdomainKokkos_);

      auto BCrowsKokkos=BCrowsKokkos_;
      Kokkos::parallel_reduce(BCrowsKokkos_.size(), KOKKOS_LAMBDA (int i, int & sum) {
        if (BCrowsKokkos(i))
	  ++sum;
	}, BCedgesLocal );
         
      auto BCdomainKokkos = BCdomainKokkos_;
      Kokkos::parallel_reduce(BCdomainKokkos_.size(), KOKKOS_LAMBDA (int i, int & sum) {
        if (BCdomainKokkos(i))
	  ++sum;
	}, BCnodesLocal);
    } else
#endif // HAVE_MUELU_KOKKOS_REFACTOR
    {
      BCrows_ = Teuchos::arcp_const_cast<bool>(Utilities::DetectDirichletRows(*SM_Matrix_,Teuchos::ScalarTraits<magnitudeType>::eps(),/*count_twos_as_dirichlet=*/true));

      if (rowSumTol > 0.)
        Utilities::ApplyRowSumCriterion(*SM_Matrix_, rowSumTol, BCrows_);

      BCcols_.resize(D0_Matrix_->getColMap()->getLocalNumElements());
      BCdomain_.resize(D0_Matrix_->getDomainMap()->getLocalNumElements());
      Utilities::DetectDirichletColsAndDomains(*D0_Matrix_,BCrows_,BCcols_,BCdomain_);

      for (auto it = BCrows_.begin(); it != BCrows_.end(); ++it)
        if (*it)
          BCedgesLocal += 1;
      for (auto it = BCdomain_.begin(); it != BCdomain_.end(); ++it)
        if (*it)
          BCnodesLocal += 1;
    }

    MueLu_sumAll(SM_Matrix_->getRowMap()->getComm(), BCedgesLocal, BCedges_);
    MueLu_sumAll(SM_Matrix_->getRowMap()->getComm(), BCnodesLocal, BCnodes_);


    allEdgesBoundary_ = Teuchos::as<Xpetra::global_size_t>(BCedges_) >= D0_Matrix_->getRangeMap()->getGlobalNumElements();
    allNodesBoundary_ = Teuchos::as<Xpetra::global_size_t>(BCnodes_) >= D0_Matrix_->getDomainMap()->getGlobalNumElements();
  }


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell_Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::removeExplicitZeros(Teuchos::ParameterList &parameterList_,
                                                                                  RCP<Matrix> & D0_Matrix_,
                                                                                  RCP<Matrix> & SM_Matrix_,
                                                                                  RCP<Matrix> & M1_Matrix_,
                                                                                  RCP<Matrix> & Ms_Matrix_) {    

    bool defaultFilter = false;

    // Remove zero entries from D0 if necessary.
    // In the construction of the prolongator we use the graph of the
    // matrix, so zero entries mess it up.
    if (parameterList_.get<bool>("refmaxwell: filter D0", true) && D0_Matrix_->getLocalMaxNumRowEntries()>2) {
      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",D0_Matrix_);
      fineLevel.setlib(D0_Matrix_->getDomainMap()->lib());
      // We expect D0 to have entries +-1, so any threshold value will do.
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",1.0e-8,/*keepDiagonal=*/false,/*expectedNNZperRow=*/2));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      D0_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());

      // If D0 has too many zeros, maybe SM and M1 do as well.
      defaultFilter = true;
    }

    if (!M1_Matrix_.is_null() && parameterList_.get<bool>("refmaxwell: filter M1", defaultFilter)) {
      RCP<Vector> diag = VectorFactory::Build(M1_Matrix_->getRowMap());
      // find a reasonable absolute value threshold
      M1_Matrix_->getLocalDiagCopy(*diag);
      magnitudeType threshold = 1.0e-8 * diag->normInf();

      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",M1_Matrix_);
      fineLevel.setlib(M1_Matrix_->getDomainMap()->lib());
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",threshold,/*keepDiagonal=*/true));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      M1_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());
    }

    if (!Ms_Matrix_.is_null() && parameterList_.get<bool>("refmaxwell: filter Ms", defaultFilter)) {
      RCP<Vector> diag = VectorFactory::Build(Ms_Matrix_->getRowMap());
      // find a reasonable absolute value threshold
      Ms_Matrix_->getLocalDiagCopy(*diag);
      magnitudeType threshold = 1.0e-8 * diag->normInf();

      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",Ms_Matrix_);
      fineLevel.setlib(Ms_Matrix_->getDomainMap()->lib());
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",threshold,/*keepDiagonal=*/true));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      Ms_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());
    }

    if (!SM_Matrix_.is_null() && parameterList_.get<bool>("refmaxwell: filter SM", defaultFilter)) {
      RCP<Vector> diag = VectorFactory::Build(SM_Matrix_->getRowMap());
      // find a reasonable absolute value threshold
      SM_Matrix_->getLocalDiagCopy(*diag);
      magnitudeType threshold = 1.0e-8 * diag->normInf();

      Level fineLevel;
      fineLevel.SetFactoryManager(null);
      fineLevel.SetLevelID(0);
      fineLevel.Set("A",SM_Matrix_);
      fineLevel.setlib(SM_Matrix_->getDomainMap()->lib());
      RCP<ThresholdAFilterFactory> ThreshFact = rcp(new ThresholdAFilterFactory("A",threshold,/*keepDiagonal=*/true));
      fineLevel.Request("A",ThreshFact.get());
      ThreshFact->Build(fineLevel);
      SM_Matrix_ = fineLevel.Get< RCP<Matrix> >("A",ThreshFact.get());
    }

  }



  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void Maxwell_Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  setMatvecParams(Matrix& A, RCP<ParameterList> matvecParams) {
    RCP<const Import> xpImporter = A.getCrsGraph()->getImporter();
    if (!xpImporter.is_null())
      xpImporter->setDistributorParameters(matvecParams);
    RCP<const Export> xpExporter = A.getCrsGraph()->getExporter();
    if (!xpExporter.is_null())
      xpExporter->setDistributorParameters(matvecParams);
  }
  
  


#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  struct StratimikosWrapperImpl {
    static RCP<Thyra::PreconditionerBase<Scalar> > setupStratimikosPreconditioner(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A,
                                                                                  RCP<ParameterList> params) {
      throw std::runtime_error("setupStratimikosPreconditioner: Requires Scalar=double");
    }
  };

  template<class LocalOrdinal, class GlobalOrdinal, class Node>
  struct StratimikosWrapperImpl<double,LocalOrdinal,GlobalOrdinal,Node> {
    static RCP<Thyra::PreconditionerBase<double> > setupStratimikosPreconditioner(RCP<Xpetra::Matrix<double,LocalOrdinal,GlobalOrdinal,Node> > A,
                                                                                  RCP<ParameterList> params) {
      typedef double Scalar;
      
      // Build Thyra linear algebra objects
      RCP<const Thyra::LinearOpBase<Scalar> > thyraA = Xpetra::ThyraUtils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::toThyra(Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(A)->getCrsMatrix());
      
      Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
      typedef Thyra::PreconditionerFactoryBase<Scalar>                                     Base;
      typedef Thyra::MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> ImplMueLu;
      linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, ImplMueLu>(), "MueLu");
#ifdef HAVE_MUELU_IFPACK2
      // Register Ifpack2 as a Stratimikos preconditioner strategy.
      typedef Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Impl;
      linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
#endif
      
      linearSolverBuilder.setParameterList(params);
      
      // Build a new "solver factory" according to the previously specified parameter list.
      // RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);
      
      auto precFactory = Thyra::createPreconditioningStrategy(linearSolverBuilder);
      auto prec = precFactory->createPrec();
      
      precFactory->initializePrec(Thyra::defaultLinearOpSource(thyraA), prec.get(), Thyra::SUPPORT_SOLVE_UNSPECIFIED);
      
      return prec;
    }
  };


  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<Thyra::PreconditionerBase<Scalar> > 
  Maxwell_Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node>::setupStratimikosPreconditioner(RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A,
                                                                                        RCP<ParameterList> params) {
    return StratimikosWrapperImpl<SC,LO,GO,NO>::setupStratimikosPreconditioner(A,params);
  }




#endif

} // namespace

#define MUELU_MAXWELL_UTILS_SHORT
#endif //ifdef MUELU_MAXWELL_UTILS_DEF_HPP
