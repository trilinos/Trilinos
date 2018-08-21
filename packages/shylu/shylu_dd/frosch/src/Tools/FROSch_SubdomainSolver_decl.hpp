//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

#ifndef _FROSCH_SUBDOMAINSOLVER_DECL_hpp
#define _FROSCH_SUBDOMAINSOLVER_DECL_hpp

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include "Epetra_LinearProblem.h"

#include "Amesos_ConfigDefs.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

#include "Amesos2.hpp"

#include <BelosXpetraAdapterOperator.hpp>
#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>


//#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

namespace FROSch {
    
    template <class SC,
    class LO ,
    class GO ,
    class NO >
    class OneLevelPreconditioner;
    
    template <class SC = typename Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC, LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC, LO, GO>::node_type>
    class SubdomainSolver : public Xpetra::Operator<SC,LO,GO,NO> {
                
    public:
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::RCP<const Map> ConstMapPtr;
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;

        typedef Teuchos::ArrayRCP<GO> GOVecPtr;

        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        typedef Epetra_CrsMatrix EpetraCrsMatrix;
        typedef Teuchos::RCP<EpetraCrsMatrix> EpetraCrsMatrixPtr;
        typedef Tpetra::CrsMatrix<SC,LO,GO,NO> TpetraCrsMatrix;
        typedef Teuchos::RCP<TpetraCrsMatrix> TpetraCrsMatrixPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        typedef Teuchos::RCP<const MultiVector> ConstMultiVectorPtr;
        typedef Epetra_MultiVector EpetraMultiVector;
        typedef Teuchos::RCP<EpetraMultiVector> EpetraMultiVectorPtr;
        typedef Tpetra::MultiVector<SC,LO,GO,NO> TpetraMultiVector;
        typedef Teuchos::RCP<TpetraMultiVector> TpetraMultiVectorPtr;
        
        typedef Teuchos::RCP<Teuchos::ParameterList> ParameterListPtr;
        
        typedef Teuchos::RCP<Epetra_LinearProblem> LinearProblemPtr;
        
        typedef Teuchos::RCP<Amesos_BaseSolver> AmesosSolverPtr;
        
        typedef Teuchos::RCP<Amesos2::Solver<EpetraCrsMatrix,EpetraMultiVector> > Amesos2SolverEpetraPtr;
        typedef Teuchos::RCP<Amesos2::Solver<TpetraCrsMatrix,TpetraMultiVector> > Amesos2SolverTpetraPtr;
        
        typedef Teuchos::RCP<MueLu::HierarchyManager<SC,LO,GO,NO> > MueLuFactoryPtr;
        typedef Teuchos::RCP<MueLu::Hierarchy<SC,LO,GO,NO> > MueLuHierarchyPtr;
        
        SubdomainSolver(CrsMatrixPtr k,
                        ParameterListPtr parameterList,
                        GOVecPtr blockCoarseSize=Teuchos::null);
        
        virtual ~SubdomainSolver();
        
        virtual int initialize();
        
        virtual int compute();
        
        // Y = alpha * A^mode * X + beta * Y
        virtual void apply(const MultiVector &x,
                           MultiVector &y,
                           Teuchos::ETransp mode=Teuchos::NO_TRANS,
                           SC alpha=Teuchos::ScalarTraits<SC>::one(),
                           SC beta=Teuchos::ScalarTraits<SC>::zero()) const;
        
        virtual ConstMapPtr getDomainMap() const;
        
        virtual ConstMapPtr getRangeMap() const;
        
        virtual void describe(Teuchos::FancyOStream &out,
                              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
        
        virtual std::string description() const;
        
        bool isInitialized() const;
        
        bool isComputed() const;
        
    protected:
        
        CrsMatrixPtr K_;
        
        ParameterListPtr ParameterList_;
        
        LinearProblemPtr EpetraLinearProblem_;
        
        AmesosSolverPtr AmesosSolver_;
        
        Amesos2SolverEpetraPtr Amesos2SolverEpetra_;
        Amesos2SolverTpetraPtr Amesos2SolverTpetra_;
        
        MueLuFactoryPtr MueLuFactory_;
        MueLuHierarchyPtr MueLuHierarchy_;
        
        Teuchos::RCP<Belos::LinearProblem<SC,Xpetra::MultiVector<SC,LO,GO,NO>,Belos::OperatorT<Xpetra::MultiVector<SC,LO,GO,NO> > > >  BelosLinearProblem_;
        Teuchos::RCP<Belos::SolverManager<SC,Xpetra::MultiVector<SC,LO,GO,NO>,Belos::OperatorT<Xpetra::MultiVector<SC,LO,GO,NO> > > > BelosSolverManager_;
        
        bool IsInitialized_;
        bool IsComputed_;        
    };
    
}

#endif
