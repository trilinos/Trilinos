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

#ifndef _FROSCH_SCHWARZOPERATOR_DECL_HPP
#define _FROSCH_SCHWARZOPERATOR_DECL_HPP

#define FROSCH_ASSERT(A,S) if(!(A)) { std::cerr<<"Assertion failed. "<<S<<std::endl; std::cout.flush(); throw std::out_of_range("Assertion.");};

#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_Export.hpp>

#include <Teuchos_DefaultSerialComm.hpp>

#include <FROSch_DDInterface_def.hpp>
#include <FROSch_EntitySet_def.hpp>

#include <FROSch_CoarseSpace_def.hpp>
#include <FROSch_InterfacePartitionOfUnity_def.hpp>
#include <FROSch_LocalPartitionOfUnityBasis_def.hpp>

#include <FROSch_SubdomainSolver_def.hpp>

// TODO: Auf const überprüfen
// TODO: #ifndef überprüfen ??????


namespace FROSch {
    
    class Solver;
    
    template <class SC = Xpetra::Operator<>::scalar_type,
              class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
              class GO = typename Xpetra::Operator<SC,LO>::global_ordinal_type,
              class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class SchwarzOperator : public Xpetra::Operator<SC,LO,GO,NO> {
        
    public:
        
        typedef Teuchos::RCP<const Teuchos::Comm<int> > CommPtr;
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::RCP<const Map> ConstMapPtr;
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;
        typedef Teuchos::ArrayRCP<MapPtrVecPtr> MapPtrVecPtr2D;
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        typedef Teuchos::ArrayRCP<MultiVectorPtr> MultiVectorPtrVecPtr;
        
        typedef Xpetra::Import<LO,GO,NO> Importer;
        typedef Teuchos::RCP<Importer> ImporterPtr;
        
        typedef Xpetra::Export<LO,GO,NO> Exporter;
        typedef Teuchos::RCP<Exporter> ExporterPtr;
        typedef Teuchos::ArrayRCP<ExporterPtr> ExporterPtrVecPtr;
        
        typedef Teuchos::RCP<Teuchos::ParameterList> ParameterListPtr;
        
        typedef Teuchos::RCP<DDInterface<SC,LO,GO,NO> > DDInterfacePtr;
        
        typedef Teuchos::RCP<EntitySet<SC,LO,GO,NO> > EntitySetPtr;
        
        typedef Teuchos::RCP<CoarseSpace<SC,LO,GO,NO> > CoarseSpacePtr;
        
        typedef Teuchos::RCP<InterfaceEntity<SC,LO,GO,NO> > InterfaceEntityPtr;
        
        typedef Teuchos::RCP<InterfacePartitionOfUnity<SC,LO,GO,NO> > InterfacePartitionOfUnityPtr;
        
        typedef Teuchos::RCP<LocalPartitionOfUnityBasis<SC,LO,GO,NO> > LocalPartitionOfUnityBasisPtr;
        
        typedef Teuchos::RCP<SchwarzOperator<SC,LO,GO,NO> > SchwarzOperatorPtr;
        typedef Teuchos::Array<SchwarzOperatorPtr> SchwarzOperatorPtrVec;
        typedef Teuchos::ArrayRCP<SchwarzOperatorPtr> SchwarzOperatorPtrVecPtr;
        
        typedef Teuchos::RCP<SubdomainSolver<SC,LO,GO,NO> > SubdomainSolverPtr;
        
        typedef unsigned UN;
        typedef Teuchos::ArrayRCP<UN> UNVecPtr;
        
        typedef Teuchos::Array<LO> LOVec;
        typedef Teuchos::ArrayRCP<LO> LOVecPtr;
        typedef Teuchos::ArrayView<LO> LOVecView;
        typedef Teuchos::ArrayRCP<LOVecPtr> LOVecPtr2D;
        
        typedef Teuchos::Array<GO> GOVec;
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        typedef Teuchos::ArrayView<GO> GOVecView;
        typedef Teuchos::ArrayView<const GO> ConstGOVecView;
        typedef Teuchos::ArrayRCP<GOVecPtr> GOVecPtr2D;
        
        typedef Teuchos::Array<SC> SCVec;
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;
        typedef Teuchos::ArrayView<const SC> ConstSCVecView;
        
        typedef Teuchos::Array<bool> BoolVec;
        typedef Teuchos::ArrayRCP<bool> BoolVecPtr;
        
        
        SchwarzOperator(CommPtr comm);
        
        SchwarzOperator(CrsMatrixPtr k,
                        ParameterListPtr parameterList);
        
        virtual ~SchwarzOperator();
        
        virtual int initialize() = 0;
        
        virtual int compute() = 0;
        
        // Y = alpha * A^mode * X + beta * Y
        virtual void apply(const MultiVector &x,
                          MultiVector &y,
                          Teuchos::ETransp mode=Teuchos::NO_TRANS,
                          SC alpha=Teuchos::ScalarTraits<SC>::one(),
                          SC beta=Teuchos::ScalarTraits<SC>::zero()) const;
        
        virtual void apply(const MultiVector &x,
                          MultiVector &y,
                          bool usePreconditionerOnly,
                          Teuchos::ETransp mode=Teuchos::NO_TRANS,
                          SC alpha=Teuchos::ScalarTraits<SC>::one(),
                          SC beta=Teuchos::ScalarTraits<SC>::zero()) const = 0;
        
        virtual ConstMapPtr getDomainMap() const;
        
        virtual ConstMapPtr getRangeMap() const;
        
        virtual void describe(Teuchos::FancyOStream &out,
                              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;
        
        virtual std::string description() const = 0;
        
        bool isInitialized() const;
        
        bool isComputed() const;
        
        
    protected:
        
        CommPtr MpiComm_;
        CommPtr SerialComm_;
        
        CrsMatrixPtr K_;
        
        ParameterListPtr ParameterList_;
        
        bool Verbose_;
        
        bool IsInitialized_;
        bool IsComputed_;
        
    };
    
}

#endif
