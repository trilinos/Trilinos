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

#ifndef _FROSCH_SCHWARZPRECONDITIONER_DECL_HPP
#define _FROSCH_SCHWARZPRECONDITIONER_DECL_HPP

#include <Xpetra_Operator.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include <FROSch_SumOperator_def.hpp>
#include <FROSch_AlgebraicOverlappingOperator_def.hpp>
#include <FROSch_GDSWCoarseOperator_def.hpp>
#include <FROSch_RGDSWCoarseOperator_def.hpp>
#include <FROSch_IPOUHarmonicCoarseOperator_def.hpp>

namespace FROSch {        
    
    template <class SC = Xpetra::Operator<>::scalar_type,
    class LO = typename Xpetra::Operator<SC>::local_ordinal_type,
    class GO = typename Xpetra::Operator<SC, LO>::global_ordinal_type,
    class NO = typename Xpetra::Operator<SC,LO,GO>::node_type>
    class SchwarzPreconditioner : public Xpetra::Operator<SC,LO,GO,NO> {
        
    public:
        
        typedef Teuchos::RCP<const Teuchos::Comm<int> > CommPtr;
        
        typedef Xpetra::Map<LO,GO,NO> Map;
        typedef Teuchos::RCP<Map> MapPtr;
        typedef Teuchos::RCP<const Map> ConstMapPtr;
        typedef Teuchos::ArrayRCP<MapPtr> MapPtrVecPtr;
        
        typedef Xpetra::Matrix<SC,LO,GO,NO> CrsMatrix;
        typedef Teuchos::RCP<CrsMatrix> CrsMatrixPtr;
        
        typedef Xpetra::MultiVector<SC,LO,GO,NO> MultiVector;
        typedef Teuchos::RCP<MultiVector> MultiVectorPtr;
        
        typedef Teuchos::RCP<Teuchos::ParameterList> ParameterListPtr;
        
        typedef Teuchos::RCP<SumOperator<SC,LO,GO,NO> > SumOperatorPtr;
        typedef Teuchos::RCP<OverlappingOperator<SC,LO,GO,NO> > OverlappingOperatorPtr;
        typedef Teuchos::RCP<AlgebraicOverlappingOperator<SC,LO,GO,NO> > AlgebraicOverlappingOperatorPtr;
        typedef Teuchos::RCP<CoarseOperator<SC,LO,GO,NO> > CoarseOperatorPtr;
        typedef Teuchos::RCP<GDSWCoarseOperator<SC,LO,GO,NO> > GDSWCoarseOperatorPtr;
        typedef Teuchos::RCP<RGDSWCoarseOperator<SC,LO,GO,NO> > RGDSWCoarseOperatorPtr;
        typedef Teuchos::RCP<IPOUHarmonicCoarseOperator<SC,LO,GO,NO> > IPOUHarmonicCoarseOperatorPtr;
        
        typedef unsigned UN;
        
        typedef Teuchos::ArrayRCP<GO> GOVecPtr;
        
        typedef Teuchos::ArrayRCP<SC> SCVecPtr;

        
        SchwarzPreconditioner(ParameterListPtr parameterList,
                              CommPtr comm);
        
        virtual ~SchwarzPreconditioner();
        
        virtual int initialize(bool useDefaultParameters = true) = 0;
        
        virtual int compute() = 0;
        
        // Y = alpha * A^mode * X + beta * Y
        virtual void apply(const MultiVector &X,
                           MultiVector &Y,
                           Teuchos::ETransp mode=Teuchos::NO_TRANS,
                           SC alpha=Teuchos::ScalarTraits<SC>::one(),
                           SC beta=Teuchos::ScalarTraits<SC>::zero()) const = 0;
        
        virtual ConstMapPtr getDomainMap() const = 0;
        
        virtual ConstMapPtr getRangeMap() const = 0;
        
        virtual void describe(Teuchos::FancyOStream &out,
                              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const = 0;
        
        virtual std::string description() const = 0;
        
        bool isInitialized() const;
        
        bool isComputed() const;
        
        
    protected:
        
        CommPtr MpiComm_;
        
        ParameterListPtr ParameterList_;
        
        bool UseTranspose_;
        bool IsInitialized_;
        bool IsComputed_;
        bool Verbose_;
        
    };
    
}

#endif
