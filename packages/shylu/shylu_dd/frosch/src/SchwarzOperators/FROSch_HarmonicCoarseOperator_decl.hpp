// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_HARMONICCOARSEOPERATOR_DECL_HPP
#define _FROSCH_HARMONICCOARSEOPERATOR_DECL_HPP

#include <FROSch_CoarseOperator_def.hpp>


namespace FROSch {

    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class HarmonicCoarseOperator : public CoarseOperator<SC,LO,GO,NO> {

    protected:

        using XMapPtr                   = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr              = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr             = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr;
        using XMapPtrVecPtr2D           = typename SchwarzOperator<SC,LO,GO,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D      = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;

        using XMatrixPtr                = typename SchwarzOperator<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr           = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtr;
        using ConstXMultiVectorPtr      = typename SchwarzOperator<SC,LO,GO,NO>::ConstXMultiVectorPtr;
        using XMultiVectorPtrVecPtr     = typename SchwarzOperator<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;

        using XCrsGraph                 = typename SchwarzOperator<SC,LO,GO,NO>::XCrsGraph;
        using GraphPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::GraphPtr;
        using ConstXCrsGraphPtr         = typename SchwarzOperator<SC,LO,GO,NO>::ConstXCrsGraphPtr;

        using ParameterListPtr          = typename SchwarzOperator<SC,LO,GO,NO>::ParameterListPtr;

        using TSerialDenseMatrixPtr     = typename SchwarzOperator<SC,LO,GO,NO>::TSerialDenseMatrixPtr;

        using TSerialQRDenseSolverPtr   = typename SchwarzOperator<SC,LO,GO,NO>::TSerialQRDenseSolverPtr;

        using CoarseSpacePtr            = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtr;
        using CoarseSpacePtrVecPtr      = typename SchwarzOperator<SC,LO,GO,NO>::CoarseSpacePtrVecPtr;

        using EntitySetPtr              = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtr;
        using EntitySetConstPtr         = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetConstPtr;
        using EntitySetPtrVecPtr        = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtrVecPtr;
        using EntitySetPtrConstVecPtr   = typename SchwarzOperator<SC,LO,GO,NO>::EntitySetPtrConstVecPtr;

        using InterfaceEntityPtr        = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtr;
        using InterfaceEntityPtrVec     = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtrVec;
        using InterfaceEntityPtrVecPtr  = typename SchwarzOperator<SC,LO,GO,NO>::InterfaceEntityPtrVecPtr;

        using SolverPtr                 = typename SchwarzOperator<SC,LO,GO,NO>::SolverPtr;
        using SolverFactoryPtr          = typename SchwarzOperator<SC,LO,GO,NO>::SolverFactoryPtr;

        using UN                        = typename SchwarzOperator<SC,LO,GO,NO>::UN;
        using UNVec                     = typename SchwarzOperator<SC,LO,GO,NO>::UNVec;
        using UNVecPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::UNVecPtr;
        using ConstUNVecView            = typename SchwarzOperator<SC,LO,GO,NO>::ConstUNVecView;

        using IntVec                    = typename SchwarzOperator<SC,LO,GO,NO>::IntVec;
        using IntVec2D                  = typename SchwarzOperator<SC,LO,GO,NO>::IntVec2D;

        using LOVec                     = typename SchwarzOperator<SC,LO,GO,NO>::LOVec;
        using LOVecPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr;
        using ConstLOVecView            = typename SchwarzOperator<SC,LO,GO,NO>::ConstLOVecView;
        using LOVecPtr2D                = typename SchwarzOperator<SC,LO,GO,NO>::LOVecPtr2D;

        using GOVec                     = typename SchwarzOperator<SC,LO,GO,NO>::GOVec;
        using GOVecView                 = typename SchwarzOperator<SC,LO,GO,NO>::GOVecView;
        using GOVec2D                   = typename SchwarzOperator<SC,LO,GO,NO>::GOVec2D;

        using SCVec                     = typename SchwarzOperator<SC,LO,GO,NO>::SCVec;
        using SCVecPtr                  = typename SchwarzOperator<SC,LO,GO,NO>::SCVecPtr;
        using ConstSCVecPtr             = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecPtr;
        using ConstSCVecView            = typename SchwarzOperator<SC,LO,GO,NO>::ConstSCVecView;

        using ConstBoolVecPtr           = typename SchwarzOperator<SC,LO,GO,NO>::ConstBoolVecPtr;

    public:

        HarmonicCoarseOperator(ConstXMatrixPtr k,
                               ParameterListPtr parameterList);

        virtual int initialize() = 0;

        XMapPtr assembleCoarseMap();

        ConstXMapPtr computeCoarseSpace(CoarseSpacePtr coarseSpace);

        #if defined(HAVE_XPETRA_TPETRA)
        template<class GOIndView, class ConstSCView, class SCView>
        struct CopyPhiViewFunctor
        {
            CopyPhiViewFunctor(UN j_in_, GOIndView indices_, ConstSCView data_in_, UN j_out_, SCView data_out_) :
            j_in(j_in_),
            indices (indices_),
            data_in (data_in_),
            j_out(j_out_),
            data_out(data_out_)
            {}

            KOKKOS_INLINE_FUNCTION
            void operator()(const int k) const {
                data_out(indices(k), j_out) = data_in(k, j_in);
            }

            UN          j_in;
            GOIndView   indices;
            ConstSCView data_in;
            UN          j_out;
            SCView      data_out;
        };


        struct ScaleTag {};
        struct CountNnzTag {};
        struct TotalNnzTag {};
        struct FillNzEntriesTag {};
        template<class indicesView, class SCView, class localRowMapType, class localMVBasisType, class RowptrType, class IndicesType, class ValuesType>
        struct detectLinearDependenciesFunctor
        {
            using Real = typename Teuchos::ScalarTraits<SC>::magnitudeType;
            using STS = Kokkos::ArithTraits<SC>;
            using RTS = Kokkos::ArithTraits<Real>;

            UN numRows;
            UN numCols;
            SCView scale;
            localMVBasisType localMVBasis;

            SC tresholdDropping;
            indicesView indicesGammaDofsAll;
            localRowMapType localRowMap;
            localRowMapType localRepeatedMap;

            RowptrType  Rowptr;
            IndicesType Indices;
            ValuesType  Values;

            // Constructor for ScaleTag
            detectLinearDependenciesFunctor(UN numRows_, UN numCols_, SCView scale_, localMVBasisType localMVBasis_) :
            numRows (numRows_),
            numCols (numCols_),
            scale (scale_),
            localMVBasis (localMVBasis_)
            {}

            // Constructor for CountNnzTag
            detectLinearDependenciesFunctor(UN numRows_, UN numCols_, localMVBasisType localMVBasis_, SC tresholdDropping_, 
                                            indicesView indicesGammaDofsAll_, localRowMapType localRowMap_, localRowMapType localRepeatedMap_,
                                            RowptrType Rowptr_) :
            numRows (numRows_),
            numCols (numCols_),
            scale (),
            localMVBasis (localMVBasis_),
            tresholdDropping (tresholdDropping_),
            indicesGammaDofsAll (indicesGammaDofsAll_),
            localRowMap (localRowMap_),
            localRepeatedMap (localRepeatedMap_),
            Rowptr (Rowptr_)
            {}

            // Constructor for FillNzEntriesTag
            detectLinearDependenciesFunctor(UN numRows_, UN numCols_, SCView scale_, localMVBasisType localMVBasis_,
                                            SC tresholdDropping_, indicesView indicesGammaDofsAll_,
                                            localRowMapType localRowMap_, localRowMapType localRepeatedMap_,
                                            RowptrType Rowptr_, IndicesType Indices_, ValuesType Values_) :
            numRows (numRows_),
            numCols (numCols_),
            scale (scale_),
            localMVBasis (localMVBasis_),
            tresholdDropping (tresholdDropping_),
            indicesGammaDofsAll (indicesGammaDofsAll_),
            localRowMap (localRowMap_),
            localRepeatedMap (localRepeatedMap_),
            Rowptr (Rowptr_),
            Indices (Indices_),
            Values (Values_)
            {}

            KOKKOS_INLINE_FUNCTION
            void operator()(const ScaleTag &, const int j) const {
                scale(j) = STS::zero();
                for (UN i = 0; i < numRows; i++) {
                    scale(j) += localMVBasis(i,j)*localMVBasis(i,j);
                }
                scale(j) = RTS::one()/RTS::sqrt(STS::abs(scale(j)));
            }

            KOKKOS_INLINE_FUNCTION
            void operator()(const CountNnzTag &, const int i) const {
                LO rowID = indicesGammaDofsAll[i];
                GO iGlobal = localRepeatedMap.getGlobalElement(rowID);
                LO iLocal = localRowMap.getLocalElement(iGlobal);
                if (iLocal!=-1) { // This should prevent duplicate entries on the interface
                    for (UN j=0; j<numCols; j++) {
                        SC valueTmp=localMVBasis(i,j);
                        if (fabs(valueTmp)>tresholdDropping) {
                            Rowptr(iLocal+1) ++;
                        }
                    }
                }
            }


            KOKKOS_INLINE_FUNCTION
            void operator()(const TotalNnzTag&, const size_t i, UN &lsum) const {
                 lsum += Rowptr[i];
            }

            KOKKOS_INLINE_FUNCTION
            void operator()(const FillNzEntriesTag &, const int i) const {
                LO rowID = indicesGammaDofsAll[i];
                GO iGlobal = localRepeatedMap.getGlobalElement(rowID);
                LO iLocal = localRowMap.getLocalElement(iGlobal);
                if (iLocal!=-1) { // This should prevent duplicate entries on the interface
                    UN nnz_i = Rowptr(iLocal);
                    for (UN j=0; j<numCols; j++) {
                        SC valueTmp=localMVBasis(i,j);
                        if (fabs(valueTmp)>tresholdDropping) {
                            Indices(nnz_i) = j; //localBasisMap.getGlobalElement(j);
                            Values(nnz_i) = valueTmp*scale(j);
                            nnz_i ++;
                        }
                    }
                }
            }
        };
        #endif

    protected:

        int intializeCoarseMap();

        int assembleInterfaceCoarseSpace();

        int addZeroCoarseSpaceBlock(ConstXMapPtr dofsMap);

        int computeVolumeFunctions(UN blockId,
                                   UN dimension,
                                   ConstXMapPtr nodesMap,
                                   ConstXMultiVectorPtr nodeList,
                                   EntitySetConstPtr interior);

        virtual XMultiVectorPtrVecPtr computeTranslations(UN blockId,
                                                          EntitySetConstPtr entitySet);

        virtual XMultiVectorPtrVecPtr computeRotations(UN blockId,
                                                       UN dimension,
                                                       ConstXMultiVectorPtr nodeList,
                                                       EntitySetConstPtr entitySet,
                                                       UN discardRotations = 0);

        virtual LOVecPtr detectLinearDependencies(GOVecView indicesGammaDofsAll,
                                                        ConstXMapPtr rowMap,
                                                        ConstXMapPtr rangeMap,
                                                        ConstXMapPtr repeatedMap,
                                                        SC tresholdDropping,
                                                        SC tresholdOrthogonalization);

        virtual XMultiVectorPtr computeExtensions(ConstXMapPtr localMap,
                                                  GOVecView indicesGammaDofsAll,
                                                  GOVecView indicesIDofsAll,
                                                  XMatrixPtr kII,
                                                  XMatrixPtr kIGamma);

        virtual ConstXMatrixPtr removeCouplingBetweenDofs(ConstXMatrixPtr matrix,
                                                          ConstXMapPtr map,
                                                          TwoDArray<int> &couplingIDsToRemove);

        virtual int buildElementNodeList();

        virtual int buildGlobalGraph(Teuchos::RCP<DDInterface<SC,LO,GO,NO> > theDDInterface_);

        virtual int buildCoarseGraph();

        virtual void extractLocalSubdomainMatrix_Symbolic();

        SolverPtr ExtensionSolver_;

        CoarseSpacePtrVecPtr InterfaceCoarseSpaces_ = CoarseSpacePtrVecPtr(0);
        CoarseSpacePtr AssembledInterfaceCoarseSpace_;

        UNVecPtr Dimensions_ = UNVecPtr(0);
        UNVecPtr DofsPerNode_ = UNVecPtr(0);

        LOVecPtr2D GammaDofs_ = LOVecPtr2D(0);
        LOVecPtr2D IDofs_ = LOVecPtr2D(0);

        ConstXMapPtrVecPtr2D DofsMaps_ = ConstXMapPtrVecPtr2D(0); // notwendig??

        ConstXMapPtr KRowMap_; // AH: What is this?

        UN NumberOfBlocks_ = 0;

        UN MaxNumNeigh_ = 0;
    };

}

#endif
