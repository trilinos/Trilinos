// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_TWOLEVELBLOCKPRECONDITIONER_DECL_HPP
#define _FROSCH_TWOLEVELBLOCKPRECONDITIONER_DECL_HPP

#include <FROSch_OneLevelPreconditioner_def.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC = double,
              class LO = int,
              class GO = DefaultGlobalOrdinal,
              class NO = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
    class TwoLevelBlockPreconditioner : public OneLevelPreconditioner<SC,LO,GO,NO> {

    protected:

        using XMapPtr                             = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtr;
        using ConstXMapPtr                        = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtr;
        using XMapPtrVecPtr                       = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtrVecPtr;
        using ConstXMapPtrVecPtr                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtrVecPtr;
        using XMapPtrVecPtr2D                     = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMapPtrVecPtr2D;
        using ConstXMapPtrVecPtr2D                = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMapPtrVecPtr2D;

        using XMatrixPtr                          = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMatrixPtr;
        using ConstXMatrixPtr                     = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMatrixPtr;

        using XMultiVectorPtr                     = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVectorPtr;
        using XMultiVectorPtrVecPtr               = typename SchwarzPreconditioner<SC,LO,GO,NO>::XMultiVectorPtrVecPtr;
        using ConstXMultiVectorPtrVecPtr          = typename SchwarzPreconditioner<SC,LO,GO,NO>::ConstXMultiVectorPtrVecPtr;

        using ParameterListPtr                    = typename SchwarzPreconditioner<SC,LO,GO,NO>::ParameterListPtr;

        using AlgebraicOverlappingOperatorPtr     = typename SchwarzPreconditioner<SC,LO,GO,NO>::AlgebraicOverlappingOperatorPtr;
        using CoarseOperatorPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::CoarseOperatorPtr;
        using GDSWCoarseOperatorPtr               = typename SchwarzPreconditioner<SC,LO,GO,NO>::GDSWCoarseOperatorPtr;
        using RGDSWCoarseOperatorPtr              = typename SchwarzPreconditioner<SC,LO,GO,NO>::RGDSWCoarseOperatorPtr;
        using IPOUHarmonicCoarseOperatorPtr       = typename SchwarzPreconditioner<SC,LO,GO,NO>::IPOUHarmonicCoarseOperatorPtr;

        using DofOrderingVecPtr                   = typename SchwarzPreconditioner<SC,LO,GO,NO>::DofOrderingVecPtr;

        using UN                                  = typename SchwarzPreconditioner<SC,LO,GO,NO>::UN;
        using UNVecPtr                            = typename SchwarzPreconditioner<SC,LO,GO,NO>::UNVecPtr;

        using LOVecPtr                            = typename SchwarzPreconditioner<SC,LO,GO,NO>::LOVecPtr;

        using GOVec                               = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVec;
        using GOVec2D                             = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVec2D;
        using GOVecPtr                            = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVecPtr;
        using GOVecPtr2D                          = typename SchwarzPreconditioner<SC,LO,GO,NO>::GOVecPtr2D;

    public:

        TwoLevelBlockPreconditioner(ConstXMatrixPtr k,
                                    ParameterListPtr parameterList);

        int initialize(UN dimension,
                       UNVecPtr dofsPerNodeVec,
                       DofOrderingVecPtr dofOrderingVec,
                       int overlap = -1,
                       ConstXMapPtrVecPtr repeatedMapVec = null,
                       ConstXMultiVectorPtrVecPtr nullSpaceBasisVec = null,
                       ConstXMultiVectorPtrVecPtr nodeListVec = null,
                       ConstXMapPtrVecPtr2D dofsMapsVec = null,
                       GOVecPtr2D dirichletBoundaryDofsVec = null);

        int initialize(UN dimension,
                       UNVecPtr dofsPerNodeVec,
                       DofOrderingVecPtr dofOrderingVec,
                       int overlap = -1,
                       ConstXMultiVectorPtrVecPtr nodeListVec = null,
                       ConstXMapPtrVecPtr repeatedMapVec = null,
                       ConstXMultiVectorPtrVecPtr nullSpaceBasisVec = null,
                       ConstXMapPtrVecPtr2D dofsMapsVec = null,
                       GOVecPtr2D dirichletBoundaryDofsVec = null);

        int compute();

        void describe(FancyOStream &out,
                      const EVerbosityLevel verbLevel=Describable::verbLevel_default) const;

        string description() const;

        int resetMatrix(ConstXMatrixPtr &k);

        int preApplyCoarse(XMultiVectorPtr &x,
                           XMultiVectorPtr &y);

    protected:

        CoarseOperatorPtr CoarseOperator_;

    private:

        /**
         * @brief Construct merged, repeated DoF map for each block.
         *
         * For block i, repeatedMapVec[i] contains all DoFs on this rank's overlapping
         * subdomain. It combines the component maps dofsMapsVec[i][j] for
         * j = 0, ..., dofsPerNodeVec[i]-1 (e.g., displacement components x, y, z).
         * These maps are repeated: the same global DoF may occur on multiple ranks.
         * AssembleSubdomainMap merges the components by node.
         *
         * The component maps use the global IDs of the full system.
         * The output maps instead use per-block numbering: the cumulative offset of
         * preceding blocks is subtracted after assembly.
         *
         * If neither repeatedMapVec nor dofsMapsVec is available, for a single-block system,
         * repeatedMapVec is constructed from the graph and row map of K_.
         * - Building a repeated map from a unique map, is not a bijevtive operation, meaning
         *   that different repeated maps can be built from a unique maps and vice versa.
         *   This process involves communication between subdomains.
         * If only dofsMapsVec is available, repeatedMapVec is built from that.
         * - This is a bijective operation, so no information is lost.
         *   This operation's communication is local to the subdomains.
         * If repeatedMapVec is available, no modification is done.
         *
         * @param[in,out] repeatedMapVec One repeated DoF map per block. Existing maps
         *     are preserved.
         * @param[in] dofsMapsVec Component DoF maps indexed by [block][component].
         *     Required for each missing output map in a system with multiple blocks.
         * @param[in] dofsPerNodeVec Number of DoFs per node in each block;
         *     its length specifies the number of blocks.
         */
        void buildRepeatedDoFMapPerBlock(ConstXMapPtrVecPtr& repeatedMapVec,
                                         const ConstXMapPtrVecPtr2D& dofsMapsVec,
                                         const UNVecPtr& dofsPerNodeVec);

        void determineDirichletBoundaryDofs(GOVecPtr2D& dirichletBoundaryDofsVec,
                                            ConstXMapPtrVecPtr repeatedMapVec,
                                            ConstXMapPtr repeatedMap,
                                            ConstXMapPtrVecPtr2D dofsMapsVec);

    };

}

#endif
