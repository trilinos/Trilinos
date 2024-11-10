// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_COARSESPACE_DEF_HPP
#define _FROSCH_COARSESPACE_DEF_HPP

#include <FROSch_CoarseSpace_decl.hpp>


namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    CoarseSpace<SC,LO,GO,NO>::CoarseSpace(CommPtr mpiComm,
                                          CommPtr serialComm) :
    MpiComm_ (mpiComm),
    SerialComm_ (serialComm)
    {

    }

    // Will man Informationen über die Subspaces als strings reingeben?
    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::addSubspace(ConstXMapPtr subspaceBasisMap,
                                              ConstXMapPtr subspaceBasisMapUnique,
                                              ConstXMultiVectorPtr subspaceBasis,
                                              UN offset)
    {
        FROSCH_ASSERT(!subspaceBasisMap.is_null(),"FROSch::CoarseSpace: subspaceBasisMap.is_null()");
        if (!subspaceBasis.is_null()) {
            FROSCH_ASSERT(subspaceBasis->getNumVectors()==subspaceBasisMap->getLocalNumElements(),"FROSch::CoarseSpace: subspaceBasis->getNumVectors()!=subspaceBasisMap->getLocalNumElements()");
        } else {
            FROSCH_ASSERT(subspaceBasisMap->getLocalNumElements()==0,"FROSch::CoarseSpace: subspaceBasisMap->getLocalNumElements()!=0");
        }

        UnassembledBasesMaps_.push_back(subspaceBasisMap);
        UnassembledBasesMapsUnique_.push_back(subspaceBasisMapUnique);
        UnassembledSubspaceBases_.push_back(subspaceBasis);
        Offsets_.push_back(offset);
        LocalSubspacesSizes_.push_back(subspaceBasisMap->getLocalNumElements());

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::assembleCoarseSpace()
    {
        FROSCH_ASSERT(UnassembledBasesMaps_.size()>0,"FROSch::CoarseSpace: UnassembledBasesMaps_.size()==0");
        FROSCH_ASSERT(UnassembledBasesMapsUnique_.size()>0,"FROSch::CoarseSpace: UnassembledBasesMapsUnique_.size()==0");
        FROSCH_ASSERT(UnassembledSubspaceBases_.size()>0,"FROSch::CoarseSpace: UnassembledSubspaceBases_.size()==0");

        UN itmp = 0;
        LOVecPtr2D partMappings;

        // BasisMap
        AssembledBasisMap_ = AssembleMaps(UnassembledBasesMaps_(),partMappings);

        // BasisMapUnique - First, we check if any of the unassembled unique maps is null. In case, we re-build a unique map
        bool buildUniqueMap = false;
        UN i=0;
        while (!buildUniqueMap && i<UnassembledBasesMapsUnique_.size()) {
            buildUniqueMap = UnassembledBasesMapsUnique_[i].is_null();
            i++;
        }
        int buildUniqueMapMax = 0;
        reduceAll(*this->MpiComm_,REDUCE_MAX,int(buildUniqueMap),ptr(&buildUniqueMapMax));

        if (buildUniqueMapMax>0) {
            FROSCH_NOTIFICATION("FROSch::CoarseSpace",this->MpiComm_->getRank()==0,"We re-build a unique map of AssembledBasisMap_.");
            AssembledBasisMapUnique_ = BuildUniqueMap<LO,GO,NO>(AssembledBasisMap_);
        } else {
            AssembledBasisMapUnique_ = AssembleMaps(UnassembledBasesMapsUnique_(),partMappings);
        }
        FROSCH_ASSERT(AssembledBasisMap_->getMaxAllGlobalIndex()==AssembledBasisMapUnique_->getMaxAllGlobalIndex(),"FROSch::CoarseSpace: AssembledBasisMap_->getMaxAllGlobalIndex()!=AssembledBasisMapUnique_->getMaxAllGlobalIndex()");
        FROSCH_ASSERT(AssembledBasisMap_->getMinAllGlobalIndex()==AssembledBasisMapUnique_->getMinAllGlobalIndex(),"FROSch::CoarseSpace: AssembledBasisMap_->getMinAllGlobalIndex()!=AssembledBasisMapUnique_->getMinAllGlobalIndex()");
        FROSCH_ASSERT(GO(AssembledBasisMapUnique_->getGlobalNumElements())==GO(AssembledBasisMapUnique_->getMaxAllGlobalIndex()+1),"FROSch::CoarseSpace: AssembledBasisMapUnique_->getGlobalNumElements()!=(AssembledBasisMapUnique_->getMaxAllGlobalIndex()+1)");

        // Basis
        if (!AssembledBasisMap_.is_null()) {
            if (AssembledBasisMap_->getGlobalNumElements()>0) { // AH 02/12/2019: Is this the right condition? Seems to work for now...
                LO totalSize = -1;
                for (UN i=0; i<UnassembledSubspaceBases_.size(); i++) {
                    if (!UnassembledSubspaceBases_[i].is_null()) totalSize = max(totalSize,LO(UnassembledSubspaceBases_[i]->getLocalLength()+Offsets_[i]));
                }
                XMapPtr serialMap = MapFactory<LO,GO,NO>::Build(AssembledBasisMap_->lib(),totalSize,0,this->SerialComm_);

                AssembledBasis_ = MultiVectorFactory<SC,LO,GO,NO >::Build(serialMap,AssembledBasisMap_->getLocalNumElements());
                #if defined(HAVE_XPETRA_TPETRA)
                if (AssembledBasis_->getMap()->lib() == UseTpetra) {
                    using execution_space = typename XMap::local_map_type::execution_space;
                    // Xpetra wrapper for Tpetra MV
                    auto assembledXTpetraMVector = rcp_dynamic_cast<TpetraMultiVector<SC,LO,GO,NO>>(AssembledBasis_, true);
                    // Tpetra MV
                    auto assembledTpetraMVector = assembledXTpetraMVector->getTpetra_MultiVector();
                    // Kokkos View
                    auto assembledView = assembledTpetraMVector->getLocalViewDevice(Tpetra::Access::ReadWrite);
                    auto assembledCols = Tpetra::getMultiVectorWhichVectors(*  assembledTpetraMVector);

                    UN itmp = 0;
                    for (UN i=0; i<UnassembledSubspaceBases_.size(); i++) {
                        if (!UnassembledSubspaceBases_[i].is_null()) {
                            const UN Offset_i = Offsets_[i];
                            const UN NumVectors_i = UnassembledSubspaceBases_[i]->getNumVectors();
                            const UN LocalLength_i = UnassembledSubspaceBases_[i]->getLocalLength();

                            FROSCH_ASSERT(NumVectors_i+itmp <= AssembledBasis_->getNumVectors(),"FROSch::CoarseSpace: NumVectors_i+itmp <= AssembledBasis_->getNumVectors()");
                            FROSCH_ASSERT(LocalLength_i+Offsets_[i] <= AssembledBasis_->getLocalLength(),"FROSch::CoarseSpace: LocalLength_i+Offsets_[i] <= AssembledBasis_");

                            Kokkos::RangePolicy<execution_space> policy (0, LocalLength_i);
                            // Xpetra wrapper for Tpetra MV
                            auto unassembledXTpetraMVector = rcp_dynamic_cast<const TpetraMultiVector<SC,LO,GO,NO>>(UnassembledSubspaceBases_[i], true);
                            // Tpetra MV
                            auto unassembledTpetraMVector = unassembledXTpetraMVector->getTpetra_MultiVector();
                            // Views
                            auto unassembledView = unassembledTpetraMVector->getLocalViewDevice(Tpetra::Access::ReadOnly);
                            auto unassembledCols = Tpetra::getMultiVectorWhichVectors(*unassembledTpetraMVector);
                            for (UN j=0; j < NumVectors_i; j++) {
                                int col_in  = unassembledTpetraMVector->isConstantStride() ? j      : unassembledCols[j];
                                int col_out =   assembledTpetraMVector->isConstantStride() ? j+itmp :   assembledCols[j+itmp];
                                Kokkos::parallel_for(
                                    "FROSch_CoarseSpace::assembleCoarseSpace", policy,
                                    KOKKOS_LAMBDA(const UN k) {
                                        assembledView(k+Offset_i, col_out) = unassembledView(k, col_in);
                                    });
                            }
                            itmp += NumVectors_i;
                        }
                    }
                    Kokkos::fence();
                } else
                #endif
                {
                    for (UN i=0; i<UnassembledSubspaceBases_.size(); i++) {
                        if (!UnassembledSubspaceBases_[i].is_null()) {
                            for (UN j=0; j<UnassembledSubspaceBases_[i]->getNumVectors(); j++) {
                                ConstSCVecPtr unassembledSubspaceBasesData = UnassembledSubspaceBases_[i]->getData(j);
                                for (UN k=0; k<UnassembledSubspaceBases_[i]->getLocalLength(); k++) {
                                    FROSCH_ASSERT(itmp<AssembledBasis_->getNumVectors(),"FROSch::CoarseSpace: itmp>=AssembledBasis_->getNumVectors()");
                                    FROSCH_ASSERT(k+Offsets_[i]<AssembledBasis_->getLocalLength(),"FROSch::CoarseSpace: k+Offsets_[i]>=AssembledBasis_->getLocalLength()");
                                    AssembledBasis_->replaceLocalValue(k+Offsets_[i],itmp,unassembledSubspaceBasesData[k]);
                                }
                                itmp++;
                            }
                        }
                    }
                }
            }
        }

        ConstXMapPtrVec emptyVec1;
        UnassembledBasesMaps_.swap(emptyVec1);

        ConstXMapPtrVec emptyVec2;
        UnassembledBasesMapsUnique_.swap(emptyVec2);

        ConstXMultiVectorPtrVec emtpyVec3;
        UnassembledSubspaceBases_.swap(emtpyVec3);

        LOVec emptyVec4;
        Offsets_.swap(emptyVec4);

        UnassembledBasesMaps_.push_back(AssembledBasisMap_);
        UnassembledBasesMapsUnique_.push_back(AssembledBasisMapUnique_);
        UnassembledSubspaceBases_.push_back(AssembledBasis_);
        Offsets_.push_back(0);

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::buildGlobalBasisMatrix(ConstXMapPtr rowMap,
                                                         ConstXMapPtr rangeMap,
                                                         ConstXMapPtr repeatedMap,
                                                         SC tresholdDropping)
    {
        FROSCH_ASSERT(!AssembledBasisMap_.is_null(),"FROSch::CoarseSpace: AssembledBasisMap_.is_null().");
        FROSCH_ASSERT(!AssembledBasis_.is_null(),"FROSch::CoarseSpace: AssembledBasis_.is_null().");

#if defined(HAVE_XPETRA_TPETRA)
        if (rowMap->lib() == UseTpetra) {
            UN numRows = AssembledBasis_->getLocalLength();
            UN numCols = AssembledBasis_->getNumVectors();

            using crsmat_type  = typename XMatrix::local_matrix_type;
            using graph_type   = typename crsmat_type::StaticCrsGraphType;
            using rowptr_type  = typename graph_type::row_map_type::non_const_type;
            using indices_type = typename graph_type::entries_type::non_const_type;
            using values_type  = typename crsmat_type::values_type::non_const_type;

            using execution_space = typename XMap::local_map_type::execution_space;
            Kokkos::RangePolicy<execution_space> policy_col (0, numCols);
            Kokkos::RangePolicy<execution_space> policy_row (0, numRows);

            auto repeatedLocalMap = repeatedMap->getLocalMap();
            auto rowLocalMap = rowMap->getLocalMap();
            auto AssembledBasisView = AssembledBasis_->getDeviceLocalView(Access::ReadOnly);

            // count number of nonzeros per row
            UN numLocalRows = rowMap->getLocalNumElements();
            rowptr_type Rowptr (Kokkos::ViewAllocateWithoutInitializing("Rowptr"), numLocalRows+1);
            Kokkos::deep_copy(Rowptr, 0);
            Kokkos::parallel_for(
                "FROSch_CoarseSpace::countGlobalBasisMatrix", policy_row,
                KOKKOS_LAMBDA(const UN i) {
                    LO iD = repeatedLocalMap.getGlobalElement(i);
                    LO lo = rowLocalMap.getLocalElement(iD);
                    if (lo != -1) {
                        for (UN j=0; j<numCols; j++) {
                            SC valueTmp=AssembledBasisView(i, j);
                            if (fabs(valueTmp) > tresholdDropping) {
                                Rowptr[lo+1] ++;
                            }
                        }
                    }
                });
            Kokkos::fence();

            // cout nnz
            UN nnz = 0; //Rowptr[numLocalRows];
            Kokkos::parallel_reduce("FROSch_CoarseSpace::fillGlobalBasisMatrix:nnz", 1+numLocalRows,
              KOKKOS_LAMBDA(const int &i, UN &lsum) { lsum += Rowptr[i]; },
              nnz);
            Kokkos::fence();

            // make it into offsets
#if KOKKOSKERNELS_VERSION >= 40199
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>
              (1+numLocalRows, Rowptr);
#else
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
              (1+numLocalRows, Rowptr);
#endif

            // fill into the local matrix
            indices_type Indices (Kokkos::ViewAllocateWithoutInitializing("Indices"), nnz);
            values_type  Values  (Kokkos::ViewAllocateWithoutInitializing("Values"),  nnz);
            auto AssembledBasisLocalMap = AssembledBasisMap_->getLocalMap();
            Kokkos::parallel_for(
                "FROSch_CoarseSpace::fillGlobalBasisMatrix", policy_row,
                KOKKOS_LAMBDA(const UN i) {
                    LO iD = repeatedLocalMap.getGlobalElement(i);
                    LO lo = rowLocalMap.getLocalElement(iD);
                    if (lo != -1) { // This should prevent duplicate entries on the interface
                        UN nnz_i = Rowptr[lo];
                        for (UN j=0; j<numCols; j++) {
                            SC valueTmp=AssembledBasisView(i, j);
                            if (fabs(valueTmp) > tresholdDropping) {
                                Values[nnz_i] = valueTmp;
                                Indices[nnz_i] = j;

                                nnz_i ++;
                            }
                        }
                    }
                });
            Kokkos::fence();

            // create local matrix
            graph_type crsgraph (Indices, Rowptr);
            crsmat_type LocalBasisMatrix = crsmat_type ("CrsMatrix", numCols, Values, crsgraph);

            /// build into GlobalBasisMatrix
            Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp (new Teuchos::ParameterList());
            params->set("sorted", false);
            GlobalBasisMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(LocalBasisMatrix,
                                                                   rowMap, AssembledBasisMap_,
                                                                   AssembledBasisMapUnique_, rangeMap,
                                                                   params);
        } else
#endif
        {
            if (rowMap->lib()==UseEpetra) {
                GlobalBasisMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(rowMap,AssembledBasisMap_->getLocalNumElements()); // Nonzeroes abhängig von dim/dofs!!!
                LO iD;
                SC valueTmp;
                for (UN i=0; i<AssembledBasis_->getLocalLength(); i++) {
                    GOVec indices;
                    SCVec values;
                    for (UN j=0; j<AssembledBasis_->getNumVectors(); j++) {
                        valueTmp=AssembledBasis_->getData(j)[i];
                        if (fabs(valueTmp)>tresholdDropping) {
                            indices.push_back(AssembledBasisMap_->getGlobalElement(j));
                            values.push_back(valueTmp);
                        }
                    }
                    iD = repeatedMap->getGlobalElement(i);

                    if (rowMap->getLocalElement(iD)!=-1) { // This should prevent duplicate entries on the interface
                        GlobalBasisMatrix_->insertGlobalValues(iD,indices(),values());
                    }
                }
            } else {
                GlobalBasisMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(rowMap,AssembledBasisMap_,AssembledBasisMap_->getLocalNumElements()); // Nonzeroes abhängig von dim/dofs!!!
                LO iD;
                SC valueTmp;
                for (UN i=0; i<AssembledBasis_->getLocalLength(); i++) {
                    LOVec indices;
                    SCVec values;
                    for (UN j=0; j<AssembledBasis_->getNumVectors(); j++) {
                        valueTmp=AssembledBasis_->getData(j)[i];
                        if (fabs(valueTmp)>tresholdDropping) {
                            indices.push_back(j);
                            values.push_back(valueTmp);
                        }
                    }
                    iD = rowMap->getLocalElement(repeatedMap->getGlobalElement(i));

                    if (iD!=-1) { // This should prevent duplicate entries on the interface
                        GlobalBasisMatrix_->insertLocalValues(iD,indices(),values());
                    }
                }
            }
            GlobalBasisMatrix_->fillComplete(AssembledBasisMapUnique_,rangeMap); //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout)); GlobalBasisMatrix_->describe(*fancy,VERB_EXTREME);
        }

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::clearCoarseSpace()
    {
        ConstXMapPtrVec emptyVec1;
        UnassembledBasesMaps_.swap(emptyVec1);

        ConstXMapPtrVec emptyVec2;
        UnassembledBasesMapsUnique_.swap(emptyVec2);

        ConstXMultiVectorPtrVec emptyVec3;
        UnassembledSubspaceBases_.swap(emptyVec3);

        AssembledBasisMap_.reset();
        AssembledBasisMapUnique_.reset();
        AssembledBasis_.reset();

        UNVec emptyVec4;
        LocalSubspacesSizes_.swap(emptyVec4);

        GlobalBasisMatrix_.reset();

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int CoarseSpace<SC,LO,GO,NO>::zeroOutBasisVectors(ConstLOVecView zeros)
    {
        FROSCH_ASSERT(!AssembledBasis_.is_null(),"FROSch::CoarseSpace: AssembledBasis_.is_null().");
        for (UN j=0; j<zeros.size(); j++) {
            AssembledBasis_->getVectorNonConst(zeros[j])->scale(ScalarTraits<SC>::zero());
        }
        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasUnassembledMaps() const
    {
        return UnassembledBasesMaps_.size()>0;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasBasisMap() const
    {
        return !AssembledBasisMap_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::ConstXMapPtr CoarseSpace<SC,LO,GO,NO>::getBasisMap() const
    {
        FROSCH_ASSERT(!AssembledBasisMap_.is_null(),"FROSch::CoarseSpace: AssembledBasisMap_.is_null().");
        return AssembledBasisMap_;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasBasisMapUnique() const
    {
        return !AssembledBasisMapUnique_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::ConstXMapPtr CoarseSpace<SC,LO,GO,NO>::getBasisMapUnique() const
    {
        FROSCH_ASSERT(!AssembledBasisMapUnique_.is_null(),"FROSch::CoarseSpace: AssembledBasisMapUnique_.is_null().");
        return AssembledBasisMapUnique_;
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasAssembledBasis() const
    {
        return !AssembledBasis_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::ConstXMultiVectorPtr CoarseSpace<SC,LO,GO,NO>::getAssembledBasis() const
    {
        FROSCH_ASSERT(!AssembledBasis_.is_null(),"FROSch::CoarseSpace: AssembledBasis_.is_null().");
        return AssembledBasis_;
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::ConstUNVecView CoarseSpace<SC,LO,GO,NO>::getLocalSubspaceSizes() const
    {
        return LocalSubspacesSizes_();
    }

    template <class SC,class LO,class GO,class NO>
    bool CoarseSpace<SC,LO,GO,NO>::hasGlobalBasisMatrix() const
    {
        return !GlobalBasisMatrix_.is_null();
    }

    template <class SC,class LO,class GO,class NO>
    typename CoarseSpace<SC,LO,GO,NO>::XMatrixPtr CoarseSpace<SC,LO,GO,NO>::getGlobalBasisMatrix() const
    {
        FROSCH_ASSERT(!GlobalBasisMatrix_.is_null(),"FROSch::CoarseSpace: GlobalBasisMatrix_.is_null().");
        return GlobalBasisMatrix_;
    }
}

#endif
