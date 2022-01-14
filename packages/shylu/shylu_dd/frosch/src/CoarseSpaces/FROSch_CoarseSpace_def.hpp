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
            FROSCH_ASSERT(subspaceBasis->getNumVectors()==subspaceBasisMap->getNodeNumElements(),"FROSch::CoarseSpace: subspaceBasis->getNumVectors()!=subspaceBasisMap->getNodeNumElements()");
        } else {
            FROSCH_ASSERT(subspaceBasisMap->getNodeNumElements()==0,"FROSch::CoarseSpace: subspaceBasisMap->getNodeNumElements()!=0");
        }

        UnassembledBasesMaps_.push_back(subspaceBasisMap);
        UnassembledBasesMapsUnique_.push_back(subspaceBasisMapUnique);
        UnassembledSubspaceBases_.push_back(subspaceBasis);
        Offsets_.push_back(offset);
        LocalSubspacesSizes_.push_back(subspaceBasisMap->getNodeNumElements());

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

                AssembledBasis_ = MultiVectorFactory<SC,LO,GO,NO >::Build(serialMap,AssembledBasisMap_->getNodeNumElements());
                #if defined(HAVE_XPETRA_KOKKOS_REFACTOR) && defined(HAVE_XPETRA_TPETRA)
                if (AssembledBasis_->getMap()->lib() == UseTpetra) {
                    UN itmp = 0;
                    for (UN i=0; i<UnassembledSubspaceBases_.size(); i++) {
                        if (!UnassembledSubspaceBases_[i].is_null()) {
                            const UN Offset_i = Offsets_[i];
                            const UN NumVectors_i = UnassembledSubspaceBases_[i]->getNumVectors();
                            const UN LocalLength_i = UnassembledSubspaceBases_[i]->getLocalLength();

                            FROSCH_ASSERT(NumVectors_i+itmp <= AssembledBasis_->getNumVectors(),"FROSch::CoarseSpace: NumVectors_i+itmp <= AssembledBasis_->getNumVectors()");
                            FROSCH_ASSERT(LocalLength_i+Offsets_[i] <= AssembledBasis_->getLocalLength(),"FROSch::CoarseSpace: LocalLength_i+Offsets_[i] <= AssembledBasis_");

                            for (UN j=0; j < NumVectors_i; j++) {
                                auto unassembledSubspaceBasesData = UnassembledSubspaceBases_[i]->getData(j).getRawPtr();
                                auto   assembledSubspaceBasesData = AssembledBasis_->getDataNonConst(itmp+j);

                                using execution_space = typename XMap::local_map_type::execution_space;
                                Kokkos::RangePolicy<execution_space> policy (0, LocalLength_i);
                                Kokkos::parallel_for(
                                    "FROSch_CoarseSpace::assembleCoarseSpace", policy,
                                    KOKKOS_LAMBDA(const UN k) {
                                        assembledSubspaceBasesData[k+Offset_i] = unassembledSubspaceBasesData[k];
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

#if defined(HAVE_XPETRA_KOKKOS_REFACTOR) && defined(HAVE_XPETRA_TPETRA)
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
            UN numLocalRows = rowMap->getNodeNumElements();
            rowptr_type Rowptr ("Rowptr", numLocalRows+1);
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

            // make it into offsets
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
              (1+numLocalRows, Rowptr);

            // fill into the local matrix
            indices_type Indices ("Indices", nnz);
            values_type  Values  ("Values",  nnz);
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
                GlobalBasisMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(rowMap,AssembledBasisMap_->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!
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
                GlobalBasisMatrix_ = MatrixFactory<SC,LO,GO,NO>::Build(rowMap,AssembledBasisMap_,AssembledBasisMap_->getNodeNumElements()); // Nonzeroes abhängig von dim/dofs!!!
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
