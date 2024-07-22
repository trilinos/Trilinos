// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _FROSCH_EXTRACTSUBMATRICES_DEF_HPP
#define _FROSCH_EXTRACTSUBMATRICES_DEF_HPP

#include <FROSch_ExtractSubmatrices_decl.hpp>

namespace FROSch {

    using namespace std;
    using namespace Teuchos;
    using namespace Xpetra;

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map)
    {
        FROSCH_DETAILTIMER_START(extractLocalSubdomainMatrixTime,"ExtractLocalSubdomainMatrix");
        RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(map,globalMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,ADD);
        //cout << *subdomainMatrix << endl;
        RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
        RCP<Map<LO,GO,NO> > localSubdomainMap = MapFactory<LO,GO,NO>::Build(map->lib(),map->getLocalNumElements(),0,SerialComm);
        RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap,globalMatrix->getGlobalMaxNumRowEntries());

        for (unsigned i=0; i<localSubdomainMap->getLocalNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<GO> indicesLocal;
                Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    GO localIndex = map->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesLocal.push_back(localIndex);
                        valuesLocal.push_back(values[j]);
                    }
                }
                localSubdomainMatrix->insertGlobalValues(i,indicesLocal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();
        return localSubdomainMatrix.getConst();
    }

    // this version just read indices without building submatrices, which is done in extractLocalSubdomainMatrix_Symbolic
    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Symbolic(RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix,       // input  : globalMatrix, re-distributed with map
                                              RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix)  // output : local submatrix
    {
        FROSCH_DETAILTIMER_START(extractLocalSubdomainMatrixTime_symbolic, "ExtractLocalSubdomainMatrix_Symbolic");
        auto subdomainMap = subdomainMatrix->getRowMap();

        const SC zero = ScalarTraits<SC>::zero();
        for (unsigned i=0; i<subdomainMap->getLocalNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(subdomainMap->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<GO> indicesLocal;
                Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    LO localIndex = subdomainMap->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesLocal.push_back(localIndex);
                        valuesLocal.push_back(zero);
                    }
                }
                localSubdomainMatrix->insertGlobalValues(i,indicesLocal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();
        return;
    }

    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Compute(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > subdomainMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > localSubdomainMatrix)
    {
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(), subdomainMatrix->getRowMap());
        ExtractLocalSubdomainMatrix_Compute(scatter, globalMatrix, subdomainMatrix, localSubdomainMatrix);
    }

    template <class SC,class LO,class GO,class NO>
    void ExtractLocalSubdomainMatrix_Compute(RCP<      Import<LO,GO,NO> >    scatter,
                                             RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > subdomainMatrix,
                                             RCP<      Matrix<SC,LO,GO,NO> > localSubdomainMatrix)
    {
        FROSCH_DETAILTIMER_START(extractLocalSubdomainMatrixTime_compute, "ExtractLocalSubdomainMatrix_Compute");
        const SC zero = ScalarTraits<SC>::zero();
        auto subdomainRowMap = subdomainMatrix->getRowMap();

        subdomainMatrix->setAllToScalar(zero);
        subdomainMatrix->resumeFill();
        subdomainMatrix->doImport(*globalMatrix, *scatter, ADD);

// "fillComplete" is quite expensive, and it seem to be cheaper to replace values each row at a time
#if 0 //defined(HAVE_XPETRA_TPETRA)
        if (globalMatrix->getRowMap()->lib() == UseTpetra) 
        {
            // NOTE: this fillComplete is expensive on GPUs
            subdomainMatrix->fillComplete();
            auto devSubdomainMap         = subdomainRowMap->getLocalMap();
            auto devSubdomainMatrix      = subdomainMatrix->getLocalMatrixDevice();
            auto devLocalSubdomainMatrix = localSubdomainMatrix->getLocalMatrixDevice();

            using UN = unsigned;
            using execution_space = typename Map<LO,GO,NO>::local_map_type::execution_space;
            UN numRows = subdomainRowMap->getLocalNumElements();
            Kokkos::RangePolicy<execution_space> policy_row (0, numRows);
            Kokkos::parallel_for(
                "FROSch :: ExtractLocalSubdomainMatrix_Compute", policy_row,
                KOKKOS_LAMBDA(const UN i) {
                    for (UN j=devSubdomainMatrix.graph.row_map(i); j<devSubdomainMatrix.graph.row_map(i+1); j++) {
                        GO globalIndex_j = devSubdomainMap.getGlobalElement(devSubdomainMatrix.graph.entries(j));
                        GO localIndex_j  = devSubdomainMap.getLocalElement(globalIndex_j);
                        if (localIndex_j>=0) {
                            // look for the same column index in localSubdomainMatrix
                            for (UN k=devLocalSubdomainMatrix.graph.row_map(i); k<devLocalSubdomainMatrix.graph.row_map(i+1); k++) {
                                if (devLocalSubdomainMatrix.graph.entries(k) == localIndex_j) {
                                    devLocalSubdomainMatrix.values[k] = devSubdomainMatrix.values[j];
                                    break;
                                }
                            }
                        }
                    }
                }
            );
            Kokkos::fence();
        } else
#endif
        {
            localSubdomainMatrix->resumeFill();

            size_t max_nnz = localSubdomainMatrix->getLocalMaxNumRowEntries();
            std::vector<LO> local_cols_vector (max_nnz);
            std::vector<SC> local_vals_vector (max_nnz);
            for (unsigned i=0; i<subdomainRowMap->getLocalNumElements(); i++) {
                ArrayView<const GO> global_indices;
                ArrayView<const SC> global_values;
                subdomainMatrix->getGlobalRowView(subdomainRowMap->getGlobalElement(i),global_indices,global_values);

                LO size = global_indices.size();
                if (size>0) {
                    // using "workspace" not to overwrite what's in localSubdomainMatrix
                    size_t new_nnz = 0;
                    ArrayView<LO> local_cols (local_cols_vector);
                    ArrayView<SC> local_vals (local_vals_vector);
                    for (LO j=0; j<size; j++) {
                        GO localIndex = subdomainRowMap->getLocalElement(global_indices[j]);
                        if (localIndex>=0) {
                            local_cols[new_nnz] = localIndex;
                            local_vals[new_nnz] = global_values[j];
                            new_nnz ++;
                        }
                    }
                    localSubdomainMatrix->replaceLocalValues(i, local_cols(0, new_nnz), local_vals(0, new_nnz));
                }
            }
            RCP<ParameterList> fillCompleteParams(new ParameterList);
            fillCompleteParams->set("No Nonlocal Changes", true);
            localSubdomainMatrix->fillComplete(fillCompleteParams);
        }

        return;
    }

    template <class SC,class LO,class GO,class NO>
    RCP<const Matrix<SC,LO,GO,NO> > ExtractLocalSubdomainMatrix(RCP<const Matrix<SC,LO,GO,NO> > globalMatrix,
                                                                RCP<const Map<LO,GO,NO> > map,
                                                                SC value)
    {
        FROSCH_DETAILTIMER_START(extractLocalSubdomainMatrixTime,"ExtractLocalSubdomainMatrix");
        RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,ADD);
        //cout << *subdomainMatrix << endl;
        RCP<const Comm<LO> > SerialComm = rcp(new MpiComm<LO>(MPI_COMM_SELF));
        RCP<Map<LO,GO,NO> > localSubdomainMap = MapFactory<LO,GO,NO>::Build(map->lib(),map->getLocalNumElements(),0,SerialComm);
        RCP<Matrix<SC,LO,GO,NO> > localSubdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(localSubdomainMap,globalMatrix->getGlobalMaxNumRowEntries());

        for (unsigned i=0; i<localSubdomainMap->getLocalNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<GO> indicesGlobal;
                Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    GO localIndex = map->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesGlobal.push_back(localIndex);
                        valuesLocal.push_back(value);
                    }
                }
                localSubdomainMatrix->insertGlobalValues(i,indicesGlobal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();
        return localSubdomainMatrix.getConst();
    }

    template <class SC,class LO,class GO,class NO>
    int UpdateLocalSubdomainMatrix(RCP<Matrix<SC,LO,GO,NO> > globalMatrix,
                                   RCP<Map<LO,GO,NO> > &map,
                                   RCP<Matrix<SC,LO,GO,NO> > &localSubdomainMatrix)
    {
        FROSCH_DETAILTIMER_START(updateLocalSubdomainMatrixTime,"UpdateLocalSubdomainMatrix");
        RCP<Matrix<SC,LO,GO,NO> > subdomainMatrix = MatrixFactory<SC,LO,GO,NO>::Build(map,2*globalMatrix->getGlobalMaxNumRowEntries());
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(globalMatrix->getRowMap(),map);
        subdomainMatrix->doImport(*globalMatrix,*scatter,ADD);

        localSubdomainMatrix->resumeFill();
        for (unsigned i=0; i<map->getLocalNumElements(); i++) {
            ArrayView<const GO> indices;
            ArrayView<const SC> values;
            subdomainMatrix->getGlobalRowView(map->getGlobalElement(i),indices,values);

            LO size = indices.size();
            if (size>0) {
                Array<LO> indicesLocal;
                Array<SC> valuesLocal;
                for (LO j=0; j<size; j++) {
                    GO localIndex = map->getLocalElement(indices[j]);
                    if (localIndex>=0) {
                        indicesLocal.push_back(localIndex);
                        valuesLocal.push_back(values[j]);
                    }
                }
                localSubdomainMatrix->replaceLocalValues(i,indicesLocal(),valuesLocal());
            }
        }
        localSubdomainMatrix->fillComplete();

        return 0;
    }

    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrices(RCP<const Matrix<SC,LO,GO,NO> > k,
                         ArrayView<GO> indI,
                         RCP<Matrix<SC,LO,GO,NO> > &kII,
                         RCP<Matrix<SC,LO,GO,NO> > &kIJ,
                         RCP<Matrix<SC,LO,GO,NO> > &kJI,
                         RCP<Matrix<SC,LO,GO,NO> > &kJJ)
    {
        FROSCH_DETAILTIMER_START(buildSubmatricesTime,"BuildSubmatrices");
        // We need four Maps
        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        RCP<Map<LO,GO,NO> > mapI = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),INVALID,indI(),0,k->getRowMap()->getComm());
        RCP<Map<LO,GO,NO> > mapILocal = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),INVALID,indI.size(),0,k->getRowMap()->getComm());

        Array<GO> indJ;
        for (unsigned i=0; i<k->getLocalNumRows(); i++) {
            if (mapI->getLocalElement(i)<0) {
                indJ.push_back(i);
            }
        }

        RCP<Map<LO,GO,NO> > mapJ = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),INVALID,indJ(),0,k->getRowMap()->getComm());
        RCP<Map<LO,GO,NO> > mapJLocal = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),INVALID,indJ.size(),0,k->getRowMap()->getComm());
        RCP<const Map<LO,GO,NO> > colMap = k->getColMap();
#if defined(HAVE_XPETRA_TPETRA)
        if (k->getRowMap()->lib() == UseTpetra) 
        {
            using crsmat_type  = typename Matrix<SC,LO,GO,NO>::local_matrix_type;
            using graph_type   = typename crsmat_type::StaticCrsGraphType;
            using rowptr_type  = typename graph_type::row_map_type::non_const_type;
            using indices_type = typename graph_type::entries_type::non_const_type;
            using values_type  = typename crsmat_type::values_type::non_const_type;
            using execution_space = typename Map<LO,GO,NO>::local_map_type::execution_space;

            using UN = unsigned;
            UN numRows = k->getLocalNumRows();
            UN numRowsI = indI.size();
            UN numRowsJ = indJ.size();

            rowptr_type RowptrII (Kokkos::ViewAllocateWithoutInitializing("RowptrII"), numRowsI+1);
            rowptr_type RowptrIJ (Kokkos::ViewAllocateWithoutInitializing("RowptrIJ"), numRowsI+1);
            rowptr_type RowptrJI (Kokkos::ViewAllocateWithoutInitializing("RowptrJI"), numRowsJ+1);
            rowptr_type RowptrJJ (Kokkos::ViewAllocateWithoutInitializing("RowptrJJ"), numRowsJ+1);

            // count nnz on each blocks
            Kokkos::deep_copy(RowptrII, 0);
            Kokkos::deep_copy(RowptrIJ, 0);
            Kokkos::deep_copy(RowptrJI, 0);
            Kokkos::deep_copy(RowptrJJ, 0);

            auto localMapI = mapI->getLocalMap();
            auto localMapJ = mapJ->getLocalMap();
            auto localColMap = colMap->getLocalMap();
            auto localK = k->getLocalMatrixDevice();
            Kokkos::RangePolicy<execution_space> policy_row (0, numRows);
            Kokkos::parallel_for(
                "FROSch_BuildSubmatrices::countNnz", policy_row,
                KOKKOS_LAMBDA(const UN i) {
                    LO tmp1=localMapI.getLocalElement(i);
                    LO tmp2=0;
                    if (tmp1>=0) {
                        for (UN j=localK.graph.row_map(i); j<localK.graph.row_map(i+1); j++) {
                            tmp2 = localMapI.getLocalElement(localColMap.getGlobalElement(localK.graph.entries(j)));
                            if (tmp2>=0) {
                                RowptrII(tmp1+1) ++;
                            } else {
                                RowptrIJ(tmp1+1) ++;
                            }
                        }
                    } else  {
                        tmp1=localMapJ.getLocalElement((GO) i);
                        for (UN j=localK.graph.row_map(i); j<localK.graph.row_map(i+1); j++) {
                            tmp2 = localMapI.getLocalElement(localColMap.getGlobalElement(localK.graph.entries(j)));
                            if (tmp2>=0) {
                                RowptrJI(tmp1+1) ++;
                            } else {
                                RowptrJJ(tmp1+1) ++;
                            }
                        }
                    }
                }
            );
            Kokkos::fence();

            // cout nnz
            UN nnzII = 0;
            UN nnzIJ = 0;
            UN nnzJI = 0;
            UN nnzJJ = 0;
            Kokkos::parallel_reduce("FROSch_CoarseSpace::BuildSubmatrices:nnzII", 1+numRowsI,
                KOKKOS_LAMBDA(const int &i, UN &lsum) { lsum += RowptrII[i]; },
                nnzII);
            Kokkos::parallel_reduce("FROSch_CoarseSpace::BuildSubmatrices:nnzIJ", 1+numRowsI,
                KOKKOS_LAMBDA(const int &i, UN &lsum) { lsum += RowptrIJ[i]; },
                nnzIJ);
            Kokkos::parallel_reduce("FROSch_CoarseSpace::BuildSubmatrices:nnzJI", 1+numRowsJ,
                KOKKOS_LAMBDA(const int &i, UN &lsum) { lsum += RowptrJI[i]; },
                nnzJI);
            Kokkos::parallel_reduce("FROSch_CoarseSpace::BuildSubmatrices:nnzJJ", 1+numRowsJ,
                KOKKOS_LAMBDA(const int &i, UN &lsum) { lsum += RowptrJJ[i]; },
                nnzJJ);

            // make it into offsets
#if KOKKOSKERNELS_VERSION >= 40199
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>
                (1+numRowsI, RowptrII);
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>
                (1+numRowsI, RowptrIJ);
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>
                (1+numRowsJ, RowptrJI);
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<execution_space>
                (1+numRowsJ, RowptrJJ);
#else
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
                (1+numRowsI, RowptrII);
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
                (1+numRowsI, RowptrIJ);
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
                (1+numRowsJ, RowptrJI);
            KokkosKernels::Impl::kk_inclusive_parallel_prefix_sum<rowptr_type, execution_space>
                (1+numRowsJ, RowptrJJ);
#endif

            // allocate kII block
            indices_type IndicesII (Kokkos::ViewAllocateWithoutInitializing("IndicesII"), nnzII);
            values_type  ValuesII  (Kokkos::ViewAllocateWithoutInitializing("ValuesII"),  nnzII);
            // allocate kIJ block
            indices_type IndicesIJ (Kokkos::ViewAllocateWithoutInitializing("IndicesIJ"), nnzIJ);
            values_type  ValuesIJ  (Kokkos::ViewAllocateWithoutInitializing("ValuesIJ"),  nnzIJ);
            // allocate kJI block
            indices_type IndicesJI (Kokkos::ViewAllocateWithoutInitializing("IndicesJI"), nnzJI);
            values_type  ValuesJI  (Kokkos::ViewAllocateWithoutInitializing("ValuesJI"),  nnzJI);
            // allocate kJJ block
            indices_type IndicesJJ (Kokkos::ViewAllocateWithoutInitializing("IndicesJJ"), nnzJJ);
            values_type  ValuesJJ  (Kokkos::ViewAllocateWithoutInitializing("ValuesJJ"),  nnzJJ);

            // fill in all the blocks
            Kokkos::parallel_for(
                "BuildSubmatrices::countNnz", policy_row,
                KOKKOS_LAMBDA(const UN i) {
                    LO tmp1=localMapI.getLocalElement((GO) i);
                    LO tmp2=0;
                    if (tmp1>=0) {
                        UN nnz_i = RowptrII[tmp1];
                        UN nnz_j = RowptrIJ[tmp1];
                        for (UN j=localK.graph.row_map(i); j<localK.graph.row_map(i+1); j++) {
                            LO colid = localColMap.getGlobalElement(localK.graph.entries(j));
                            tmp2 = localMapI.getLocalElement(colid);
                            if (tmp2>=0) {
                                IndicesII(nnz_i) = tmp2;
                                ValuesII(nnz_i)  = localK.values[j];
                                nnz_i++;
                            } else {
                                tmp2 = localMapJ.getLocalElement(colid);
                                IndicesIJ(nnz_j) = tmp2;
                                ValuesIJ(nnz_j)  = localK.values[j];
                                nnz_j++;
                            }
                        }
                    } else  {
                        tmp1=localMapJ.getLocalElement((GO) i);
                        UN nnz_i = RowptrJI[tmp1];
                        UN nnz_j = RowptrJJ[tmp1];
                        for (UN j=localK.graph.row_map(i); j<localK.graph.row_map(i+1); j++) {
                            LO colid = localColMap.getGlobalElement(localK.graph.entries(j));
                            tmp2 = localMapI.getLocalElement(colid);
                            if (tmp2>=0) {
                                IndicesJI(nnz_i) = tmp2;
                                ValuesJI(nnz_i)  = localK.values[j];
                                nnz_i++;
                            } else {
                                tmp2 = localMapJ.getLocalElement(colid);
                                IndicesJJ(nnz_j) = tmp2;
                                ValuesJJ(nnz_j)  = localK.values[j];
                                nnz_j++;
                            }
                        }
                    }
                }
            );
            Kokkos::fence();

            Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp (new Teuchos::ParameterList());
            params->set("sorted", false);

            // create kII
            graph_type crsgraphII (IndicesII, RowptrII);
            crsmat_type LocalII = crsmat_type ("CrsMatrix", numRowsI, ValuesII, crsgraphII);
            kII = MatrixFactory<SC,LO,GO,NO>::Build(LocalII, mapILocal, mapILocal, mapILocal, mapILocal,
                                                    params);
            // create kIJ
            graph_type crsgraphIJ (IndicesIJ, RowptrIJ);
            crsmat_type LocalIJ = crsmat_type ("CrsMatrix", numRowsI, ValuesIJ, crsgraphIJ);
            kIJ = MatrixFactory<SC,LO,GO,NO>::Build(LocalIJ, mapILocal, mapJLocal, mapJLocal, mapILocal,
                                                    params);
            // create kJI
            graph_type crsgraphJI (IndicesJI, RowptrJI);
            crsmat_type LocalJI = crsmat_type ("CrsMatrix", numRowsJ, ValuesJI, crsgraphJI);
            kJI = MatrixFactory<SC,LO,GO,NO>::Build(LocalJI, mapJLocal, mapILocal, mapILocal, mapJLocal,
                                                    params);
            // create kJJ
            graph_type crsgraphJJ (IndicesJJ, RowptrJJ);
            crsmat_type LocalJJ = crsmat_type ("CrsMatrix", numRowsJ, ValuesJJ, crsgraphJJ);
            kJJ = MatrixFactory<SC,LO,GO,NO>::Build(LocalJJ, mapJLocal, mapJLocal, mapJLocal, mapJLocal,
                                                    params);
        } else
#endif
        {
            kII = MatrixFactory<SC,LO,GO,NO>::Build(mapILocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
            kIJ = MatrixFactory<SC,LO,GO,NO>::Build(mapILocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indJ.size()));
            kJI = MatrixFactory<SC,LO,GO,NO>::Build(mapJLocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
            kJJ = MatrixFactory<SC,LO,GO,NO>::Build(mapJLocal,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indJ.size()));

            for (unsigned i=0; i<k->getLocalNumRows(); i++) {
                ArrayView<const LO> indices;
                ArrayView<const SC> values;

                k->getLocalRowView(i,indices,values);
                //cout << numEntries << endl;
                Array<GO> indicesI;
                Array<SC> valuesI;
                Array<GO> indicesJ;
                Array<SC> valuesJ;
                LO tmp1=mapI->getLocalElement(i);
                LO tmp2=0;
                if (tmp1>=0) {
                    for (LO j=0; j<indices.size(); j++) {
                        tmp2 = mapI->getLocalElement(colMap->getGlobalElement(indices[j]));
                        if (tmp2>=0) {
                            indicesI.push_back(tmp2);
                            valuesI.push_back(values[j]);
                        } else {
                            indicesJ.push_back(mapJ->getLocalElement(colMap->getGlobalElement(indices[j])));
                            valuesJ.push_back(values[j]);
                        }
                    }
                    kII->insertGlobalValues(tmp1,indicesI(),valuesI());
                    kIJ->insertGlobalValues(tmp1,indicesJ(),valuesJ());
                } else  {
                    tmp1=mapJ->getLocalElement((GO) i);
                    for (LO j=0; j<indices.size(); j++) {
                        tmp2 = mapI->getLocalElement(colMap->getGlobalElement(indices[j]));
                        if (tmp2>=0) {
                            indicesI.push_back(tmp2);
                            valuesI.push_back(values[j]);
                        } else {
                            indicesJ.push_back(mapJ->getLocalElement(colMap->getGlobalElement(indices[j])));
                            valuesJ.push_back(values[j]);
                        }
                    }
                    kJI->insertGlobalValues(tmp1,indicesI(),valuesI());
                    kJJ->insertGlobalValues(tmp1,indicesJ(),valuesJ());
                }
            }

            kII->fillComplete(mapILocal,mapILocal);
            kIJ->fillComplete(mapJLocal,mapILocal);
            kJI->fillComplete(mapILocal,mapJLocal);
            kJJ->fillComplete(mapJLocal,mapJLocal);
        }

        return 0;
    }


    template <class SC,class LO,class GO,class NO>
    int BuildSubmatrix(RCP<Matrix<SC,LO,GO,NO> > k,
                       ArrayView<GO> indI,
                       RCP<Matrix<SC,LO,GO,NO> > &kII)
    {
        FROSCH_DETAILTIMER_START(buildSubmatrixTime,"BuildSubmatrix");
        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));
        RCP<Map<LO,GO,NO> > mapI = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),INVALID,indI(),0,k->getRowMap()->getComm());

        kII = MatrixFactory<SC,LO,GO,NO>::Build(mapI,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        GO maxGID = mapI->getMaxAllGlobalIndex();
        GO minGID = mapI->getMinAllGlobalIndex();
        for (unsigned i=0; i<k->getLocalNumRows(); i++) {
            ArrayView<const LO> indices;
            ArrayView<const SC> values;

            k->getLocalRowView(i,indices,values);

            Array<GO> indicesI;
            Array<SC> valuesI;

            LO tmp1=mapI->getLocalElement(k->getRowMap()->getGlobalElement(i));
            GO tmp2=0;
            if (tmp1>=0) {
                for (LO j=0; j<indices.size(); j++) {
                    tmp2 = k->getColMap()->getGlobalElement(indices[j]);
                    if (minGID<=tmp2 && tmp2<=maxGID) {
                        indicesI.push_back(tmp2);
                        valuesI.push_back(values[j]);
                    }
                }
                kII->insertGlobalValues(mapI->getGlobalElement(tmp1),indicesI(),valuesI());
            }
        }
        kII->fillComplete(mapI,mapI);

        return 0;
    }

    template <class LO,class GO,class NO>
    int BuildSubgraph(RCP<CrsGraph<LO,GO,NO> > k,
                      ArrayView<GO> indI,
                      RCP<CrsGraph<LO,GO,NO> > &kII)
    {
        FROSCH_DETAILTIMER_START(buildSubgraphTime,"BuildSubgraph");
        const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid();
        //RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));
        RCP<Map<LO,GO,NO> > mapI = MapFactory<LO,GO,NO>::Build(k->getRowMap()->lib(),INVALID,indI(),0,k->getRowMap()->getComm());

        kII = CrsGraphFactory<LO,GO,NO>::Build(mapI,min((LO) k->getGlobalMaxNumRowEntries(),(LO) indI.size()));
        GO maxGID = mapI->getMaxAllGlobalIndex();
        GO minGID = mapI->getMinAllGlobalIndex();
        for (unsigned i=0; i<k->getLocalNumRows(); i++) {
            ArrayView<const LO> indices;

            k->getLocalRowView(i,indices);

            Array<GO> indicesI;

            LO tmp1=mapI->getLocalElement(k->getRowMap()->getGlobalElement(i));
            GO tmp2=0;
            if (tmp1>=0) {
                for (LO j=0; j<indices.size(); j++) {
                    tmp2 = k->getColMap()->getGlobalElement(indices[j]);
                    if (minGID<=tmp2 && tmp2<=maxGID) {
                        indicesI.push_back(tmp2);
                    }
                }
                kII->insertGlobalValues(mapI->getGlobalElement(tmp1),indicesI());
            }
        }
        kII->fillComplete(mapI,mapI);

        return 0;
    }
}

#endif
