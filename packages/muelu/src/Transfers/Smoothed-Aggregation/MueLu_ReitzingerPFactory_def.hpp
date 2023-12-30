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
#ifndef MUELU_REITZINGERPFACTORY_DEF_HPP
#define MUELU_REITZINGERPFACTORY_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportUtils.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_IO.hpp>

#include "MueLu_ReitzingerPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: enable");
  SET_VALID_ENTRY("repartition: use subcommunicators");
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("tentative: constant column sums");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  validParamList->set<RCP<const FactoryBase> >("D0", Teuchos::null, "Generating factory of the matrix D0");
  validParamList->set<RCP<const FactoryBase> >("NodeAggMatrix", Teuchos::null, "Generating factory of the matrix NodeAggMatrix");
  validParamList->set<RCP<const FactoryBase> >("Pnodal", Teuchos::null, "Generating factory of the matrix P");
  validParamList->set<RCP<const FactoryBase> >("NodeImporter", Teuchos::null, "Generating factory of the matrix NodeImporter");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "A");
  Input(fineLevel, "D0");
  Input(fineLevel, "NodeAggMatrix");
  Input(coarseLevel, "NodeAggMatrix");
  Input(coarseLevel, "Pnodal");
  //    Input(coarseLevel, "NodeImporter");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);
  using Teuchos::arcp_const_cast;
  using MT                    = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using XMM                   = Xpetra::MatrixMatrix<SC, LO, GO, NO>;
  Teuchos::FancyOStream& out0 = GetBlackHole();
  const ParameterList& pL     = GetParameterList();

  bool update_communicators = pL.get<bool>("repartition: enable") && pL.get<bool>("repartition: use subcommunicators");

  // If these are set correctly we assume that the nodal P contains only ones
  bool nodal_p_is_all_ones = !pL.get<bool>("tentative: constant column sums") && !pL.get<bool>("tentative: calculate qr");

  RCP<Matrix> EdgeMatrix = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<Matrix> D0         = Get<RCP<Matrix> >(fineLevel, "D0");
  RCP<Matrix> NodeMatrix = Get<RCP<Matrix> >(fineLevel, "NodeAggMatrix");
  RCP<Matrix> Pn         = Get<RCP<Matrix> >(coarseLevel, "Pnodal");

  const GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();
  const LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // This needs to be an Operator because if NodeMatrix gets repartitioned away, we get an Operator on the level
  RCP<Operator> CoarseNodeMatrix = Get<RCP<Operator> >(coarseLevel, "NodeAggMatrix");
  int MyPID                      = EdgeMatrix.is_null() ? -1 : EdgeMatrix->getRowMap()->getComm()->getRank();

  // Matrix matrix params
  RCP<ParameterList> mm_params = rcp(new ParameterList);
  ;
  if (pL.isSublist("matrixmatrix: kernel params"))
    mm_params->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

  // Normalize P
  if (!nodal_p_is_all_ones) {
    // The parameters told us the nodal P isn't all ones, so we make a copy that is.
    GetOStream(Runtime0) << "ReitzingerPFactory::BuildP(): Assuming Pn is not normalized" << std::endl;
    RCP<Matrix> Pn_old = Pn;

    Pn = Xpetra::MatrixFactory<SC, LO, GO, NO>::Build(Pn->getCrsGraph());
    Pn->setAllToScalar(Teuchos::ScalarTraits<SC>::one());
    Pn->fillComplete(Pn->getDomainMap(), Pn->getRangeMap());
  } else {
    // The parameters claim P is all ones.
    GetOStream(Runtime0) << "ReitzingerPFactory::BuildP(): Assuming Pn is normalized" << std::endl;
  }

  // TODO: We need to make sure Pn isn't normalized.  Right now this has to be done explicitly by the user

  // TODO: We need to look through and see which of these really need importers and which ones don't

  /* Generate the D0 * Pn matrix and its transpose */
  RCP<Matrix> D0_Pn, PnT_D0T, D0_Pn_nonghosted;
  Teuchos::Array<int> D0_Pn_col_pids;
  {
    RCP<Matrix> dummy;
    SubFactoryMonitor m2(*this, "Generate D0*Pn", coarseLevel);
    D0_Pn = XMM::Multiply(*D0, false, *Pn, false, dummy, out0, true, true, "D0*Pn", mm_params);

    // We don't want this guy getting accidently used later
    if (!mm_params.is_null()) mm_params->remove("importer", false);

    // Save this so we don't need to do the multiplication again later
    D0_Pn_nonghosted = D0_Pn;

    // Get owning PID information on columns for tie-breaking
    if (!D0_Pn->getCrsGraph()->getImporter().is_null()) {
      Xpetra::ImportUtils<LO, GO, NO> utils;
      utils.getPids(*D0_Pn->getCrsGraph()->getImporter(), D0_Pn_col_pids, false);
    } else {
      D0_Pn_col_pids.resize(D0_Pn->getCrsGraph()->getColMap()->getLocalNumElements(), MyPID);
    }
  }

  {
    // Get the transpose
    SubFactoryMonitor m2(*this, "Transpose D0*Pn", coarseLevel);
    PnT_D0T = Utilities::Transpose(*D0_Pn, true);
  }

  // We really need a ghosted version of D0_Pn here.
  // The reason is that if there's only one fine edge between two coarse nodes, somebody is going
  // to own the associated coarse edge.  The sum/sign rule doesn't guarantee the fine owner is the coarse owner.
  // So you can wind up with a situation that only guy who *can* register the coarse edge isn't the sum/sign
  // owner.  Adding more ghosting fixes that.
  if (!PnT_D0T->getCrsGraph()->getImporter().is_null()) {
    RCP<const Import> Importer     = PnT_D0T->getCrsGraph()->getImporter();
    RCP<const CrsMatrix> D0_Pn_crs = rcp_dynamic_cast<const CrsMatrixWrap>(D0_Pn)->getCrsMatrix();
    RCP<Matrix> D0_Pn_new          = rcp(new CrsMatrixWrap(CrsMatrixFactory::Build(D0_Pn_crs, *Importer, D0_Pn->getDomainMap(), Importer->getTargetMap())));
    D0_Pn                          = D0_Pn_new;
    // Get owning PID information on columns for tie-breaking
    if (!D0_Pn->getCrsGraph()->getImporter().is_null()) {
      Xpetra::ImportUtils<LO, GO, NO> utils;
      utils.getPids(*D0_Pn->getCrsGraph()->getImporter(), D0_Pn_col_pids, false);
    } else {
      D0_Pn_col_pids.resize(D0_Pn->getCrsGraph()->getColMap()->getLocalNumElements(), MyPID);
    }
  }

  // FIXME: This is using deprecated interfaces
  ArrayView<const LO> colind_E, colind_N;
  ArrayView<const SC> values_E, values_N;

  size_t Ne = EdgeMatrix->getLocalNumRows();
  size_t Nn = NodeMatrix->getLocalNumRows();

  // Upper bound on local number of coarse edges
  size_t max_edges = (NodeMatrix->getLocalNumEntries() + Nn + 1) / 2;
  ArrayRCP<size_t> D0_rowptr(Ne + 1);
  ArrayRCP<LO> D0_colind(max_edges);
  ArrayRCP<SC> D0_values(max_edges);
  D0_rowptr[0] = 0;

  LO current = 0;
  LO Nnc     = PnT_D0T->getRowMap()->getLocalNumElements();

  // Get the node maps for D0_coarse
  RCP<const Map> ownedCoarseNodeMap           = Pn->getDomainMap();
  RCP<const Map> ownedPlusSharedCoarseNodeMap = D0_Pn->getCrsGraph()->getColMap();

  for (LO i = 0; i < (LO)Nnc; i++) {
    LO local_column_i = ownedPlusSharedCoarseNodeMap->getLocalElement(PnT_D0T->getRowMap()->getGlobalElement(i));

    // FIXME: We don't really want an std::map here.  This is just a first cut implementation
    using value_type = bool;
    std::map<LO, value_type> ce_map;

    // FIXME: This is using deprecated interfaces
    PnT_D0T->getLocalRowView(i, colind_E, values_E);

    for (LO j = 0; j < (LO)colind_E.size(); j++) {
      // NOTE: Edges between procs will be via handled via the a version
      // of ML's odd/even rule
      // For this to function correctly, we make two assumptions:
      //  (a) The processor that owns a fine edge owns at least one of the attached nodes.
      //  (b) Aggregation is uncoupled.

      // TODO: Add some debug code to check the assumptions

      // Check to see if we own this edge and continue if we don't
      GO edge_gid = PnT_D0T->getColMap()->getGlobalElement(colind_E[j]);
      LO j_row    = D0_Pn->getRowMap()->getLocalElement(edge_gid);
      int pid0, pid1;
      D0_Pn->getLocalRowView(j_row, colind_N, values_N);

      // Skip incomplete rows
      if (colind_N.size() != 2) continue;

      pid0 = D0_Pn_col_pids[colind_N[0]];
      pid1 = D0_Pn_col_pids[colind_N[1]];
      //        printf("[%d] Row %d considering edge (%d)%d -> (%d)%d\n",MyPID,global_i,colind_N[0],D0_Pn->getColMap()->getGlobalElement(colind_N[0]),colind_N[1],D0_Pn->getColMap()->getGlobalElement(colind_N[1]));

      // Check to see who owns these nodes
      // If the sum of owning procs is odd, the lower ranked proc gets it

      bool zero_matches     = pid0 == MyPID;
      bool one_matches      = pid1 == MyPID;
      bool keep_shared_edge = false, own_both_nodes = false;
      if (zero_matches && one_matches) {
        own_both_nodes = true;
      } else {
        bool sum_is_even  = (pid0 + pid1) % 2 == 0;
        bool i_am_smaller = MyPID == std::min(pid0, pid1);
        if (sum_is_even && i_am_smaller) keep_shared_edge = true;
        if (!sum_is_even && !i_am_smaller) keep_shared_edge = true;
      }
      //        printf("[%d] - matches %d/%d keep_shared = %d own_both = %d\n",MyPID,(int)zero_matches,(int)one_matches,(int)keep_shared_edge,(int)own_both_nodes);
      if (!keep_shared_edge && !own_both_nodes) continue;

      // We're doing this in GID space, but only because it allows us to explain
      // the edge orientation as "always goes from lower GID to higher GID".  This could
      // be done entirely in LIDs, but then the ordering is a little more confusing.
      // This could be done in local indices later if we need the extra performance.
      for (LO k = 0; k < (LO)colind_N.size(); k++) {
        LO my_colind = colind_N[k];
        if (my_colind != LO_INVALID && ((keep_shared_edge && my_colind != local_column_i) || (own_both_nodes && my_colind > local_column_i))) {
          ce_map.emplace(std::make_pair(my_colind, true));
        }
      }  // end for k < colind_N.size()
    }    // end for j < colind_E.size()

    // std::map is sorted, so we'll just iterate through this
    for (auto iter = ce_map.begin(); iter != ce_map.end(); iter++) {
      LO col = iter->first;
      if (col == local_column_i) {
        continue;
      }

      // NOTE: "i" here might not be a valid local column id, so we read it from the map
      D0_colind[current] = local_column_i;
      D0_values[current] = -1;
      current++;
      D0_colind[current] = col;
      D0_values[current] = 1;
      current++;
      D0_rowptr[current / 2] = current;
    }

  }  // end for i < Nn

  LO num_coarse_edges = current / 2;
  D0_rowptr.resize(num_coarse_edges + 1);
  D0_colind.resize(current);
  D0_values.resize(current);

  // We're assuming that if the coarse NodeMatrix has no nodes on a rank, the coarse edge guy won't either.
  // We check that here.
  TEUCHOS_TEST_FOR_EXCEPTION((num_coarse_edges > 0 && CoarseNodeMatrix.is_null()) ||
                                 (num_coarse_edges == 0 && !CoarseNodeMatrix.is_null()),
                             Exceptions::RuntimeError, "MueLu::ReitzingerPFactory: Mismatched num_coarse_edges and NodeMatrix repartition.");

  // Count the total number of edges
  // NOTE: Since we solve the ownership issue above, this should do what we want
  RCP<const Map> ownedCoarseEdgeMap = Xpetra::MapFactory<LO, GO, NO>::Build(EdgeMatrix->getRowMap()->lib(), GO_INVALID, num_coarse_edges, EdgeMatrix->getRowMap()->getIndexBase(), EdgeMatrix->getRowMap()->getComm());

  // Create the coarse D0
  RCP<CrsMatrix> D0_coarse;
  {
    SubFactoryMonitor m2(*this, "Build D0", coarseLevel);
    // FIXME: We can be smarter with memory here
    // TODO: Is there a smarter way to get this importer?
    D0_coarse = CrsMatrixFactory::Build(ownedCoarseEdgeMap, ownedPlusSharedCoarseNodeMap, 0);
    TEUCHOS_TEST_FOR_EXCEPTION(D0_coarse.is_null(), Exceptions::RuntimeError, "MueLu::ReitzingerPFactory: CrsMatrixFatory failed.");

    // FIXME: Deprecated code
    ArrayRCP<size_t> ia;
    ArrayRCP<LO> ja;
    ArrayRCP<SC> val;
    D0_coarse->allocateAllValues(current, ia, ja, val);
    std::copy(D0_rowptr.begin(), D0_rowptr.end(), ia.begin());
    std::copy(D0_colind.begin(), D0_colind.end(), ja.begin());
    std::copy(D0_values.begin(), D0_values.end(), val.begin());
    D0_coarse->setAllValues(ia, ja, val);

#if 0
      {
        char fname[80];
        printf("[%d] D0: ia.size() = %d ja.size() = %d\n",MyPID,(int)ia.size(),(int)ja.size());
        printf("[%d] D0: ia  :",MyPID);
        for(int i=0; i<(int)ia.size(); i++)
          printf("%d ",(int)ia[i]);
        printf("\n[%d] D0: global ja  :",MyPID);
        for(int i=0; i<(int)ja.size(); i++)
          printf("%d ",(int)ownedPlusSharedCoarseNodeMap->getGlobalElement(ja[i]));
        printf("\n[%d] D0: local ja  :",MyPID);
        for(int i=0; i<(int)ja.size(); i++)
          printf("%d ",(int)ja[i]);
        printf("\n");

        sprintf(fname,"D0_global_ja_%d_%d.dat",MyPID,fineLevel.GetLevelID());
        FILE * f = fopen(fname,"w");
        for(int i=0; i<(int)ja.size(); i++)
          fprintf(f,"%d ",(int)ownedPlusSharedCoarseNodeMap->getGlobalElement(ja[i]));
        fclose(f);

        sprintf(fname,"D0_local_ja_%d_%d.dat",MyPID,fineLevel.GetLevelID());
        f = fopen(fname,"w");
        for(int i=0; i<(int)ja.size(); i++)
          fprintf(f,"%d ",(int)ja[i]);
        fclose(f);
        
      }
#endif
    D0_coarse->expertStaticFillComplete(ownedCoarseNodeMap, ownedCoarseEdgeMap);
  }
  RCP<Matrix> D0_coarse_m         = rcp(new CrsMatrixWrap(D0_coarse));
  RCP<Teuchos::FancyOStream> fout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  // Create the Pe matrix, but with the extra entries.  From ML's notes:
  /* The general idea is that the matrix                              */
  /*                        T_h P_n T_H^*                             */
  /* is almost Pe. If we make sure that P_n contains 1's and -1's, the*/
  /* matrix triple product will yield a matrix with +/- 1 and +/- 2's.*/
  /* If we remove all the 1's and divide the 2's by 2. we arrive at Pe*/
  RCP<Matrix> Pe;
  {
    SubFactoryMonitor m2(*this, "Generate Pe (pre-fix)", coarseLevel);

    RCP<Matrix> dummy;
    RCP<Matrix> Pn_D0cT = XMM::Multiply(*Pn, false, *D0_coarse_m, true, dummy, out0, true, true, "Pn*D0c'", mm_params);

    // We don't want this guy getting accidently used later
    if (!mm_params.is_null()) mm_params->remove("importer", false);

    Pe = XMM::Multiply(*D0, false, *Pn_D0cT, false, dummy, out0, true, true, "D0*(Pn*D0c')", mm_params);

    // TODO: Something like this *might* work.  But this specifically, doesn't
    // Pe = XMM::Multiply(*D0_Pn_nonghosted,false,*D0_coarse_m,true,dummy,out0,true,true,"(D0*Pn)*D0c'",mm_params);
  }

  /* Weed out the +/- entries, shrinking the matrix as we go */
  {
    SubFactoryMonitor m2(*this, "Generate Pe (post-fix)", coarseLevel);
    Pe->resumeFill();
    SC one     = Teuchos::ScalarTraits<SC>::one();
    MT two     = 2 * Teuchos::ScalarTraits<MT>::one();
    SC zero    = Teuchos::ScalarTraits<SC>::zero();
    SC neg_one = -one;

    RCP<const CrsMatrix> Pe_crs = rcp_dynamic_cast<const CrsMatrixWrap>(Pe)->getCrsMatrix();
    TEUCHOS_TEST_FOR_EXCEPTION(Pe_crs.is_null(), Exceptions::RuntimeError, "MueLu::ReitzingerPFactory: Pe is not a crs matrix.");
    ArrayRCP<const size_t> rowptr_const;
    ArrayRCP<const LO> colind_const;
    ArrayRCP<const SC> values_const;
    Pe_crs->getAllValues(rowptr_const, colind_const, values_const);
    ArrayRCP<size_t> rowptr = arcp_const_cast<size_t>(rowptr_const);
    ArrayRCP<LO> colind     = arcp_const_cast<LO>(colind_const);
    ArrayRCP<SC> values     = arcp_const_cast<SC>(values_const);
    LO ct                   = 0;
    LO lower                = rowptr[0];
    for (LO i = 0; i < (LO)Ne; i++) {
      for (size_t j = lower; j < rowptr[i + 1]; j++) {
        if (values[j] == one || values[j] == neg_one || values[j] == zero) {
          // drop this guy
        } else {
          colind[ct] = colind[j];
          values[ct] = values[j] / two;
          ct++;
        }
      }
      lower         = rowptr[i + 1];
      rowptr[i + 1] = ct;
    }
    rowptr[Ne] = ct;
    colind.resize(ct);
    values.resize(ct);
    rcp_const_cast<CrsMatrix>(Pe_crs)->setAllValues(rowptr, colind, values);

    Pe->fillComplete(Pe->getDomainMap(), Pe->getRangeMap());
  }

  /* Check commuting property */
  CheckCommutingProperty(*Pe, *D0_coarse_m, *D0, *Pn);

  /*  If we're repartitioning here, we need to cut down the communicators */
  // NOTE: We need to do this *after* checking the commuting property, since
  // that's going to need to fineLevel's communicators, not the repartitioned ones
  if (update_communicators) {
    // NOTE: We can only do D0 here.  We have to do Ke_coarse=(Re Ke_fine Pe) in RebalanceAcFactory
    RCP<const Teuchos::Comm<int> > newComm;
    if (!CoarseNodeMatrix.is_null()) newComm = CoarseNodeMatrix->getDomainMap()->getComm();
    RCP<const Map> newMap = Xpetra::MapFactory<LO, GO, NO>::copyMapWithNewComm(D0_coarse_m->getRowMap(), newComm);
    D0_coarse_m->removeEmptyProcessesInPlace(newMap);

    // The "in place" still leaves a dummy matrix here.  That needs to go
    if (newMap.is_null()) D0_coarse_m = Teuchos::null;

    Set(coarseLevel, "InPlaceMap", newMap);
  }

  /* Set output on the level */
  Set(coarseLevel, "P", Pe);
  Set(coarseLevel, "Ptent", Pe);

  Set(coarseLevel, "D0", D0_coarse_m);
  coarseLevel.Set("D0", D0_coarse_m, NoFactory::get());
  coarseLevel.AddKeepFlag("D0", NoFactory::get(), MueLu::Final);
  coarseLevel.RemoveKeepFlag("D0", NoFactory::get(), MueLu::UserData);

#if 0
  {
    int numProcs = Pe->getRowMap()->getComm()->getSize();
    char fname[80];

    sprintf(fname,"Pe_%d_%d.mat",numProcs,fineLevel.GetLevelID());  Xpetra::IO<SC,LO,GO,NO>::Write(fname,*Pe);
    sprintf(fname,"Pn_%d_%d.mat",numProcs,fineLevel.GetLevelID());  Xpetra::IO<SC,LO,GO,NO>::Write(fname,*Pn);
    sprintf(fname,"D0c_%d_%d.mat",numProcs,fineLevel.GetLevelID());  Xpetra::IO<SC,LO,GO,NO>::Write(fname,*D0_coarse_m);
    sprintf(fname,"D0f_%d_%d.mat",numProcs,fineLevel.GetLevelID());  Xpetra::IO<SC,LO,GO,NO>::Write(fname,*D0);
  }
#endif

}  // end Build

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CheckCommutingProperty(const Matrix& Pe, const Matrix& D0_c, const Matrix& D0_f, const Matrix& Pn) const {
  if (IsPrint(Statistics0)) {
    using XMM = Xpetra::MatrixMatrix<SC, LO, GO, NO>;
    using MT  = typename Teuchos::ScalarTraits<SC>::magnitudeType;
    SC one    = Teuchos::ScalarTraits<SC>::one();
    SC zero   = Teuchos::ScalarTraits<SC>::zero();

    RCP<Matrix> dummy;
    Teuchos::FancyOStream& out0 = GetBlackHole();
    RCP<Matrix> left            = XMM::Multiply(Pe, false, D0_c, false, dummy, out0);
    RCP<Matrix> right           = XMM::Multiply(D0_f, false, Pn, false, dummy, out0);

    // We need a non-FC matrix for the add, sadly
    RCP<CrsMatrix> sum_c  = CrsMatrixFactory::Build(left->getRowMap(), left->getLocalMaxNumRowEntries() + right->getLocalMaxNumRowEntries());
    RCP<Matrix> summation = rcp(new CrsMatrixWrap(sum_c));
    XMM::TwoMatrixAdd(*left, false, one, *summation, zero);
    XMM::TwoMatrixAdd(*right, false, -one, *summation, one);

    MT norm = summation->getFrobeniusNorm();
    GetOStream(Statistics0) << "CheckCommutingProperty: ||Pe D0_c - D0_f Pn || = " << norm << std::endl;
  }

}  // end CheckCommutingProperty

}  // namespace MueLu

#define MUELU_REITZINGERPFACTORY_SHORT
#endif  // MUELU_REITZINGERPFACTORY_DEF_HPP
