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
#ifndef MUELU_CLASSICALPFACTORY_DEF_HPP
#define MUELU_CLASSICALPFACTORY_DEF_HPP

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportUtils.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_StridedMapFactory.hpp>

#include <Teuchos_OrdinalTraits.hpp>

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_ClassicalPFactory_decl.hpp"
#include "MueLu_ClassicalMapFactory.hpp"
#include "MueLu_AmalgamationInfo.hpp"
#include "MueLu_GraphBase.hpp"

//#define CMS_DEBUG
//#define CMS_DUMP

namespace {

template <class SC>
int Sign(SC val) {
  using STS                           = typename Teuchos::ScalarTraits<SC>;
  typename STS::magnitudeType MT_ZERO = Teuchos::ScalarTraits<typename STS::magnitudeType>::zero();
  if (STS::real(val) > MT_ZERO)
    return 1;
  else if (STS::real(val) < MT_ZERO)
    return -1;
  else
    return 0;
}

}  // namespace

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("aggregation: deterministic");
  SET_VALID_ENTRY("aggregation: coloring algorithm");
  SET_VALID_ENTRY("aggregation: classical scheme");

  // To know if we need BlockNumber
  SET_VALID_ENTRY("aggregation: drop scheme");
  {
    typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
    validParamList->getEntry("aggregation: classical scheme").setValidator(rcp(new validatorType(Teuchos::tuple<std::string>("direct", "ext+i", "classical modified"), "aggregation: classical scheme")));
  }

#undef SET_VALID_ENTRY
  validParamList->set<RCP<const FactoryBase> >("A", Teuchos::null, "Generating factory of the matrix A");
  //    validParamList->set< RCP<const FactoryBase> >("UnAmalgamationInfo", Teuchos::null, "Generating factory of UnAmalgamationInfo");
  validParamList->set<RCP<const FactoryBase> >("Graph", null, "Generating factory of the graph");
  validParamList->set<RCP<const FactoryBase> >("DofsPerNode", null, "Generating factory for variable \'DofsPerNode\', usually the same as for \'Graph\'");
  validParamList->set<RCP<const FactoryBase> >("CoarseMap", Teuchos::null, "Generating factory of the CoarseMap");
  validParamList->set<RCP<const FactoryBase> >("FC Splitting", Teuchos::null, "Generating factory of the FC Splitting");
  validParamList->set<RCP<const FactoryBase> >("BlockNumber", Teuchos::null, "Generating factory for Block Number");
  //    validParamList->set< RCP<const FactoryBase> >("Nullspace",      Teuchos::null, "Generating factory of the nullspace");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& /* coarseLevel */) const {
  Input(fineLevel, "A");
  Input(fineLevel, "Graph");
  Input(fineLevel, "DofsPerNode");
  //    Input(fineLevel, "UnAmalgamationInfo");
  Input(fineLevel, "CoarseMap");
  Input(fineLevel, "FC Splitting");

  const ParameterList& pL = GetParameterList();
  std::string drop_algo   = pL.get<std::string>("aggregation: drop scheme");
  if (drop_algo.find("block diagonal") != std::string::npos || drop_algo == "signed classical") {
    Input(fineLevel, "BlockNumber");
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);
  // using STS = Teuchos::ScalarTraits<SC>;

  // We start by assuming that someone did a reasonable strength of connection
  // algorithm before we start to get our Graph, DofsPerNode and UnAmalgamationInfo

  // We begin by getting a MIS (from a graph coloring) and then at that point we need
  // to start generating entries for the prolongator.
  RCP<const Matrix> A                              = Get<RCP<Matrix> >(fineLevel, "A");
  RCP<const Map> ownedCoarseMap                    = Get<RCP<const Map> >(fineLevel, "CoarseMap");
  RCP<const LocalOrdinalVector> owned_fc_splitting = Get<RCP<LocalOrdinalVector> >(fineLevel, "FC Splitting");
  RCP<const GraphBase> graph                       = Get<RCP<GraphBase> >(fineLevel, "Graph");
  //    LO nDofsPerNode                 = Get<LO>(fineLevel, "DofsPerNode");
  //    RCP<AmalgamationInfo> amalgInfo = Get< RCP<AmalgamationInfo> >     (fineLevel, "UnAmalgamationInfo");
  RCP<const Import> Importer = A->getCrsGraph()->getImporter();
  Xpetra::UnderlyingLib lib  = ownedCoarseMap->lib();

  //    RCP<MultiVector> fineNullspace = Get< RCP<MultiVector> > (fineLevel, "Nullspace");
  RCP<Matrix> P;
  //    SC SC_ZERO = STS::zero();
  LO LO_INVALID           = Teuchos::OrdinalTraits<LO>::invalid();
  const point_type C_PT   = ClassicalMapFactory::C_PT;
  const point_type F_PT   = ClassicalMapFactory::F_PT;
  const ParameterList& pL = GetParameterList();

  // FIXME: This guy doesn't work right now for NumPDEs != 1
  TEUCHOS_TEST_FOR_EXCEPTION(A->GetFixedBlockSize() != 1, Exceptions::RuntimeError, "ClassicalPFactory: Multiple PDEs per node not supported yet");

  // NOTE: Let's hope we never need to deal with this case
  TEUCHOS_TEST_FOR_EXCEPTION(!A->getRowMap()->isSameAs(*A->getDomainMap()), Exceptions::RuntimeError, "ClassicalPFactory: MPI Ranks > 1 not supported yet");

  // Do we need ghosts rows of A and myPointType?
  std::string scheme   = pL.get<std::string>("aggregation: classical scheme");
  bool need_ghost_rows = false;
  if (scheme == "ext+i")
    need_ghost_rows = true;
  else if (scheme == "direct")
    need_ghost_rows = false;
  else if (scheme == "classical modified")
    need_ghost_rows = true;
  // NOTE: ParameterList validator will check this guy so we don't really need an "else" here

  if (GetVerbLevel() & Statistics1) {
    GetOStream(Statistics1) << "ClassicalPFactory: scheme = " << scheme << std::endl;
  }

  // Ghost the FC splitting and grab the data (if needed)
  RCP<const LocalOrdinalVector> fc_splitting;
  ArrayRCP<const LO> myPointType;
  if (Importer.is_null()) {
    fc_splitting = owned_fc_splitting;
  } else {
    RCP<LocalOrdinalVector> fc_splitting_nonconst = LocalOrdinalVectorFactory::Build(A->getCrsGraph()->getColMap());
    fc_splitting_nonconst->doImport(*owned_fc_splitting, *Importer, Xpetra::INSERT);
    fc_splitting = fc_splitting_nonconst;
  }
  myPointType = fc_splitting->getData(0);

  /* Ghost A (if needed) */
  RCP<const Matrix> Aghost;
  RCP<const LocalOrdinalVector> fc_splitting_ghost;
  ArrayRCP<const LO> myPointType_ghost;
  RCP<const Import> remoteOnlyImporter;
  if (need_ghost_rows && !Importer.is_null()) {
    ArrayView<const LO> remoteLIDs = Importer->getRemoteLIDs();
    size_t numRemote               = Importer->getNumRemoteIDs();
    Array<GO> remoteRows(numRemote);
    for (size_t i = 0; i < numRemote; i++)
      remoteRows[i] = Importer->getTargetMap()->getGlobalElement(remoteLIDs[i]);

    RCP<const Map> remoteRowMap = MapFactory::Build(lib, Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(), remoteRows(),
                                                    A->getDomainMap()->getIndexBase(), A->getDomainMap()->getComm());

    remoteOnlyImporter        = Importer->createRemoteOnlyImport(remoteRowMap);
    RCP<const CrsMatrix> Acrs = rcp_dynamic_cast<const CrsMatrixWrap>(A)->getCrsMatrix();
    RCP<CrsMatrix> Aghost_crs = CrsMatrixFactory::Build(Acrs, *remoteOnlyImporter, A->getDomainMap(), remoteOnlyImporter->getTargetMap());
    Aghost                    = rcp(new CrsMatrixWrap(Aghost_crs));
    // CAG: It seems that this isn't actually needed?
    // We also may need need to ghost myPointType for Aghost
    // RCP<const Import> Importer2 = Aghost->getCrsGraph()->getImporter();
    // if(Importer2.is_null()) {
    //   RCP<LocalOrdinalVector> fc_splitting_ghost_nonconst = LocalOrdinalVectorFactory::Build(Aghost->getColMap());
    //   fc_splitting_ghost_nonconst->doImport(*owned_fc_splitting,*Importer,Xpetra::INSERT);
    //   fc_splitting_ghost = fc_splitting_ghost_nonconst;
    //   myPointType_ghost  = fc_splitting_ghost->getData(0);
    // }
  }

  /* Generate the ghosted Coarse map using the "Tuminaro maneuver" (if needed)*/
  RCP<const Map> coarseMap;
  if (Importer.is_null())
    coarseMap = ownedCoarseMap;
  else {
    // Generate a domain vector with the coarse ID's as entries for C points
    GhostCoarseMap(*A, *Importer, myPointType, ownedCoarseMap, coarseMap);
  }

  /* Get the block number, if we need it (and ghost it) */
  RCP<LocalOrdinalVector> BlockNumber;
  std::string drop_algo = pL.get<std::string>("aggregation: drop scheme");
  if (drop_algo.find("block diagonal") != std::string::npos || drop_algo == "signed classical") {
    RCP<LocalOrdinalVector> OwnedBlockNumber;
    OwnedBlockNumber = Get<RCP<LocalOrdinalVector> >(fineLevel, "BlockNumber");
    if (Importer.is_null())
      BlockNumber = OwnedBlockNumber;
    else {
      BlockNumber = LocalOrdinalVectorFactory::Build(A->getColMap());
      BlockNumber->doImport(*OwnedBlockNumber, *Importer, Xpetra::INSERT);
    }
  }

#if defined(CMS_DEBUG) || defined(CMS_DUMP)
  {
    RCP<const CrsMatrix> Acrs = rcp_dynamic_cast<const CrsMatrixWrap>(A)->getCrsMatrix();
    int rank                  = A->getRowMap()->getComm()->getRank();
    printf("[%d] A local size = %dx%d\n", rank, (int)Acrs->getRowMap()->getLocalNumElements(), (int)Acrs->getColMap()->getLocalNumElements());

    printf("[%d] graph local size = %dx%d\n", rank, (int)graph->GetDomainMap()->getLocalNumElements(), (int)graph->GetImportMap()->getLocalNumElements());

    std::ofstream ofs(std::string("dropped_graph_") + std::to_string(fineLevel.GetLevelID()) + std::string("_") + std::to_string(rank) + std::string(".dat"), std::ofstream::out);
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(ofs));
    graph->print(*fancy, Debug);
    std::string out_fc    = std::string("fc_splitting_") + std::to_string(fineLevel.GetLevelID()) + std::string("_") + std::to_string(rank) + std::string(".dat");
    std::string out_block = std::string("block_number_") + std::to_string(fineLevel.GetLevelID()) + std::string("_") + std::to_string(rank) + std::string(".dat");

    // We don't support writing LO vectors in Xpetra (boo!) so....
    using real_type             = typename Teuchos::ScalarTraits<SC>::magnitudeType;
    using RealValuedMultiVector = typename Xpetra::MultiVector<real_type, LO, GO, NO>;
    typedef Xpetra::MultiVectorFactory<real_type, LO, GO, NO> RealValuedMultiVectorFactory;

    RCP<RealValuedMultiVector> mv = RealValuedMultiVectorFactory::Build(fc_splitting->getMap(), 1);
    ArrayRCP<real_type> mv_data   = mv->getDataNonConst(0);

    // FC Splitting
    ArrayRCP<const LO> fc_data = fc_splitting->getData(0);
    for (LO i = 0; i < (LO)fc_data.size(); i++)
      mv_data[i] = Teuchos::as<real_type>(fc_data[i]);
    Xpetra::IO<real_type, LO, GO, NO>::Write(out_fc, *mv);

    // Block Number
    if (!BlockNumber.is_null()) {
      RCP<RealValuedMultiVector> mv2 = RealValuedMultiVectorFactory::Build(BlockNumber->getMap(), 1);
      ArrayRCP<real_type> mv_data2   = mv2->getDataNonConst(0);
      ArrayRCP<const LO> b_data      = BlockNumber->getData(0);
      for (LO i = 0; i < (LO)b_data.size(); i++) {
        mv_data2[i] = Teuchos::as<real_type>(b_data[i]);
      }
      Xpetra::IO<real_type, LO, GO, NO>::Write(out_block, *mv2);
    }
  }

#endif

  /* Generate reindexing arrays */
  // Note: cpoint2pcol is ghosted if myPointType is
  // NOTE: Since the ghosted coarse column map follows the ordering of
  // the fine column map, this *should* work, because it is in local indices.
  // pcol2cpoint - Takes a LCID for P and turns in into an LCID for A.
  // cpoint2pcol - Takes a LCID for A --- if it is a C Point --- and turns in into an LCID for P.
  Array<LO> cpoint2pcol(myPointType.size(), LO_INVALID);
  Array<LO> pcol2cpoint(coarseMap->getLocalNumElements(), LO_INVALID);
  LO num_c_points = 0;
  LO num_f_points = 0;
  for (LO i = 0; i < (LO)myPointType.size(); i++) {
    if (myPointType[i] == C_PT) {
      cpoint2pcol[i] = num_c_points;
      num_c_points++;
    } else if (myPointType[i] == F_PT)
      num_f_points++;
  }
  for (LO i = 0; i < (LO)cpoint2pcol.size(); i++) {
    if (cpoint2pcol[i] != LO_INVALID)
      pcol2cpoint[cpoint2pcol[i]] = i;
  }

  // Generate edge strength flags (this will make everything easier later)
  // These do *not* need to be ghosted (unlike A)
  Teuchos::Array<size_t> eis_rowptr;
  Teuchos::Array<bool> edgeIsStrong;
  {
    SubFactoryMonitor sfm(*this, "Strength Flags", coarseLevel);
    GenerateStrengthFlags(*A, *graph, eis_rowptr, edgeIsStrong);
  }

  // Phase 3: Generate the P matrix
  RCP<const Map> coarseColMap    = coarseMap;
  RCP<const Map> coarseDomainMap = ownedCoarseMap;
  if (scheme == "ext+i") {
    SubFactoryMonitor sfm(*this, "Ext+i Interpolation", coarseLevel);
    Coarsen_Ext_Plus_I(*A, Aghost, *graph, coarseColMap, coarseDomainMap, num_c_points, num_f_points, myPointType(), myPointType_ghost(), cpoint2pcol, pcol2cpoint, eis_rowptr, edgeIsStrong, BlockNumber, P);
  } else if (scheme == "direct") {
    SubFactoryMonitor sfm(*this, "Direct Interpolation", coarseLevel);
    Coarsen_Direct(*A, Aghost, *graph, coarseColMap, coarseDomainMap, num_c_points, num_f_points, myPointType(), myPointType_ghost(), cpoint2pcol, pcol2cpoint, eis_rowptr, edgeIsStrong, BlockNumber, P);
  } else if (scheme == "classical modified") {
    SubFactoryMonitor sfm(*this, "Classical Modified Interpolation", coarseLevel);
    Coarsen_ClassicalModified(*A, Aghost, *graph, coarseColMap, coarseDomainMap, num_c_points, num_f_points, myPointType(), myPointType_ghost(), cpoint2pcol, pcol2cpoint, eis_rowptr, edgeIsStrong, BlockNumber, remoteOnlyImporter, P);
  }
  // NOTE: ParameterList validator will check this guy so we don't really need an "else" here

#ifdef CMS_DEBUG
  Xpetra::IO<SC, LO, GO, NO>::Write("classical_p.mat." + std::to_string(coarseLevel.GetLevelID()), *P);
  // Xpetra::IO<SC,LO,GO,NO>::WriteLocal("classical_p.mat."+std::to_string(coarseLevel.GetLevelID()), *P);
#endif

  // Save output
  Set(coarseLevel, "P", P);
  RCP<const CrsGraph> pg = P->getCrsGraph();
  Set(coarseLevel, "P Graph", pg);

  // RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap, fineNullspace->getNumVectors());
  //     P->apply(*fineNullspace, *coarseNullspace, Teuchos::TRANS, Teuchos::ScalarTraits<SC>::one(), Teuchos::ScalarTraits<SC>::zero());
  //     Set(coarseLevel, "Nullspace", coarseNullspace);

  if (IsPrint(Statistics1)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics1) << PerfUtils::PrintMatrixInfo(*P, "P", params);
  }
}
/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Coarsen_ClassicalModified(const Matrix& A, const RCP<const Matrix>& Aghost, const GraphBase& graph, RCP<const Map>& coarseColMap, RCP<const Map>& coarseDomainMap, LO num_c_points, LO num_f_points, const Teuchos::ArrayView<const LO>& myPointType, const Teuchos::ArrayView<const LO>& myPointType_ghost, const Teuchos::Array<LO>& cpoint2pcol, const Teuchos::Array<LO>& pcol2cpoint, Teuchos::Array<size_t>& eis_rowptr, Teuchos::Array<bool>& edgeIsStrong, RCP<LocalOrdinalVector>& BlockNumber, RCP<const Import> remoteOnlyImporter, RCP<Matrix>& P) const {
  /* ============================================================= */
  /* Phase 3 : Classical Modified Interpolation                    */
  /* De Sterck, Falgout, Nolting and Yang. "Distance-two           */
  /* interpolation for parallel algebraic multigrid", NLAA 2008    */
  /* 15:115-139                                                    */
  /* ============================================================= */
  /* Definitions:                                                        */
  /* F = F-points                                                        */
  /* C = C-points                                                        */
  /* N_i = non-zero neighbors of node i                                  */
  /* S_i = {j\in N_i | j strongly influences i } [strong neighbors of i] */
  /* F_i^s = F \cap S_i [strong F-neighbors of i]                        */
  /* C_i^s = C \cap S_i [strong C-neighbors of i]                        */

  /* N_i^w = N_i\ (F_i^s \cup C_i^s) [weak neighbors of i]               */
  /*         This guy has a typo.  The paper had a \cap instead of \cup  */
  /*         I would note that this set can contain both F-points and    */
  /*         C-points.  They're just weak neighbors of this guy.         */
  /*         Note that N_i^w \cup F_i^s \cup C_i^s = N_i by construction */

  /* \bar{a}_ij = {    0, if sign(a_ij) == sign(a_ii)                    */
  /*              { a_ij, otherwise                                      */

  /* F_i^s\star = {k\in N_i | C_i^s \cap C_k^s = \emptyset}              */
  /*              [set of F-neighbors of i that do not share a strong    */
  /*               C-neighbor with i]                                    */

  /* Rewritten Equation (9) on p. 120                                    */
  /* \tilde{a}_ii =  (a_ij + \sum_{k\in{N_i^w \cup F_i^s\star}} a_ik     */
  /*                                                                     */
  /* f_ij = \sum_{k\in{F_i^s\setminusF_i^s*}} \frac{a_ik \bar{a}_kj}{\sum_{m\inC_i^s \bar{a}_km}}    */
  /*                                                                     */
  /* w_ij = \frac{1}{\tilde{a}_ii} ( a_ij + f_ij)  for all j in C_i^s    */

  //    const point_type F_PT = ClassicalMapFactory::F_PT;
  const point_type C_PT         = ClassicalMapFactory::C_PT;
  const point_type DIRICHLET_PT = ClassicalMapFactory::DIRICHLET_PT;
  using STS                     = typename Teuchos::ScalarTraits<SC>;
  LO LO_INVALID                 = Teuchos::OrdinalTraits<LO>::invalid();
  //    size_t ST_INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  SC SC_ZERO = STS::zero();
#ifdef CMS_DEBUG
  int rank = A.getRowMap()->getComm()->getRank();
#endif

  // Get the block number if we have it.
  ArrayRCP<const LO> block_id;
  if (!BlockNumber.is_null())
    block_id = BlockNumber->getData(0);

  // Initial (estimated) allocation
  // NOTE: If we only used Tpetra, then we could use these guys as is, but because Epetra, we can't, so there
  // needs to be a copy below.
  size_t Nrows                         = A.getLocalNumRows();
  double c_point_density               = (double)num_c_points / (num_c_points + num_f_points);
  double mean_strong_neighbors_per_row = (double)graph.GetNodeNumEdges() / graph.GetNodeNumVertices();
  //    double mean_neighbors_per_row = (double)A.getLocalNumEntries() / Nrows;
  double nnz_per_row_est = c_point_density * mean_strong_neighbors_per_row;

  size_t nnz_est = std::max(Nrows, std::min((size_t)A.getLocalNumEntries(), (size_t)(nnz_per_row_est * Nrows)));
  Array<size_t> tmp_rowptr(Nrows + 1);
  Array<LO> tmp_colind(nnz_est);

  // Algorithm (count+realloc)
  // For each row, i,
  // - Count the number of elements in \hat{C}_j, aka [C-neighbors and C-neighbors of strong F-neighbors of i]
  size_t ct = 0;
  for (LO row = 0; row < (LO)Nrows; row++) {
    size_t row_start = eis_rowptr[row];
    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    std::set<LO> C_hat;
    if (myPointType[row] == DIRICHLET_PT) {
      // Dirichlet points get ignored completely
    } else if (myPointType[row] == C_PT) {
      // C-Points get a single 1 in their row
      C_hat.insert(cpoint2pcol[row]);
    } else {
      // F-Points have a more complicated interpolation

      // C-neighbors of row
      A.getLocalRowView(row, indices, vals);
      for (LO j = 0; j < indices.size(); j++)
        if (myPointType[indices[j]] == C_PT && edgeIsStrong[row_start + j])
          C_hat.insert(cpoint2pcol[indices[j]]);
    }  // end else

    // Realloc if needed
    if (ct + (size_t)C_hat.size() > (size_t)tmp_colind.size()) {
      tmp_colind.resize(std::max(ct + (size_t)C_hat.size(), (size_t)2 * tmp_colind.size()));
    }

    // Copy
    std::copy(C_hat.begin(), C_hat.end(), tmp_colind.begin() + ct);
    ct += C_hat.size();
    tmp_rowptr[row + 1] = tmp_rowptr[row] + C_hat.size();
  }
  // Resize down
  tmp_colind.resize(tmp_rowptr[Nrows]);

  RCP<CrsMatrix> Pcrs = CrsMatrixFactory::Build(A.getRowMap(), coarseColMap, 0);
  ArrayRCP<size_t> P_rowptr;
  ArrayRCP<LO> P_colind;
  ArrayRCP<SC> P_values;

  if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
    P_rowptr = Teuchos::arcpFromArray(tmp_rowptr);
    P_colind = Teuchos::arcpFromArray(tmp_colind);
    P_values.resize(P_rowptr[Nrows]);
  } else {
    // Make the matrix and then get the graph out of it (necessary for Epetra)
    // NOTE: The lack of an Epetra_CrsGraph::ExpertStaticFillComplete makes this
    // impossible to do the obvious way
    Pcrs->allocateAllValues(tmp_rowptr[Nrows], P_rowptr, P_colind, P_values);
    std::copy(tmp_rowptr.begin(), tmp_rowptr.end(), P_rowptr.begin());
    std::copy(tmp_colind.begin(), tmp_colind.end(), P_colind.begin());
    Pcrs->setAllValues(P_rowptr, P_colind, P_values);
    Pcrs->expertStaticFillComplete(/*domain*/ coarseDomainMap, /*range*/ A.getDomainMap());
  }

  // Generate a remote-ghosted version of the graph (if we're in parallel)
  RCP<const CrsGraph> Pgraph;
  RCP<const CrsGraph> Pghost;
  // TODO: We might want to be more efficient here and actually use
  // Pgraph in the matrix constructor.
  ArrayRCP<const size_t> Pghost_rowptr;
  ArrayRCP<const LO> Pghost_colind;
  if (!remoteOnlyImporter.is_null()) {
    if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
      RCP<CrsGraph> tempPgraph = CrsGraphFactory::Build(A.getRowMap(), coarseColMap, P_rowptr, P_colind);
      tempPgraph->fillComplete(coarseDomainMap, A.getDomainMap());
      Pgraph = tempPgraph;
    } else {
      // Epetra does not have a graph constructor that uses rowptr and colind.
      Pgraph = Pcrs->getCrsGraph();
    }
    TEUCHOS_ASSERT(!Pgraph.is_null());
    Pghost = CrsGraphFactory::Build(Pgraph, *remoteOnlyImporter, Pgraph->getDomainMap(), remoteOnlyImporter->getTargetMap());
    Pghost->getAllIndices(Pghost_rowptr, Pghost_colind);
  }

  // Gustavson-style perfect hashing
  ArrayRCP<LO> Acol_to_Pcol(A.getColMap()->getLocalNumElements(), LO_INVALID);

  // Get a quick reindexing array from Pghost LCIDs to P LCIDs
  ArrayRCP<LO> Pghostcol_to_Pcol;
  if (!Pghost.is_null()) {
    Pghostcol_to_Pcol.resize(Pghost->getColMap()->getLocalNumElements(), LO_INVALID);
    for (LO i = 0; i < (LO)Pghost->getColMap()->getLocalNumElements(); i++)
      Pghostcol_to_Pcol[i] = Pgraph->getColMap()->getLocalElement(Pghost->getColMap()->getGlobalElement(i));
  }  // end Pghost

  // Get a quick reindexing array from Aghost LCIDs to A LCIDs
  ArrayRCP<LO> Aghostcol_to_Acol;
  if (!Aghost.is_null()) {
    Aghostcol_to_Acol.resize(Aghost->getColMap()->getLocalNumElements(), LO_INVALID);
    for (LO i = 0; i < (LO)Aghost->getColMap()->getLocalNumElements(); i++)
      Aghostcol_to_Acol[i] = A.getColMap()->getLocalElement(Aghost->getColMap()->getGlobalElement(i));
  }  // end Aghost

  // Algorithm (numeric)
  for (LO i = 0; i < (LO)Nrows; i++) {
    if (myPointType[i] == DIRICHLET_PT) {
      // Dirichlet points get ignored completely
#ifdef CMS_DEBUG
      // DEBUG
      printf("[%d] ** A(%d,:) is a Dirichlet-Point.\n", rank, i);
#endif
    } else if (myPointType[i] == C_PT) {
      // C Points get a single 1 in their row
      P_values[P_rowptr[i]] = Teuchos::ScalarTraits<SC>::one();
#ifdef CMS_DEBUG
      // DEBUG
      printf("[%d] ** A(%d,:) is a C-Point.\n", rank, i);
#endif
    } else {
      // F Points get all of the fancy stuff
#ifdef CMS_DEBUG
      // DEBUG
      printf("[%d] ** A(%d,:) is a F-Point.\n", rank, i);
#endif

      // Get all of the relevant information about this row
      ArrayView<const LO> A_indices_i, A_indices_k;
      ArrayView<const SC> A_vals_i, A_vals_k;
      A.getLocalRowView(i, A_indices_i, A_vals_i);
      size_t row_start = eis_rowptr[i];

      ArrayView<const LO> P_indices_i = P_colind.view(P_rowptr[i], P_rowptr[i + 1] - P_rowptr[i]);
      ArrayView<SC> P_vals_i          = P_values.view(P_rowptr[i], P_rowptr[i + 1] - P_rowptr[i]);

      // FIXME: Do we need this?
      for (LO j = 0; j < (LO)P_vals_i.size(); j++)
        P_vals_i[j] = SC_ZERO;

      // Stash the hash:  Flag any strong C-points with their index into P_colind
      // NOTE:  We'll consider any points that are LO_INVALID or less than P_rowptr[i] as not strong C-points
      for (LO j = 0; j < (LO)P_indices_i.size(); j++) {
        Acol_to_Pcol[pcol2cpoint[P_indices_i[j]]] = P_rowptr[i] + j;
      }

      // Loop over the entries in the row
      SC first_denominator = SC_ZERO;
#ifdef CMS_DEBUG
      SC a_ii = SC_ZERO;
#endif
      for (LO k0 = 0; k0 < (LO)A_indices_i.size(); k0++) {
        LO k      = A_indices_i[k0];
        SC a_ik   = A_vals_i[k0];
        LO pcol_k = Acol_to_Pcol[k];

        if (k == i) {
          // Case A: Diagonal value (add to first denominator)
          // FIXME:  Add BlockNumber matching here
          first_denominator += a_ik;
#ifdef CMS_DEBUG
          a_ii = a_ik;
          printf("- A(%d,%d) is the diagonal\n", i, k);
#endif

        } else if (myPointType[k] == DIRICHLET_PT) {
          // Case B: Ignore dirichlet connections completely
#ifdef CMS_DEBUG
          printf("- A(%d,%d) is a Dirichlet point\n", i, k);
#endif

        } else if (pcol_k != LO_INVALID && pcol_k >= (LO)P_rowptr[i]) {
          // Case C: a_ik is strong C-Point (goes directly into the weight)
          P_values[pcol_k] += a_ik;
#ifdef CMS_DEBUG
          printf("- A(%d,%d) is a strong C-Point\n", i, k);
#endif
        } else if (!edgeIsStrong[row_start + k0]) {
          // Case D: Weak non-Dirichlet neighbor (add to first denominator)
          if (block_id.size() == 0 || block_id[i] == block_id[k]) {
            first_denominator += a_ik;
#ifdef CMS_DEBUG
            printf("- A(%d,%d) is weak adding to diagonal(%d,%d) (%6.4e)\n", i, k, block_id[i], block_id[k], a_ik);
          } else {
            printf("- A(%d,%d) is weak but does not match blocks (%d,%d), discarding\n", i, k, block_id[i], block_id[k]);
#endif
          }

        } else {  // Case E
          // Case E: Strong F-Point (adds to the first denominator if we don't share a
          // a strong C-Point with i; adds to the second denominator otherwise)
#ifdef CMS_DEBUG
          printf("- A(%d,%d) is a strong F-Point\n", i, k);
#endif

          SC a_kk               = SC_ZERO;
          SC second_denominator = SC_ZERO;
          int sign_akk          = 0;

          if (k < (LO)Nrows) {
            // Grab the diagonal a_kk
            A.getLocalRowView(k, A_indices_k, A_vals_k);
            for (LO m0 = 0; m0 < (LO)A_indices_k.size(); m0++) {
              LO m = A_indices_k[m0];
              if (k == m) {
                a_kk = A_vals_k[m0];
                break;
              }
            }  // end for A_indices_k

            // Compute the second denominator term
            sign_akk = Sign(a_kk);
            for (LO m0 = 0; m0 < (LO)A_indices_k.size(); m0++) {
              LO m = A_indices_k[m0];
              if (m != k && Acol_to_Pcol[A_indices_k[m0]] >= (LO)P_rowptr[i]) {
                SC a_km = A_vals_k[m0];
                second_denominator += (Sign(a_km) == sign_akk ? SC_ZERO : a_km);
              }
            }  // end for A_indices_k

            // Now we have the second denominator, for this particular strong F point.
            // So we can now add the sum to the w_ij components for the P values
            if (second_denominator != SC_ZERO) {
              for (LO j0 = 0; j0 < (LO)A_indices_k.size(); j0++) {
                LO j = A_indices_k[j0];
                // NOTE: Row k should be in fis_star, so I should have to check for diagonals here
                //                  printf("Acol_to_Pcol[%d] = %d P_values.size() = %d\n",j,Acol_to_Pcol[j],(int)P_values.size());
                if (Acol_to_Pcol[j] >= (LO)P_rowptr[i]) {
                  SC a_kj         = A_vals_k[j0];
                  SC sign_akj_val = sign_akk == Sign(a_kj) ? SC_ZERO : a_kj;
                  P_values[Acol_to_Pcol[j]] += a_ik * sign_akj_val / second_denominator;
#ifdef CMS_DEBUG
                  printf("- - Unscaled P(%d,A-%d) += %6.4e = %5.4e\n", i, j, a_ik * sign_akj_val / second_denominator, P_values[Acol_to_Pcol[j]]);
#endif
                }
              }  // end for A_indices_k
            }    // end if second_denominator != 0
            else {
#ifdef CMS_DEBUG
              printf("- - A(%d,%d) second denominator is zero\n", i, k);
#endif
              if (block_id.size() == 0 || block_id[i] == block_id[k])
                first_denominator += a_ik;
            }  // end else second_denominator != 0
          }    // end if k < Nrows
          else {
            // Ghost row
            LO kless = k - Nrows;
            // Grab the diagonal a_kk
            // NOTE: ColMap is not locally fitted to the RowMap
            // so we need to check GIDs here
            Aghost->getLocalRowView(kless, A_indices_k, A_vals_k);
            GO k_g = Aghost->getRowMap()->getGlobalElement(kless);
            for (LO m0 = 0; m0 < (LO)A_indices_k.size(); m0++) {
              GO m_g = Aghost->getColMap()->getGlobalElement(A_indices_k[m0]);
              if (k_g == m_g) {
                a_kk = A_vals_k[m0];
                break;
              }
            }  // end for A_indices_k

            // Compute the second denominator term
            sign_akk = Sign(a_kk);
            for (LO m0 = 0; m0 < (LO)A_indices_k.size(); m0++) {
              GO m_g    = Aghost->getColMap()->getGlobalElement(A_indices_k[m0]);
              LO mghost = A_indices_k[m0];            // Aghost LCID
              LO m      = Aghostcol_to_Acol[mghost];  // A's LID (could be LO_INVALID)
              if (m_g != k_g && m != LO_INVALID && Acol_to_Pcol[m] >= (LO)P_rowptr[i]) {
                SC a_km = A_vals_k[m0];
                second_denominator += (Sign(a_km) == sign_akk ? SC_ZERO : a_km);
              }
            }  // end for A_indices_k

            // Now we have the second denominator, for this particular strong F point.
            // So we can now add the sum to the w_ij components for the P values
            if (second_denominator != SC_ZERO) {
              for (LO j0 = 0; j0 < (LO)A_indices_k.size(); j0++) {
                LO jghost = A_indices_k[j0];            // Aghost LCID
                LO j      = Aghostcol_to_Acol[jghost];  // A's LID (could be LO_INVALID)
                // NOTE: Row k should be in fis_star, so I should have to check for diagonals here
                if ((j != LO_INVALID) && (Acol_to_Pcol[j] >= (LO)P_rowptr[i])) {
                  SC a_kj         = A_vals_k[j0];
                  SC sign_akj_val = sign_akk == Sign(a_kj) ? SC_ZERO : a_kj;
                  P_values[Acol_to_Pcol[j]] += a_ik * sign_akj_val / second_denominator;
#ifdef CMS_DEBUG
                  printf("- - Unscaled P(%d,A-%d) += %6.4e\n", i, j, a_ik * sign_akj_val / second_denominator);
#endif
                }

              }  // end for A_indices_k
            }    // end if second_denominator != 0
            else {
#ifdef CMS_DEBUG
              printf("- - A(%d,%d) second denominator is zero\n", i, k);
#endif
              if (block_id.size() == 0 || block_id[i] == block_id[k])
                first_denominator += a_ik;
            }  // end else second_denominator != 0
          }    // end else k < Nrows
        }      // end else Case A,...,E

      }  // end for A_indices_i

      // Now, downscale by the first_denominator
      if (first_denominator != SC_ZERO) {
        for (LO j0 = 0; j0 < (LO)P_indices_i.size(); j0++) {
#ifdef CMS_DEBUG
          SC old_pij = P_vals_i[j0];
          P_vals_i[j0] /= -first_denominator;
          printf("P(%d,%d) = %6.4e = %6.4e / (%6.4e + %6.4e)\n", i, P_indices_i[j0], P_vals_i[j0], old_pij, a_ii, first_denominator - a_ii);
#else
          P_vals_i[j0] /= -first_denominator;
#endif
        }  // end for P_indices_i
      }    // end if first_denominator != 0

    }  // end else C-Point

  }  // end if i < Nrows

  // Finish up
  Pcrs->setAllValues(P_rowptr, P_colind, P_values);
  Pcrs->expertStaticFillComplete(/*domain*/ coarseDomainMap, /*range*/ A.getDomainMap());
  // Wrap from CrsMatrix to Matrix
  P = rcp(new CrsMatrixWrap(Pcrs));

}  // end Coarsen_ClassicalModified

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Coarsen_Direct(const Matrix& A, const RCP<const Matrix>& Aghost, const GraphBase& graph, RCP<const Map>& coarseColMap, RCP<const Map>& coarseDomainMap, LO num_c_points, LO num_f_points, const Teuchos::ArrayView<const LO>& myPointType, const Teuchos::ArrayView<const LO>& myPointType_ghost, const Teuchos::Array<LO>& cpoint2pcol, const Teuchos::Array<LO>& pcol2cpoint, Teuchos::Array<size_t>& eis_rowptr, Teuchos::Array<bool>& edgeIsStrong, RCP<LocalOrdinalVector>& BlockNumber, RCP<Matrix>& P) const {
  /* ============================================================= */
  /* Phase 3 : Direct Interpolation                                */
  /* We do not use De Sterck, Falgout, Nolting and Yang (2008)     */
  /* here.  Instead we follow:                                     */
  /* Trottenberg, Oosterlee and Schueller, Multigrid, 2001.        */
  /* with some modifications inspirted by PyAMG                    */
  /* ============================================================= */
  /* Definitions:                                                        */
  /* F = F-points                                                        */
  /* C = C-points                                                        */
  /* N_i = non-zero neighbors of node i                                  */
  /* S_i = {j\in N_i | j strongly influences i } [strong neighbors of i] */
  /* F_i^s = F \cap S_i [strong F-neighbors of i]                        */
  /* C_i^s = C \cap S_i [strong C-neighbors of i]                        */
  /* P_i = Set of interpolatory variables for row i [here = C_i^s]       */

  /* (A.2.17) from p. 426                                                */
  /* a_ij^- = {  a_ij,  if a_ij < 0                                      */
  /*          {     0,  otherwise                                        */
  /* a_ij^+ = {  a_ij,  if a_ij > 0                                      */
  /*          {     0,  otherwise                                        */
  /* P_i^- =  P_i \cap {k | a_ij^- != 0 and a_ij^- = a_ij}               */
  /*          [strong C-neighbors with negative edges]                   */
  /* P_i^+ =  P_i \cap {k | a_ij^+ != 0 and a_ij^+ = a_ij}               */
  /*          [strong C-neighbors with positive edges]                   */

  /* de Sterck et al., gives us this:                                                      */
  /* Rewritten Equation (6) on p. 119                                                      */
  /* w_ij = - a_ji / a_ii \frac{\sum_{k\in N_i} a_ik} {\sum k\inC_i^s} a_ik},   j\in C_i^s */

  /* Trottenberg et al. (A.7.6) and (A.7.7) on p. 479 gives this:                          */
  /* alpha_i = \frac{ \sum_{j\in N_i} a_ij^- }{ \sum_{k\in P_i} a_ik^- }                   */
  /* beta_i  = \frac{ \sum_{j\in N_i} a_ij^+ }{ \sum_{k\in P_i} a_ik^+ }                   */
  /* w_ik    = { - alpha_i (a_ik / a_ii),   if k\in P_i^-                                  */
  /*           { -  beta_i (a_ik / a_ii),   if k\in P_i^+                                  */
  /* NOTE: The text says to modify, if  P_i^+ is zero but it isn't entirely clear how that */
  /* works.  We'll follow the PyAMG implementation in a few important ways.                */

  const point_type C_PT         = ClassicalMapFactory::C_PT;
  const point_type DIRICHLET_PT = ClassicalMapFactory::DIRICHLET_PT;

  // Initial (estimated) allocation
  // NOTE: If we only used Tpetra, then we could use these guys as is, but because Epetra, we can't, so there
  // needs to be a copy below.
  using STS                            = typename Teuchos::ScalarTraits<SC>;
  using MT                             = typename STS::magnitudeType;
  using MTS                            = typename Teuchos::ScalarTraits<MT>;
  size_t Nrows                         = A.getLocalNumRows();
  double c_point_density               = (double)num_c_points / (num_c_points + num_f_points);
  double mean_strong_neighbors_per_row = (double)graph.GetNodeNumEdges() / graph.GetNodeNumVertices();
  //    double mean_neighbors_per_row = (double)A.getLocalNumEntries() / Nrows;
  double nnz_per_row_est = c_point_density * mean_strong_neighbors_per_row;

  size_t nnz_est = std::max(Nrows, std::min((size_t)A.getLocalNumEntries(), (size_t)(nnz_per_row_est * Nrows)));
  SC SC_ZERO     = STS::zero();
  MT MT_ZERO     = MTS::zero();
  Array<size_t> tmp_rowptr(Nrows + 1);
  Array<LO> tmp_colind(nnz_est);

  // Algorithm (count+realloc)
  // For each row, i,
  // - Count the number of elements in \hat{C}_j, aka [C-neighbors and C-neighbors of strong F-neighbors of i]
  size_t ct = 0;
  for (LO row = 0; row < (LO)Nrows; row++) {
    size_t row_start = eis_rowptr[row];
    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    std::set<LO> C_hat;
    if (myPointType[row] == DIRICHLET_PT) {
      // Dirichlet points get ignored completely
    } else if (myPointType[row] == C_PT) {
      // C-Points get a single 1 in their row
      C_hat.insert(cpoint2pcol[row]);
    } else {
      // F-Points have a more complicated interpolation

      // C-neighbors of row
      A.getLocalRowView(row, indices, vals);
      for (LO j = 0; j < indices.size(); j++)
        if (myPointType[indices[j]] == C_PT && edgeIsStrong[row_start + j])
          C_hat.insert(cpoint2pcol[indices[j]]);
    }  // end else

    // Realloc if needed
    if (ct + (size_t)C_hat.size() > (size_t)tmp_colind.size()) {
      tmp_colind.resize(std::max(ct + (size_t)C_hat.size(), (size_t)2 * tmp_colind.size()));
    }

    // Copy
    std::copy(C_hat.begin(), C_hat.end(), tmp_colind.begin() + ct);
    ct += C_hat.size();
    tmp_rowptr[row + 1] = tmp_rowptr[row] + C_hat.size();
  }
  // Resize down
  tmp_colind.resize(tmp_rowptr[Nrows]);

  // Allocate memory & copy
  P                   = rcp(new CrsMatrixWrap(A.getRowMap(), coarseColMap, 0));
  RCP<CrsMatrix> PCrs = rcp_dynamic_cast<CrsMatrixWrap>(P)->getCrsMatrix();
  ArrayRCP<size_t> P_rowptr;
  ArrayRCP<LO> P_colind;
  ArrayRCP<SC> P_values;

#ifdef CMS_DEBUG
  printf("CMS: Allocating P w/ %d nonzeros\n", (int)tmp_rowptr[Nrows]);
#endif
  PCrs->allocateAllValues(tmp_rowptr[Nrows], P_rowptr, P_colind, P_values);
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_rowptr.size() != P_rowptr.size(), Exceptions::RuntimeError, "ClassicalPFactory: Allocation size error (rowptr)");
  TEUCHOS_TEST_FOR_EXCEPTION(tmp_colind.size() != P_colind.size(), Exceptions::RuntimeError, "ClassicalPFactory: Allocation size error (colind)");
  // FIXME:  This can be short-circuited for Tpetra, if we decide we care
  for (LO i = 0; i < (LO)Nrows + 1; i++)
    P_rowptr[i] = tmp_rowptr[i];
  for (LO i = 0; i < (LO)tmp_rowptr[Nrows]; i++)
    P_colind[i] = tmp_colind[i];

  // Algorithm (numeric)
  for (LO i = 0; i < (LO)Nrows; i++) {
    if (myPointType[i] == DIRICHLET_PT) {
      // Dirichlet points get ignored completely
#ifdef CMS_DEBUG
      // DEBUG
      printf("** A(%d,:) is a Dirichlet-Point.\n", i);
#endif
    } else if (myPointType[i] == C_PT) {
      // C Points get a single 1 in their row
      P_values[P_rowptr[i]] = Teuchos::ScalarTraits<SC>::one();
#ifdef CMS_DEBUG
      // DEBUG
      printf("** A(%d,:) is a C-Point.\n", i);
#endif
    } else {
      /* Trottenberg et al. (A.7.6) and (A.7.7) on p. 479 gives this:                          */
      /* alpha_i = \frac{ \sum_{j\in N_i} a_ij^- }{ \sum_{k\in P_i} a_ik^- }                   */
      /* beta_i  = \frac{ \sum_{j\in N_i} a_ij^+ }{ \sum_{k\in P_i} a_ik^+ }                   */
      /* w_ik    = { - alpha_i (a_ik / a_ii),   if k\in P_i^-                                  */
      /*           { -  beta_i (a_ik / a_ii),   if k\in P_i^+                                  */
      ArrayView<const LO> A_indices_i, A_incides_k;
      ArrayView<const SC> A_vals_i, A_indices_k;
      A.getLocalRowView(i, A_indices_i, A_vals_i);
      size_t row_start = eis_rowptr[i];

      ArrayView<LO> P_indices_i = P_colind.view(P_rowptr[i], P_rowptr[i + 1] - P_rowptr[i]);
      ArrayView<SC> P_vals_i    = P_values.view(P_rowptr[i], P_rowptr[i + 1] - P_rowptr[i]);

#ifdef CMS_DEBUG
      // DEBUG
      {
        char mylabel[5] = "FUCD";
        char sw[3]      = "ws";
        printf("** A(%d,:) = ", i);
        for (LO j = 0; j < (LO)A_indices_i.size(); j++) {
          printf("%6.4e(%d-%c%c) ", A_vals_i[j], A_indices_i[j], mylabel[1 + myPointType[A_indices_i[j]]], sw[(int)edgeIsStrong[row_start + j]]);
        }
        printf("\n");
      }
#endif

      SC a_ii          = SC_ZERO;
      SC pos_numerator = SC_ZERO, neg_numerator = SC_ZERO;
      SC pos_denominator = SC_ZERO, neg_denominator = SC_ZERO;
      // Find the diagonal and compute the sum ratio
      for (LO j = 0; j < (LO)A_indices_i.size(); j++) {
        SC a_ik = A_vals_i[j];
        LO k    = A_indices_i[j];

        // Diagonal
        if (i == k) {
          a_ii = a_ik;
        }
        // Only strong C-neighbors are in the denomintor
        if (myPointType[k] == C_PT && edgeIsStrong[row_start + j]) {
          if (STS::real(a_ik) > MT_ZERO)
            pos_denominator += a_ik;
          else
            neg_denominator += a_ik;
        }

        // All neighbors are in the numerator
        // NOTE: As per PyAMG, this does not include the diagonal
        if (i != k) {
          if (STS::real(a_ik) > MT_ZERO)
            pos_numerator += a_ik;
          else
            neg_numerator += a_ik;
        }
      }
      SC alpha = (neg_denominator == MT_ZERO) ? SC_ZERO : (neg_numerator / neg_denominator);
      SC beta  = (pos_denominator == MT_ZERO) ? SC_ZERO : (pos_numerator / pos_denominator);
      alpha /= -a_ii;
      beta /= -a_ii;

      // Loop over the entries
      for (LO p_j = 0; p_j < (LO)P_indices_i.size(); p_j++) {
        LO P_col = pcol2cpoint[P_indices_i[p_j]];
        SC a_ij  = SC_ZERO;

        // Find A_ij (if it is there)
        // FIXME: We can optimize this if we assume sorting
        for (LO a_j = 0; a_j < (LO)A_indices_i.size(); a_j++) {
          if (A_indices_i[a_j] == P_col) {
            a_ij = A_vals_i[a_j];
            break;
          }
        }
        SC w_ij = (STS::real(a_ij) < 0) ? (alpha * a_ij) : (beta * a_ij);
#ifdef CMS_DEBUG
        SC alpha_or_beta = (STS::real(a_ij) < 0) ? alpha : beta;
        printf("P(%d,%d/%d) =  - %6.4e  * %6.4e  = %6.4e\n", i, P_indices_i[p_j], pcol2cpoint[P_indices_i[p_j]], alpha_or_beta, a_ij, w_ij);
#endif
        P_vals_i[p_j] = w_ij;
      }  // end for A_indices_i
    }    // end else C_PT
  }      // end for Numrows

  // Finish up
  PCrs->setAllValues(P_rowptr, P_colind, P_values);
  PCrs->expertStaticFillComplete(/*domain*/ coarseDomainMap, /*range*/ A.getDomainMap());
}

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    Coarsen_Ext_Plus_I(const Matrix& A, const RCP<const Matrix>& Aghost, const GraphBase& graph, RCP<const Map>& coarseColMap, RCP<const Map>& coarseDomainMap, LO num_c_points, LO num_f_points, const Teuchos::ArrayView<const LO>& myPointType, const Teuchos::ArrayView<const LO>& myPointType_ghost, const Teuchos::Array<LO>& cpoint2pcol, const Teuchos::Array<LO>& pcol2cpoint, Teuchos::Array<size_t>& eis_rowptr, Teuchos::Array<bool>& edgeIsStrong, RCP<LocalOrdinalVector>& BlockNumber, RCP<Matrix>& P) const {
  /* ============================================================= */
  /* Phase 3 : Extended+i Interpolation                            */
  /* De Sterck, Falgout, Nolting and Yang. "Distance-two           */
  /* interpolation for parallel algebraic multigrid", NLAA 2008    */
  /* 15:115-139                                                    */
  /* ============================================================= */
  /* Definitions:                                                        */
  /* F = F-points                                                        */
  /* C = C-points                                                        */
  /* N_i = non-zero neighbors of node i                                  */
  /* S_i = {j\in N_i | j strongly influences i } [strong neighbors of i] */
  /* F_i^s = F \cap S_i [strong F-neighbors of i]                        */
  /* C_i^s = C \cap S_i [strong C-neighbors of i]                        */
  /* N_i^w = N_i\ (F_i^s \cup C_i^s) [weak neighbors of i]               */
  /*         This guy has a typo.  The paper had a \cap instead of \cup  */
  /*         I would note that this set can contain both F-points and    */
  /*         C-points.  They're just weak neighbors of this guy.         */
  /*         Note that N_i^w \cup F_i^s \cup C_i^s = N_i by construction */

  /* \hat{C}_i = C_i \cup (\bigcup_{j\inF_i^s} C_j)                      */
  /*         [C-neighbors and C-neighbors of strong F-neighbors of i]    */
  /*                                                                     */

  /* \bar{a}_ij = {    0, if sign(a_ij) == sign(a_ii)                    */
  /*              { a_ij, otherwise                                      */

  /* Rewritten Equation (19) on p. 123                                   */
  /* f_ik = \frac{\bar{a}_kj}{\sum{l\in \hat{C}_i\cup {i}} \bar{a}_kl    */
  /* w_ij = -\tilde{a}_ii^{-1} (a_ij + \sum_{k\inF_i^s} a_ik f_ik        */
  /*         for j in \hat{C}_i                                          */

  /* Rewritten Equation (20) on p. 124 [for the lumped diagonal]                                  */
  /* g_ik = \frac{\bar{a}_ki}{\sum{l\in \hat{C}_i\cup {i}} \bar{a}_kl                             */
  /* \tilde{a}_ii = a_ii + \sum_{n\inN_i^w\setminus \hat{C}_i} a_in + \sum_{k\inF_i^s} a_ik g_ik  */
  TEUCHOS_TEST_FOR_EXCEPTION(1, std::runtime_error, "ClassicalPFactory: Ext+i not implemented");
}

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GenerateStrengthFlags(const Matrix& A, const GraphBase& graph, Teuchos::Array<size_t>& eis_rowptr, Teuchos::Array<bool>& edgeIsStrong) const {
  // To make this easier, we'll create a bool array equal to the nnz in the matrix
  // so we know whether each edge is strong or not.  This will save us a bunch of
  // trying to match the graph and matrix later
  size_t Nrows = A.getLocalNumRows();
  eis_rowptr.resize(Nrows + 1);

  if (edgeIsStrong.size() == 0) {
    // Preferred
    edgeIsStrong.resize(A.getLocalNumEntries(), false);
  } else {
    edgeIsStrong.resize(A.getLocalNumEntries(), false);
    edgeIsStrong.assign(A.getLocalNumEntries(), false);
  }

  eis_rowptr[0] = 0;
  for (LO i = 0; i < (LO)Nrows; i++) {
    LO rowstart = eis_rowptr[i];
    ArrayView<const LO> A_indices;
    ArrayView<const SC> A_values;
    A.getLocalRowView(i, A_indices, A_values);
    LO A_size = (LO)A_indices.size();

    ArrayView<const LO> G_indices = graph.getNeighborVertices(i);
    LO G_size                     = (LO)G_indices.size();

    // Both of these guys should be in the same (sorted) order, but let's check
    bool is_ok = true;
    for (LO j = 0; j < A_size - 1; j++)
      if (A_indices[j] >= A_indices[j + 1]) {
        is_ok = false;
        break;
      }
    for (LO j = 0; j < G_size - 1; j++)
      if (G_indices[j] >= G_indices[j + 1]) {
        is_ok = false;
        break;
      }
    TEUCHOS_TEST_FOR_EXCEPTION(!is_ok, Exceptions::RuntimeError, "ClassicalPFactory: Exected A and Graph to be sorted");

    // Now cycle through and set the flags - if the edge is in G it is strong
    for (LO g_idx = 0, a_idx = 0; g_idx < G_size; g_idx++) {
      LO col = G_indices[g_idx];
      while (A_indices[a_idx] != col && a_idx < A_size) a_idx++;
      if (a_idx == A_size) {
        is_ok = false;
        break;
      }
      edgeIsStrong[rowstart + a_idx] = true;
    }

    eis_rowptr[i + 1] = eis_rowptr[i] + A_size;
  }
}

/* ************************************************************************* */
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ClassicalPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    GhostCoarseMap(const Matrix& A, const Import& Importer, const ArrayRCP<const LO> myPointType, const RCP<const Map>& coarseMap, RCP<const Map>& coarseColMap) const {
  const point_type C_PT                = ClassicalMapFactory::C_PT;
  const GO GO_INVALID                  = Teuchos::OrdinalTraits<GO>::invalid();
  RCP<GlobalOrdinalVector> d_coarseIds = GlobalOrdinalVectorFactory::Build(A.getRowMap());
  ArrayRCP<GO> d_data                  = d_coarseIds->getDataNonConst(0);
  LO ct                                = 0;

  for (LO i = 0; i < (LO)d_data.size(); i++) {
    if (myPointType[i] == C_PT) {
      d_data[i] = coarseMap->getGlobalElement(ct);
      ct++;
    } else
      d_data[i] = GO_INVALID;
  }

  // Ghost this guy
  RCP<GlobalOrdinalVector> c_coarseIds = GlobalOrdinalVectorFactory::Build(A.getColMap());
  c_coarseIds->doImport(*d_coarseIds, Importer, Xpetra::INSERT);

  // If we assume that A is in Aztec ordering, then any subset of A's unknowns will
  // be in Aztec ordering as well, which means we can just condense these guys down
  // Overallocate, count and view
  ArrayRCP<GO> c_data = c_coarseIds->getDataNonConst(0);

  Array<GO> c_gids(c_data.size());
  LO count = 0;

  for (LO i = 0; i < (LO)c_data.size(); i++) {
    if (c_data[i] != GO_INVALID) {
      c_gids[count] = c_data[i];
      count++;
    }
  }
  // FIXME: Assumes scalar PDE
  std::vector<size_t> stridingInfo_(1);
  stridingInfo_[0]   = 1;
  GO domainGIDOffset = 0;

  coarseColMap = StridedMapFactory::Build(coarseMap->lib(),
                                          Teuchos::OrdinalTraits<Xpetra::global_size_t>::invalid(),
                                          c_gids.view(0, count),
                                          coarseMap->getIndexBase(),
                                          stridingInfo_,
                                          coarseMap->getComm(),
                                          domainGIDOffset);
}

}  // namespace MueLu

#define MUELU_CLASSICALPFACTORY_SHORT
#endif  // MUELU_CLASSICALPFACTORY_DEF_HPP
