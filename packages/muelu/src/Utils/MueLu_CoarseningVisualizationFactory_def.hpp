// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COARSENINGVISUALIZATIONFACTORY_DEF_HPP_
#define MUELU_COARSENINGVISUALIZATIONFACTORY_DEF_HPP_

#include <Xpetra_Matrix.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include "MueLu_CoarseningVisualizationFactory_decl.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = VisualizationHelpers::GetValidParameterList();

  validParamList->set<int>("visualization: start level", 0, "visualize only levels with level ids greater or equal than start level");  // Remove me?

  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Prolongator factory. The user has to declare either P or Ptent but not both at the same time.");
  validParamList->set<RCP<const FactoryBase> >("Ptent", Teuchos::null, "Tentative prolongator factory. The user has to declare either P or Ptent as input but not both at the same time");
  validParamList->set<RCP<const FactoryBase> >("Coordinates", Teuchos::null, "Factory for Coordinates.");
  validParamList->set<RCP<const FactoryBase> >("Graph", Teuchos::null, "Factory for Graph.");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  this->Input(fineLevel, "Coordinates");

  const ParameterList &pL = this->GetParameterList();
  TEUCHOS_TEST_FOR_EXCEPTION(pL.isParameter("P") && GetFactory("P") != Teuchos::null &&
                                 pL.isParameter("Ptent") && GetFactory("Ptent") != Teuchos::null,
                             Exceptions::RuntimeError,
                             "You must not declare both P and Ptent. Use only once for visualization.");
  TEUCHOS_TEST_FOR_EXCEPTION(GetFactory("P") == Teuchos::null && GetFactory("Ptent") == Teuchos::null, Exceptions::RuntimeError,
                             "You have to either declare P or Ptent for visualization, but not both.");

  if (GetFactory("P") != Teuchos::null && GetFactory("Ptent") == Teuchos::null)
    this->Input(coarseLevel, "P");
  else if (GetFactory("Ptent") != Teuchos::null && GetFactory("P") == Teuchos::null)
    this->Input(coarseLevel, "Ptent");

  if (pL.get<bool>("visualization: fine graph edges"))
    Input(fineLevel, "Graph");
#if 0
    if(pL.get<bool>("visualization: coarse graph edges")) {
      Input(coarseLevel, "Coordinates");
      Input(coarseLevel, "Graph");
    }
#endif
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void CoarseningVisualizationFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &fineLevel, Level &coarseLevel) const {
  RCP<LWGraph> fineGraph  = Teuchos::null;
  RCP<Matrix> P           = Teuchos::null;
  const ParameterList &pL = this->GetParameterList();
  if (this->GetFactory("P") != Teuchos::null && this->GetFactory("Ptent") == Teuchos::null)
    P = Get<RCP<Matrix> >(coarseLevel, "P");
  if (GetFactory("Ptent") != Teuchos::null && GetFactory("P") == Teuchos::null)
    P = Get<RCP<Matrix> >(coarseLevel, "Ptent");

  RCP<const Teuchos::Comm<int> > comm = P->getRowMap()->getComm();

  LocalOrdinal dofsPerNode        = 1;
  LocalOrdinal stridedRowOffset   = 0;
  RCP<const StridedMap> strRowMap = Teuchos::null;
  if (P->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(P->getRowMap("stridedMaps")) != Teuchos::null) {
    strRowMap            = Teuchos::rcp_dynamic_cast<const StridedMap>(P->getRowMap("stridedMaps"));
    LocalOrdinal blockid = strRowMap->getStridedBlockId();
    if (blockid > -1) {
      std::vector<size_t> stridingInfo = strRowMap->getStridingData();
      for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
        stridedRowOffset += stridingInfo[j];
      dofsPerNode = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
    } else {
      dofsPerNode = strRowMap->getFixedBlockSize();
    }
    GetOStream(Runtime1) << "CoarseningVisualizationFactory::Build():"
                         << " #dofs per node = " << dofsPerNode << std::endl;
  }

  LocalOrdinal columnsPerNode        = dofsPerNode;
  LocalOrdinal stridedColumnOffset   = 0;
  RCP<const StridedMap> strDomainMap = Teuchos::null;
  if (P->IsView("stridedMaps") && Teuchos::rcp_dynamic_cast<const StridedMap>(P->getRowMap("stridedMaps")) != Teuchos::null) {
    strDomainMap         = Teuchos::rcp_dynamic_cast<const StridedMap>(P->getColMap("stridedMaps"));
    LocalOrdinal blockid = strDomainMap->getStridedBlockId();

    if (blockid > -1) {
      std::vector<size_t> stridingInfo = strDomainMap->getStridingData();
      for (size_t j = 0; j < Teuchos::as<size_t>(blockid); j++)
        stridedColumnOffset += stridingInfo[j];
      columnsPerNode = Teuchos::as<LocalOrdinal>(stridingInfo[blockid]);
    } else {
      columnsPerNode = strDomainMap->getFixedBlockSize();
    }
    GetOStream(Runtime1) << "CoarseningVisualizationFactory::Build():"
                         << " #columns per node = " << columnsPerNode << std::endl;
  }

  // TODO add support for overlapping aggregates
  TEUCHOS_TEST_FOR_EXCEPTION(strDomainMap->getLocalNumElements() != P->getColMap()->getLocalNumElements(), Exceptions::RuntimeError,
                             "CoarseningVisualization only supports non-overlapping transfers");

  // number of local "aggregates"
  LocalOrdinal numLocalAggs = strDomainMap->getLocalNumElements() / columnsPerNode;
  std::vector<std::set<LocalOrdinal> > localAggs(numLocalAggs);

  // do loop over all local rows of prolongator and extract column information
  for (LO row = 0; row < Teuchos::as<LO>(P->getRowMap()->getLocalNumElements()); ++row) {
    ArrayView<const LO> indices;
    ArrayView<const SC> vals;
    P->getLocalRowView(row, indices, vals);

    for (typename ArrayView<const LO>::iterator c = indices.begin(); c != indices.end(); ++c) {
      localAggs[(*c) / columnsPerNode].insert(row / dofsPerNode);
    }
  }

  // determine number of "aggs" per proc and calculate local "agg" offset...
  std::vector<int> myLocalAggsPerProc(comm->getSize(), 0);
  std::vector<int> numLocalAggsPerProc(comm->getSize(), 0);
  myLocalAggsPerProc[comm->getRank()] = numLocalAggs;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, comm->getSize(), &myLocalAggsPerProc[0], &numLocalAggsPerProc[0]);

  LocalOrdinal myAggOffset = 0;
  for (int i = 0; i < comm->getRank(); ++i) {
    myAggOffset += numLocalAggsPerProc[i];
  }

  /*for (LocalOrdinal i = 0; i < numLocalAggs; ++i) {

    std::cout << "PROC: " << comm->getRank() << " Local aggregate: " << i + myAggOffset << " with nodes: ";
    for( typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
      std::cout << *it << ", ";
    }
    std::cout << std::endl;
  }*/

  // get fine level coordinate information
  Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > coords = Get<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > >(fineLevel, "Coordinates");

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::as<LO>(P->getRowMap()->getLocalNumElements()) / dofsPerNode != Teuchos::as<LocalOrdinal>(coords->getLocalLength()), Exceptions::RuntimeError,
                             "Number of fine level nodes in coordinates is inconsistent with dof based information");

  // communicate fine level coordinates if necessary
  if (pL.get<bool>("visualization: fine graph edges")) {
    fineGraph = Get<RCP<LWGraph> >(fineLevel, "Graph");

    RCP<Import> coordImporter                                                                                                          = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(coords->getMap(), fineGraph->GetImportMap());
    RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > ghostedCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>::Build(fineGraph->GetImportMap(), coords->getNumVectors());
    ghostedCoords->doImport(*coords, *coordImporter, Xpetra::INSERT);
    coords = ghostedCoords;
  }

  Teuchos::RCP<const Map> nodeMap = coords->getMap();

  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> xCoords = Teuchos::arcp_reinterpret_cast<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>(coords->getData(0));
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> yCoords = Teuchos::arcp_reinterpret_cast<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>(coords->getData(1));
  Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> zCoords = Teuchos::null;
  if (coords->getNumVectors() == 3) {
    zCoords = Teuchos::arcp_reinterpret_cast<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>(coords->getData(2));
  }

  // determine number of nodes on fine level
  LocalOrdinal numFineNodes = Teuchos::as<LocalOrdinal>(coords->getLocalLength());

  // create vertex2AggId array
  ArrayRCP<LocalOrdinal> vertex2AggId(numFineNodes, -1);
  for (LocalOrdinal i = 0; i < numLocalAggs; ++i) {
    // TODO: check if entry = -1
    for (typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
      vertex2AggId[*it] = i;
    }
  }

  // we have no information which node is root and which is not
  // we could either look at the entries in P again or build some new root nodes
  // assuming that root nodes are usually at the center of the aggregate
  std::vector<bool> isRoot(numFineNodes, false);
  for (LocalOrdinal i = 0; i < numLocalAggs; ++i) {
    typename Teuchos::ScalarTraits<Scalar>::coordinateType xCenter = 0.0;
    typename Teuchos::ScalarTraits<Scalar>::coordinateType yCenter = 0.0;
    typename Teuchos::ScalarTraits<Scalar>::coordinateType zCenter = 0.0;

    // loop over all nodes in aggregate i and determine center coordinates of aggregate
    for (typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
      xCenter += xCoords[*it];
      yCenter += yCoords[*it];
      if (coords->getNumVectors() == 3) zCenter += zCoords[*it];
    }
    xCenter /= Teuchos::as<LocalOrdinal>(localAggs[i].size());
    yCenter /= Teuchos::as<LocalOrdinal>(localAggs[i].size());
    zCenter /= Teuchos::as<LocalOrdinal>(localAggs[i].size());

    // loop over all nodes in aggregate i and find node which is closest to aggregate center
    LocalOrdinal rootCandidate                                         = -1;
    typename Teuchos::ScalarTraits<Scalar>::coordinateType minDistance = Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::coordinateType>::one() / Teuchos::ScalarTraits<typename Teuchos::ScalarTraits<Scalar>::coordinateType>::sfmin();
    for (typename std::set<LocalOrdinal>::iterator it = localAggs[i].begin(); it != localAggs[i].end(); ++it) {
      typename Teuchos::ScalarTraits<Scalar>::coordinateType tempx = xCenter - xCoords[*it];
      typename Teuchos::ScalarTraits<Scalar>::coordinateType tempy = yCenter - yCoords[*it];
      typename Teuchos::ScalarTraits<Scalar>::coordinateType tempz = 0.0;
      if (coords->getNumVectors() == 3) tempz = zCenter - zCoords[*it];
      typename Teuchos::ScalarTraits<Scalar>::coordinateType mydistance = 0.0;
      mydistance += tempx * tempx;
      mydistance += tempy * tempy;
      mydistance += tempz * tempz;
      mydistance = sqrt(mydistance);
      if (mydistance <= minDistance) {
        minDistance   = mydistance;
        rootCandidate = *it;
      }
    }

    isRoot[rootCandidate] = true;
  }

  std::vector<LocalOrdinal> vertices;
  std::vector<LocalOrdinal> geomSize;
  int vizLevel = pL.get<int>("visualization: start level");
  if (vizLevel <= fineLevel.GetLevelID()) {
    std::string aggStyle = pL.get<std::string>("visualization: style");
    if (aggStyle == "Point Cloud")
      this->doPointCloud(vertices, geomSize, numLocalAggs, numFineNodes);
    else if (aggStyle == "Jacks")
      this->doJacks(vertices, geomSize, numLocalAggs, numFineNodes, isRoot, vertex2AggId);
    else if (aggStyle == "Convex Hulls") {
      // TODO do a smarter distinction and check the number of z levels...
      // loop over all coordinates and collect coordinate components in sets...
      if (coords->getNumVectors() == 3)
        this->doConvexHulls3D(vertices, geomSize, numLocalAggs, numFineNodes, isRoot, vertex2AggId, xCoords, yCoords, zCoords);
      else if (coords->getNumVectors() == 2)
        this->doConvexHulls2D(vertices, geomSize, numLocalAggs, numFineNodes, isRoot, vertex2AggId, xCoords, yCoords);
    } else {
      GetOStream(Warnings0) << "   Warning: Unrecognized agg style.\nPossible values are Point Cloud, Jacks, Convex Hulls.\nDefaulting to Point Cloud." << std::endl;
      aggStyle = "Point Cloud";
      this->doPointCloud(vertices, geomSize, numLocalAggs, numFineNodes);
    }
  }

  // write out fine edge information
  if (pL.get<bool>("visualization: fine graph edges")) {
    TEUCHOS_TEST_FOR_EXCEPTION(fineGraph == Teuchos::null, Exceptions::RuntimeError,
                               "Could not get information about fine graph.");

    std::vector<LocalOrdinal> fine_edges_vertices;
    std::vector<LocalOrdinal> fine_edges_geomSize;
    this->doGraphEdges(fine_edges_vertices, fine_edges_geomSize, fineGraph, xCoords, yCoords, zCoords);

    std::string fEdgeFineFile = this->getFileName(comm->getSize(), comm->getRank(), fineLevel.GetLevelID(), pL);
    std::string fEdgeFile     = fEdgeFineFile.insert(fEdgeFineFile.rfind(".vtu"), "-finegraph");
    std::ofstream edgeStream(fEdgeFile.c_str());

    std::vector<int> uniqueFineEdges = this->makeUnique(fine_edges_vertices);
    this->writeFileVTKOpening(edgeStream, uniqueFineEdges, fine_edges_geomSize);
    this->writeFileVTKNodes(edgeStream, uniqueFineEdges, nodeMap);
    this->writeFileVTKData(edgeStream, uniqueFineEdges, myAggOffset, vertex2AggId, comm->getRank());
    this->writeFileVTKCoordinates(edgeStream, uniqueFineEdges, xCoords, yCoords, zCoords, coords->getNumVectors());
    this->writeFileVTKCells(edgeStream, uniqueFineEdges, fine_edges_vertices, fine_edges_geomSize);
    this->writeFileVTKClosing(edgeStream);
    edgeStream.close();
  }

  // communicate fine level coordinates if necessary
#if 0  // we don't have access to the coarse graph
    if (pL.get<bool>("visualization: coarse graph edges")) {
      RCP<LWGraph> coarseGraph = Get<RCP<LWGraph> >(coarseLevel, "Graph");

      Teuchos::RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > coarsecoords = Get<RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > >(coarseLevel, "Coordinates");

      RCP<Import> coarsecoordImporter = Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(coarsecoords->getMap(), coarseGraph->GetImportMap());
      RCP<Xpetra::MultiVector<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node> > coarseghostedCoords = Xpetra::MultiVectorFactory<typename Teuchos::ScalarTraits<Scalar>::coordinateType, LocalOrdinal, GlobalOrdinal, Node>::Build(coarseGraph->GetImportMap(), coarsecoords->getNumVectors());
      coarseghostedCoords->doImport(*coarsecoords, *coarsecoordImporter, Xpetra::INSERT);
      coarsecoords = coarseghostedCoords;

      Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> cx = Teuchos::arcp_reinterpret_cast<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>(coarsecoords->getData(0));
      Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> cy = Teuchos::arcp_reinterpret_cast<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>(coarsecoords->getData(1));
      Teuchos::ArrayRCP<const typename Teuchos::ScalarTraits<Scalar>::coordinateType> cz = Teuchos::null;
      if(coarsecoords->getNumVectors() == 3) {
        cz = Teuchos::arcp_reinterpret_cast<const typename Teuchos::ScalarTraits<Scalar>::coordinateType>(coarsecoords->getData(2));
      }

      Teuchos::RCP<const Map> coarsenodeMap = coarsecoords->getMap();

      std::vector<LocalOrdinal> coarse_edges_vertices;
      std::vector<LocalOrdinal> coarse_edges_geomSize;
      this->doGraphEdges(coarse_edges_vertices, coarse_edges_geomSize, coarseGraph, cx, cy, cz);

      std::string cEdgeFineFile = this->getFileName(comm->getSize(), comm->getRank(), coarseLevel.GetLevelID(), pL);
      std::string cEdgeFile = cEdgeFineFile.insert(cEdgeFineFile.rfind(".vtu"), "-coarsegraph");
      std::ofstream cedgeStream(cEdgeFile.c_str());

      std::vector<int> uniqueCoarseEdges = this->makeUnique(coarse_edges_vertices);
      this->writeFileVTKOpening(cedgeStream, uniqueCoarseEdges, coarse_edges_geomSize);
      this->writeFileVTKNodes(cedgeStream, uniqueCoarseEdges, coarsenodeMap);
      //this->writeFileVTKData(edgeStream, uniqueCoarseEdges, myAggOffset, vertex2AggId, comm->getRank());
      this->writeFileVTKCoordinates(cedgeStream, uniqueCoarseEdges, cx, cy, cz, coarsecoords->getNumVectors());
      this->writeFileVTKCells(cedgeStream, uniqueCoarseEdges, coarse_edges_vertices, coarse_edges_geomSize);
      this->writeFileVTKClosing(cedgeStream);
      cedgeStream.close();
    }
#endif

  if (pL.get<int>("visualization: start level") <= fineLevel.GetLevelID()) {
    // write out coarsening information
    std::string filenameToWrite = this->getFileName(comm->getSize(), comm->getRank(), fineLevel.GetLevelID(), pL);
    std::ofstream fout(filenameToWrite.c_str());

    std::vector<int> uniqueFine = this->makeUnique(vertices);
    this->writeFileVTKOpening(fout, uniqueFine, geomSize);
    this->writeFileVTKNodes(fout, uniqueFine, nodeMap);
    this->writeFileVTKData(fout, uniqueFine, myAggOffset, vertex2AggId, comm->getRank());
    this->writeFileVTKCoordinates(fout, uniqueFine, xCoords, yCoords, zCoords, coords->getNumVectors());
    this->writeFileVTKCells(fout, uniqueFine, vertices, geomSize);
    this->writeFileVTKClosing(fout);
    fout.close();

    // create pvtu file
    if (comm->getRank() == 0) {
      std::string pvtuFilename = this->getPVTUFileName(comm->getSize(), comm->getRank(), fineLevel.GetLevelID(), pL);
      std::string baseFname    = this->getBaseFileName(comm->getSize(), fineLevel.GetLevelID(), pL);
      std::ofstream pvtu(pvtuFilename.c_str());
      this->writePVTU(pvtu, baseFname, comm->getSize(), pL.get<bool>("visualization: fine graph edges"));
      pvtu.close();
    }
  }

  if (comm->getRank() == 0 && pL.get<bool>("visualization: build colormap")) {
    this->buildColormap();
  }
}
}  // namespace MueLu

#endif /* MUELU_AGGREGATIONEXPORTFACTORY_DEF_HPP_ */
