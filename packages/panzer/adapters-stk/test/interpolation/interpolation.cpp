// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <string>
#include <iostream>
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_DefaultComm.hpp"
#include <Teuchos_RCP.hpp>
// #include "Panzer_Traits.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_Interpolation.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_STK_MeshFactory.hpp"
#include "Panzer_STK_LineMeshFactory.hpp"
#include "Panzer_STK_SquareQuadMeshFactory.hpp"
#include "Panzer_STK_SquareTriMeshFactory.hpp"
#include "Panzer_STK_CubeHexMeshFactory.hpp"
#include "Panzer_STK_CubeTetMeshFactory.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Teko_Utilities.hpp"

using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::RCP;
using Teuchos::rcpFromRef;

// #define PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE = 1

namespace panzer {

  Teuchos::RCP<panzer_stk::STK_Interface> getMesh(Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
                                                  const std::string &mesh_type) {

    int x_elements = 4;
    int y_elements = 2;
    int z_elements = 2;
    // build mesh
    Teuchos::RCP<panzer_stk::STK_MeshFactory> mesh_factory;
    if (mesh_type == "tet")
      mesh_factory = rcp(new panzer_stk::CubeTetMeshFactory());
    else if (mesh_type == "hex")
      mesh_factory = rcp(new panzer_stk::CubeHexMeshFactory());
    else if (mesh_type == "tri")
      mesh_factory = rcp(new panzer_stk::SquareTriMeshFactory());
    else if (mesh_type == "quad")
      mesh_factory = rcp(new panzer_stk::SquareQuadMeshFactory());
    else if (mesh_type == "line")
      mesh_factory = rcp(new panzer_stk::LineMeshFactory());
    else
      throw;
    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());
    pl->set<int>("X Elements", x_elements);
    if (mesh_type != "line")
      pl->set<int>("Y Elements",y_elements);
    if (mesh_type == "tet" || mesh_type == "hex")
      pl->set<int>("Z Elements",z_elements);
    if (mesh_type != "line") {
      pl->set<int>("X Procs", -1);
      pl->set<int>("Y Procs", -1);
    }
    if (mesh_type == "tet" || mesh_type == "hex")
      pl->set<int>("Z Procs",-1);
    mesh_factory->setParameterList(pl);
    auto mesh = mesh_factory->buildUncommitedMesh((*comm->getRawMpiComm())());
    mesh_factory->completeMeshConstruction(*mesh,(*comm->getRawMpiComm())());
    return mesh;
  }

  Teuchos::RCP<const panzer::FieldPattern> buildFieldPattern(const std::string basis_type, int basis_order, const shards::CellTopology & cell_topology)
  {

    // build a geometric pattern from a single basis
    RCP<Intrepid2::Basis<PHX::exec_space, double, double> > basis = panzer::createIntrepid2Basis<PHX::exec_space, double, double>(basis_type, basis_order, cell_topology);
    RCP<const panzer::FieldPattern> pattern = rcp(new panzer::Intrepid2FieldPattern(basis));
    return pattern;
  }

  Teuchos::RCP<const panzer::BlockedDOFManager> buildBlockedIndexer64(Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
                                                                      const std::string &mesh_type,
                                                                      const std::string &space,
                                                                      const std::vector<int>& orders)
  {
    Teuchos::RCP<panzer_stk::STK_Interface> mesh = getMesh(comm, mesh_type);
    const Teuchos::RCP<panzer::ConnManager> connManager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::BlockedDOFManager> indexer = rcp(new panzer::BlockedDOFManager());
    indexer->setConnManager(connManager, (*comm->getRawMpiComm())());
    indexer->setOrientationsRequired(true);

    std::vector<shards::CellTopology> topologies;
    connManager->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];

    const int numBlocks = orders.size();
    std::vector<std::vector<std::string> > fieldOrder(numBlocks);

    for (size_t j = 0; j<orders.size(); ++j) {
      int order = orders[j];
      Teuchos::RCP<const FieldPattern> pattern = buildFieldPattern(space, order, topology);
      indexer->addField(std::to_string(order), pattern);
      fieldOrder[j].push_back(std::to_string(order));
    }

    indexer->setFieldOrder(fieldOrder);
    indexer->buildGlobalUnknowns();

    return indexer;
  }

  Teuchos::RCP<const panzer::BlockedDOFManager> buildDeRhamBlockedIndexer64(Teuchos::RCP<const Teuchos::MpiComm<int> > &comm,
                                                                            const std::string &mesh_type,
                                                                            int order)
  {
    Teuchos::RCP<panzer_stk::STK_Interface> mesh = getMesh(comm, mesh_type);
    const Teuchos::RCP<panzer::ConnManager> connManager = rcp(new panzer_stk::STKConnManager(mesh));

    RCP<panzer::BlockedDOFManager> indexer = rcp(new panzer::BlockedDOFManager());
    indexer->setConnManager(connManager, (*comm->getRawMpiComm())());

    std::vector<shards::CellTopology> topologies;
    connManager->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];
    int dim = static_cast<int>(topology.getDimension());

    TEUCHOS_ASSERT(dim == 2 || dim == 3);

    const int numBlocks = dim + 1;
    std::vector<std::vector<std::string> > fieldOrder(numBlocks);

    int j = 0;

    // HGRAD
    Teuchos::RCP<const FieldPattern> patternHGrad = buildFieldPattern("HGrad", order, topology);
    indexer->addField("HGRAD", patternHGrad);
    fieldOrder[j++].push_back("HGRAD");

    // HCURL
    Teuchos::RCP<const FieldPattern> patternHCURL = buildFieldPattern("HCurl", order, topology);
    indexer->addField("HCURL", patternHCURL);
    fieldOrder[j++].push_back("HCURL");

    if (dim == 3) {
      // HDIV
      Teuchos::RCP<const FieldPattern> patternHDIV = buildFieldPattern("HDiv", order, topology);
      indexer->addField("HDIV", patternHDIV);
      fieldOrder[j++].push_back("HDIV");
    }

    // HVOL
    Teuchos::RCP<const FieldPattern> patternHVOL = buildFieldPattern("HVol", order-1, topology);
    indexer->addField("HVOL", patternHVOL);
    fieldOrder[j++].push_back("HVOL");

    indexer->setFieldOrder(fieldOrder);
    indexer->buildGlobalUnknowns();

    return indexer;
  }


  Teuchos::RCP<panzer::DOFManager> getDOFManager(RCP<const panzer::BlockedDOFManager> blockedDOFMngr, const std::string basisName) {
    using UGI = panzer::GlobalIndexer;
    std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
    int fieldNum = blockedDOFMngr->getFieldNum(basisName);
    int blockIndex = blockedDOFMngr->getFieldBlock(fieldNum);
    RCP<panzer::DOFManager> ugi = rcp_dynamic_cast<panzer::DOFManager>(fieldDOFMngrs[blockIndex], true);
    return ugi;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Construct FE spaces of different order and check that the product of
  // interpolations between them matches the interpolation between lowest
  // and highest order space.
  void checkInterpolations(const std::string &mesh_type,
                           const std::string &space,
                           const std::vector<int> &orders,
                           Teuchos::FancyOStream& out,
                           bool &success)
  {

    TEUCHOS_ASSERT(orders.size() >= 3);

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(comm, mesh_type, space, orders);

    RCP<const panzer::ConnManager> conn = blkIndexer->getConnManager();
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);

    std::vector<Teko::LinearOp> interpolations;

    for (size_t j = 0; j < orders.size() - 1; ++j) {
      out << "Building interpolation " << orders[j+1] << " " << orders[j] << std::endl;
      interpolations.push_back(buildInterpolation(conn,
                                                  getDOFManager(blkIndexer, std::to_string(orders[j])),
                                                  getDOFManager(blkIndexer, std::to_string(orders[j+1])),
                                                  std::to_string(orders[j]), std::to_string(orders[j+1]),
                                                  Intrepid2::OPERATOR_VALUE));
#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
      Teko::writeMatrix("interpolation_"+mesh_type+"_"+space+"_"+std::to_string(orders[j])+std::to_string(orders[j+1]), interpolations[j]);
#endif
    }

    out << "Building interpolation " << orders[orders.size()-1] << " " << orders[0] << std::endl;
    auto interp = buildInterpolation(conn,
                                     getDOFManager(blkIndexer, std::to_string(orders[0])),
                                     getDOFManager(blkIndexer, std::to_string(orders[orders.size()-1])),
                                     std::to_string(orders[0]),
                                     std::to_string(orders[orders.size()-1]),
                                     Intrepid2::OPERATOR_VALUE);
#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
    Teko::writeMatrix("interpolation_"+mesh_type+"_"+space+"_"+std::to_string(orders[0])+std::to_string(orders[orders.size()-1]), interp);
#endif

    Teko::LinearOp diff, temp;
    double diffNorm;

    temp = interpolations[0];
    for (size_t j = 1; j<interpolations.size(); ++j) {
      temp = Teko::explicitMultiply(interpolations[j], temp);
    }

    diff = Teko::explicitAdd(interp, Teko::scale(-1.0, temp));
    diffNorm = Teko::infNorm(diff);
    out << diffNorm << std::endl;
    TEST_ASSERT(diffNorm < 1.0e-8);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, line_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("line", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tri_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("tri", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tri_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("tri", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tri_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkInterpolations("tri", "HVol", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, quad_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("quad", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, quad_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("quad", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, quad_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkInterpolations("quad", "HVol", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tet_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("tet", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tet_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("tet", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tet_hdiv_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("tet", "HDiv", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, tet_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkInterpolations("tet", "HVol", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, hex_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("hex", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, hex_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("hex", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, hex_hdiv_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkInterpolations("hex", "HDiv", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tInterpolation, hex_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkInterpolations("hex", "HVol", orders, out, success);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Construct the discrete deRham complex and check that the range of
  // derivatives are in the kernel of the next derivative.
  void checkDerivatives(const int basis_order,
                        const std::string &mesh_type,
                        Teuchos::FancyOStream& out,
                        bool &success)
  {

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    RCP<const panzer::BlockedDOFManager> blkIndexer = buildDeRhamBlockedIndexer64(comm, mesh_type, basis_order);

    RCP<const panzer::ConnManager> conn = blkIndexer->getConnManager();
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];
    int dim = static_cast<int>(topology.getDimension());

    Teko::LinearOp grad;
    Teko::LinearOp curl;
    Teko::LinearOp div;

    out << "Build gradient\n";
    grad = buildInterpolation(conn,
                              getDOFManager(blkIndexer, "HGRAD"),
                              getDOFManager(blkIndexer, "HCURL"),
                              "HGRAD", "HCURL", Intrepid2::OPERATOR_GRAD);
#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
    Teko::writeMatrix("grad_"+mesh_type+"_"+std::to_string(basis_order), grad);
#endif
    if (dim == 2) {
      out << "Build curl\n";
      curl = buildInterpolation(conn,
                                getDOFManager(blkIndexer, "HCURL"),
                                getDOFManager(blkIndexer, "HVOL"),
                                "HCURL", "HVOL", Intrepid2::OPERATOR_CURL);
#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
      Teko::writeMatrix("curl_"+mesh_type+"_"+std::to_string(basis_order), curl);
#endif
    }
    else if (dim == 3) {
      out << "Build curl\n";
      curl = buildInterpolation(conn,
                                getDOFManager(blkIndexer, "HCURL"),
                                getDOFManager(blkIndexer, "HDIV"),
                                "HCURL", "HDIV", Intrepid2::OPERATOR_CURL);
#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
      Teko::writeMatrix("curl_"+mesh_type+"_"+std::to_string(basis_order), curl);
#endif
      out << "Build div\n";
      div = buildInterpolation(conn,
                               getDOFManager(blkIndexer, "HDIV"),
                               getDOFManager(blkIndexer, "HVOL"),
                               "HDIV", "HVOL", Intrepid2::OPERATOR_DIV);
#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
      Teko::writeMatrix("div_"+mesh_type+"_"+std::to_string(basis_order), div);
#endif
    }

    Teko::LinearOp diff;
    double diffNorm;

    diff = Teko::explicitMultiply(curl, grad);
    diffNorm = Teko::infNorm(diff);
    out << "|curl * grad| = " << diffNorm << std::endl;
    TEST_ASSERT(diffNorm < 1.0e-8);

    if (dim == 3) {
      diff = Teko::explicitMultiply(div, curl);
      diffNorm = Teko::infNorm(diff);
      out << "|div * curl| = " << diffNorm << std::endl;
      TEST_ASSERT(diffNorm < 1.0e-8);
    }

  }

  TEUCHOS_UNIT_TEST(tDerivatives, tri_order1)
  {
    checkDerivatives( /*basis_order=*/1, "tri", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, quad_order1)
  {
    checkDerivatives( /*basis_order=*/1, "quad", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, tet_order1)
  {
    checkDerivatives( /*basis_order=*/1, "tet", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, hex_order1)
  {
    checkDerivatives( /*basis_order=*/1, "hex", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, tri_order2)
  {
    checkDerivatives( /*basis_order=*/2, "tri", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, quad_order2)
  {
    checkDerivatives( /*basis_order=*/2, "quad", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, tet_order2)
  {
    checkDerivatives( /*basis_order=*/2, "tet", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, hex_order2)
  {
    checkDerivatives( /*basis_order=*/2, "hex", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, tri_order3)
  {
    checkDerivatives( /*basis_order=*/3, "tri", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, quad_order3)
  {
    checkDerivatives( /*basis_order=*/3, "quad", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, tet_order3)
  {
    checkDerivatives( /*basis_order=*/3, "tet", out, success);
  }

  TEUCHOS_UNIT_TEST(tDerivatives, hex_order3)
  {
    checkDerivatives( /*basis_order=*/3, "hex", out, success);
  }


  void diffOperators(const Teko::LinearOp &A,
                     const Teko::LinearOp &B,
                     Teuchos::FancyOStream& out,
                     bool &success) {
    using Scalar = double;
    using ThyLinOp = Thyra::TpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal>;
    using CrsMatrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>;
    using Operator = Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal>;
    using MV = Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>;

    const RCP<const ThyLinOp>  opA = rcp_dynamic_cast<const ThyLinOp>(A, true);
    const RCP<const CrsMatrix> tpA = rcp_dynamic_cast<const CrsMatrix>(opA->getConstTpetraOperator(),true);
    const RCP<const ThyLinOp>  opB = rcp_dynamic_cast<const ThyLinOp>(B, true);
    const RCP<const Operator> tpB = opB->getConstTpetraOperator();

    auto tp_domainmap = tpA->getDomainMap();
    auto tp_rangemap  = tpA->getRangeMap();

    auto testX  = rcp(new MV(tp_domainmap, 1));
    auto testY1 = rcp(new MV(tp_rangemap, 1));
    auto testY2 = rcp(new MV(tp_rangemap, 1));
    testX->randomize();
    testY1->putScalar(1.);
    testY2->putScalar(1.);

    double diffNorm;

    tpA->apply(*testX, *testY1, Teuchos::NO_TRANS, 3.0, 2.0);
    tpB->apply(*testX, *testY2, Teuchos::NO_TRANS, 3.0, 2.0);
    testY1->update(-1.0,*testY2,1.0);
    diffNorm = testY1->getVector(0)->norm2();
    out << "norm difference for 3 * M * X + 2 Y: " << diffNorm << std::endl;
    TEST_ASSERT(diffNorm < 1.0e-8);

    tpA->apply(*testX, *testY1);
    tpB->apply(*testX, *testY2);
    testY1->update(-1.0, *testY2, 1.0);
    diffNorm = testY1->getVector(0)->norm2();
    out << "norm difference for M * X: " << diffNorm << std::endl;
    TEST_ASSERT(diffNorm < 1.0e-8);

    testX  = rcp(new MV(tp_rangemap, 1));
    testY1 = rcp(new MV(tp_domainmap, 1));
    testY2 = rcp(new MV(tp_domainmap, 1));
    testX->randomize();
    testY1->putScalar(1.);
    testY2->putScalar(1.);

    tpA->apply(*testX, *testY1, Teuchos::TRANS, 3.0, 2.0);
    tpB->apply(*testX, *testY2, Teuchos::TRANS, 3.0, 2.0);
    testY1->update(-1.0, *testY2, 1.0);
    diffNorm = testY1->getVector(0)->norm2();
    out << "norm difference for 3 * M^T * X + 2 Y: " << diffNorm << std::endl;
    TEST_ASSERT(diffNorm < 1.0e-8);

    tpA->apply(*testX, *testY1, Teuchos::TRANS);
    tpB->apply(*testX, *testY2, Teuchos::TRANS);

#ifdef PANZER_INTERPOLATION_WRITE_MATRICES_TO_FILE
    static int counter = 0;
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> >::writeDenseFile("X_" + std::to_string(counter)+".mm", *testX);
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> >::writeDenseFile("Y1_" + std::to_string(counter)+".mm", *testY1);
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> >::writeDenseFile("Y2_" + std::to_string(counter)+".mm", *testY2);
    ++counter;
#endif

    testY1->update(-1.0, *testY2, 1.0);
    diffNorm = testY1->getVector(0)->norm2();
    out << "norm difference for M^T * X: " << diffNorm << std::endl;
    TEST_ASSERT(diffNorm < 1.0e-8);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Construct FE spaces of different order and check that assembled and
  // matrix-free operators match.
  void checkMatrixFreeInterpolations(const std::string &mesh_type,
                                     const std::string &space,
                                     const std::vector<int> &orders,
                                     Teuchos::FancyOStream& out,
                                     bool &success)
  {

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    RCP<const panzer::BlockedDOFManager> blkIndexer = buildBlockedIndexer64(comm, mesh_type, space, orders);

    RCP<const panzer::ConnManager> conn = blkIndexer->getConnManager();
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);

    const size_t worksetSize = 2000;
    const bool forceVectorial = false;
    const bool useTpetra = true;

    for (size_t j = 0; j < orders.size() - 1; ++j) {
      out << "Building interpolation " << orders[j+1] << " " << orders[j] << std::endl;
      auto interpolation = buildInterpolation(conn,
                                              getDOFManager(blkIndexer, std::to_string(orders[j])),
                                              getDOFManager(blkIndexer, std::to_string(orders[j+1])),
                                              std::to_string(orders[j]), std::to_string(orders[j+1]),
                                              Intrepid2::OPERATOR_VALUE,
                                              worksetSize,
                                              forceVectorial,
                                              useTpetra,
                                              /*matrixFree=*/false);
      auto matrixFreeinterpolation = buildInterpolation(conn,
                                                        getDOFManager(blkIndexer, std::to_string(orders[j])),
                                                        getDOFManager(blkIndexer, std::to_string(orders[j+1])),
                                                        std::to_string(orders[j]), std::to_string(orders[j+1]),
                                                        Intrepid2::OPERATOR_VALUE,
                                                        worksetSize,
                                                        forceVectorial,
                                                        useTpetra,
                                                        /*matrixFree=*/true);
      diffOperators(interpolation, matrixFreeinterpolation, out, success);
    }

    out << "Building interpolation " << orders[orders.size()-1] << " " << orders[0] << std::endl;
    auto interpolation = buildInterpolation(conn,
                                            getDOFManager(blkIndexer, std::to_string(orders[0])),
                                            getDOFManager(blkIndexer, std::to_string(orders[orders.size()-1])),
                                            std::to_string(orders[0]),
                                            std::to_string(orders[orders.size()-1]),
                                            Intrepid2::OPERATOR_VALUE,
                                            worksetSize,
                                            forceVectorial,
                                            useTpetra,
                                            /*matrixFree=*/false);
    auto matrixFreeinterpolation = buildInterpolation(conn,
                                                      getDOFManager(blkIndexer, std::to_string(orders[0])),
                                                      getDOFManager(blkIndexer, std::to_string(orders[orders.size()-1])),
                                                      std::to_string(orders[0]),
                                                      std::to_string(orders[orders.size()-1]),
                                                      Intrepid2::OPERATOR_VALUE,
                                                      worksetSize,
                                                      forceVectorial,
                                                      useTpetra,
                                                      /*matrixFree=*/true);
    diffOperators(interpolation, matrixFreeinterpolation, out, success);
  }


  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, line_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("line", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tri_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("tri", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tri_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("tri", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tri_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkMatrixFreeInterpolations("tri", "HVol", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, quad_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("quad", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, quad_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("quad", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, quad_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkMatrixFreeInterpolations("quad", "HVol", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tet_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("tet", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tet_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("tet", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tet_hdiv_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("tet", "HDiv", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, tet_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkMatrixFreeInterpolations("tet", "HVol", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, hex_hgrad_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("hex", "HGrad", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, hex_hcurl_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("hex", "HCurl", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, hex_hdiv_123)
  {
    std::vector<int> orders = {1, 2, 3};
    checkMatrixFreeInterpolations("hex", "HDiv", orders, out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeInterpolation, hex_hvol_012)
  {
    std::vector<int> orders = {0, 1, 2};
    checkMatrixFreeInterpolations("hex", "HVol", orders, out, success);
  }


  ////////////////////////////////////////////////////////////////////////////
  // Construct the discrete deRham complex and check that assembled and
  // matrix-free operators match.
  void checkMatrixFreeDerivatives(const int basis_order,
                                  const std::string &mesh_type,
                                  Teuchos::FancyOStream& out,
                                  bool &success)
  {

    Teuchos::RCP<const Teuchos::MpiComm<int> > comm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(Teuchos::DefaultComm<int>::getComm());

    RCP<const panzer::BlockedDOFManager> blkIndexer = buildDeRhamBlockedIndexer64(comm, mesh_type, basis_order);

    RCP<const panzer::ConnManager> conn = blkIndexer->getConnManager();
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];
    int dim = static_cast<int>(topology.getDimension());

    const size_t worksetSize = 2000;
    const bool forceVectorial = false;
    const bool useTpetra = true;

    out << "Build gradient\n";
    auto grad = buildInterpolation(conn,
                                   getDOFManager(blkIndexer, "HGRAD"),
                                   getDOFManager(blkIndexer, "HCURL"),
                                   "HGRAD", "HCURL", Intrepid2::OPERATOR_GRAD,
                                   worksetSize,
                                   forceVectorial,
                                   useTpetra,
                                   /*matrixFree=*/false);
    auto grad_mf = buildInterpolation(conn,
                                      getDOFManager(blkIndexer, "HGRAD"),
                                      getDOFManager(blkIndexer, "HCURL"),
                                      "HGRAD", "HCURL", Intrepid2::OPERATOR_GRAD,
                                      worksetSize,
                                      forceVectorial,
                                      useTpetra,
                                      /*matrixFree=*/true);
    diffOperators(grad, grad_mf, out, success);

    if (dim == 2) {
      out << "Build curl\n";
      auto curl = buildInterpolation(conn,
                                     getDOFManager(blkIndexer, "HCURL"),
                                     getDOFManager(blkIndexer, "HVOL"),
                                     "HCURL", "HVOL", Intrepid2::OPERATOR_CURL,
                                     worksetSize,
                                     forceVectorial,
                                     useTpetra,
                                     /*matrixFree=*/false);
      auto curl_mf = buildInterpolation(conn,
                                        getDOFManager(blkIndexer, "HCURL"),
                                        getDOFManager(blkIndexer, "HVOL"),
                                        "HCURL", "HVOL", Intrepid2::OPERATOR_CURL,
                                        worksetSize,
                                        forceVectorial,
                                        useTpetra,
                                        /*matrixFree=*/true);
      diffOperators(curl, curl_mf, out, success);
    }
    else if (dim == 3) {
      out << "Build curl\n";
      auto curl = buildInterpolation(conn,
                                     getDOFManager(blkIndexer, "HCURL"),
                                     getDOFManager(blkIndexer, "HDIV"),
                                     "HCURL", "HDIV", Intrepid2::OPERATOR_CURL,
                                     worksetSize,
                                     forceVectorial,
                                     useTpetra,
                                     /*matrixFree=*/false);
      auto curl_mf = buildInterpolation(conn,
                                        getDOFManager(blkIndexer, "HCURL"),
                                        getDOFManager(blkIndexer, "HDIV"),
                                        "HCURL", "HDIV", Intrepid2::OPERATOR_CURL,
                                        worksetSize,
                                        forceVectorial,
                                        useTpetra,
                                        /*matrixFree=*/true);
      diffOperators(curl, curl_mf, out, success);

      out << "Build div\n";
      auto div = buildInterpolation(conn,
                                    getDOFManager(blkIndexer, "HDIV"),
                                    getDOFManager(blkIndexer, "HVOL"),
                                    "HDIV", "HVOL", Intrepid2::OPERATOR_DIV,
                                    worksetSize,
                                    forceVectorial,
                                    useTpetra,
                                    /*matrixFree=*/false);
      auto div_mf = buildInterpolation(conn,
                                       getDOFManager(blkIndexer, "HDIV"),
                                       getDOFManager(blkIndexer, "HVOL"),
                                       "HDIV", "HVOL", Intrepid2::OPERATOR_DIV,
                                       worksetSize,
                                       forceVectorial,
                                       useTpetra,
                                       /*matrixFree=*/true);
      diffOperators(div, div_mf, out, success);
    }
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, tri_order1)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/1, "tri", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, quad_order1)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/1, "quad", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, tet_order1)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/1, "tet", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, hex_order1)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/1, "hex", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, tri_order2)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/2, "tri", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, quad_order2)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/2, "quad", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, tet_order2)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/2, "tet", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, hex_order2)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/2, "hex", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, tri_order3)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/3, "tri", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, quad_order3)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/3, "quad", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, tet_order3)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/3, "tet", out, success);
  }

  TEUCHOS_UNIT_TEST(tMatrixFreeDerivatives, hex_order3)
  {
    checkMatrixFreeDerivatives( /*basis_order=*/3, "hex", out, success);
  }

}
