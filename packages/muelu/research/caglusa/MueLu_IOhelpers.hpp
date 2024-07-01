#ifndef MUELU_IOHELPERS_HPP
#define MUELU_IOHELPERS_HPP

#include <Teuchos_ParameterList.hpp>
#include <Xpetra_IO.hpp>

namespace MueLu {

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
struct IOhelpers {
  static Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Read(const std::string& filename,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rowMap,
       RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > colMap,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > domainMap = Teuchos::null,
       const RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > rangeMap  = Teuchos::null,
       const bool callFillComplete                                                = true,
       const bool binary                                                          = false,
       const bool readLocal                                                       = false) {
    using IO = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A;
    if (readLocal)
      A = IO::ReadLocal(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary);
    else
      A = IO::Read(filename, rowMap, colMap, domainMap, rangeMap, callFillComplete, binary);
    return A;
  }

  static Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Read(std::string& filename,
       RCP<const Teuchos::Comm<int> >& comm) {
    Teuchos::ParameterList hierarchicalParams;
    Teuchos::updateParametersFromXmlFileAndBroadcast(filename, Teuchos::Ptr<Teuchos::ParameterList>(&hierarchicalParams), *comm);
    auto op = Read(hierarchicalParams, comm);
    return op;
  }

  static Teuchos::RCP<Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Read(Teuchos::ParameterList& hierarchicalParams,
       RCP<const Teuchos::Comm<int> >& comm) {
    using HOp                 = Xpetra::HierarchicalOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
    using blocked_matrix_type = typename HOp::blocked_matrix_type;
    using blocked_map_type    = typename blocked_matrix_type::blocked_map_type;
    using matrix_type         = typename HOp::matrix_type;
    using map_type            = typename HOp::map_type;
    using lo_vec_type         = typename blocked_map_type::lo_vec_type;

    auto lib = Xpetra::UseTpetra;
    RCP<HOp> op;
    RCP<const map_type> map, near_colmap, clusterCoeffMap, ghosted_clusterCoeffMap, clusterMap, ghosted_clusterMap;
    RCP<matrix_type> nearField, basisMatrix, kernelApproximations, kernelBlockGraph;

    std::vector<RCP<blocked_matrix_type> > transferMatrices;
    RCP<lo_vec_type> clusterSizes;
    RCP<blocked_map_type> blockedClusterMap, ghosted_blockedClusterMap;
    RCP<blocked_matrix_type> blockKernelApproximations;

    const bool readBinary = hierarchicalParams.get<bool>("read binary", false);
    const bool readLocal  = hierarchicalParams.get<bool>("read local", false);

    using IO = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    // row, domain and range map of the operator
    map = IO::ReadMap(hierarchicalParams.get<std::string>("map"), lib, comm, readBinary);
    // colmap of near field
    near_colmap = IO::ReadMap(hierarchicalParams.get<std::string>("near colmap"), lib, comm, readBinary);
    if (hierarchicalParams.isType<std::string>("coefficient map")) {
      // 1-to-1 map for the cluster coefficients
      clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("coefficient map"), lib, comm, readBinary);
      // overlapping map for the cluster coefficients
      ghosted_clusterCoeffMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted coefficient map"), lib, comm, readBinary);
      // 1-to-1 map for the clusters
      clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("cluster map"), lib, comm, readBinary);
      // overlapping map for the clusters
      ghosted_clusterMap = IO::ReadMap(hierarchicalParams.get<std::string>("ghosted cluster map"), lib, comm, readBinary);

      // blocked cluster map
      clusterSizes      = Xpetra::IO<Scalar, LocalOrdinal, GlobalOrdinal, Node>::ReadMultiVectorLO(hierarchicalParams.get<std::string>("gid_cluster_to_gid_coeff"), clusterMap)->getVectorNonConst(0);
      blockedClusterMap = rcp(new blocked_map_type(clusterCoeffMap, clusterSizes));
    }

    // near field interactions
    nearField = Read(hierarchicalParams.get<std::string>("near field matrix"), map, near_colmap, map, map, true, readBinary, readLocal);

    if (hierarchicalParams.isType<std::string>("coefficient map")) {
      // far field basis expansion coefficients
      basisMatrix = IOhelpers::Read(hierarchicalParams.get<std::string>("basis expansion coefficient matrix"), map, clusterCoeffMap, clusterCoeffMap, map, true, readBinary, readLocal);

      // far field interactions
      kernelApproximations = IOhelpers::Read(hierarchicalParams.get<std::string>("far field interaction matrix"), clusterCoeffMap, ghosted_clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, true, readBinary, readLocal);
      // block graph of far field interactions
      kernelBlockGraph = IOhelpers::Read(hierarchicalParams.get<std::string>("far field interaction matrix") + ".block", clusterMap, ghosted_clusterMap, clusterMap, clusterMap, true, readBinary, readLocal);

      {
        auto import                           = kernelBlockGraph->getCrsGraph()->getImporter();
        RCP<lo_vec_type> ghosted_clusterSizes = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(ghosted_clusterMap);
        ghosted_clusterSizes->doImport(*clusterSizes, *import, Xpetra::INSERT);
        ghosted_blockedClusterMap = rcp(new blocked_map_type(ghosted_clusterCoeffMap, ghosted_clusterSizes));
      }

      blockKernelApproximations = rcp(new blocked_matrix_type(kernelApproximations, kernelBlockGraph, blockedClusterMap, ghosted_blockedClusterMap));

      // Transfer matrices
      auto transfersList = hierarchicalParams.sublist("shift coefficient matrices");
      for (int i = 0; i < transfersList.numParams(); i++) {
        std::string filename = transfersList.get<std::string>(std::to_string(i));
        auto transferPoint   = IOhelpers::Read(filename, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, clusterCoeffMap, true, readBinary, readLocal);
        auto transferBlock   = IOhelpers::Read(filename + ".block", clusterMap, clusterMap, clusterMap, clusterMap, true, readBinary, readLocal);
        auto transfer        = rcp(new blocked_matrix_type(transferPoint, transferBlock, blockedClusterMap));
        transferMatrices.push_back(transfer);
      }
    }

    RCP<Teuchos::ParameterList> params;
    if (hierarchicalParams.isSublist("params")) {
      params = rcp(new Teuchos::ParameterList(hierarchicalParams.sublist("params")));
    }

    if (hierarchicalParams.isType<std::string>("coefficient map")) {
      op = rcp(new HOp(nearField, blockKernelApproximations, basisMatrix, transferMatrices, params));

      return op;
    } else
      return nearField;
  }
};

}  // namespace MueLu

#endif
