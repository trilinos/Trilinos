#ifndef XPETRA_REGION_MATRIX_IMPL_HPP_
#define XPETRA_REGION_MATRIX_IMPL_HPP_

// Xpetra
#include "Xpetra_EpetraCrsMatrix_fwd.hpp"
#include "Xpetra_IO.hpp"
#include "Xpetra_RegionMatrix_decl.hpp"
#include "Xpetra_RegionUtils_impl.hpp"

template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::RegionMatrix(
    const std::string& matrixFileName,
    Teuchos::RCP<const RegionManager> regionManager,
    Teuchos::RCP<const Teuchos::Comm<int> > comm)
    : comm_(comm),
      regionManager_(regionManager),
//      regionMatricesInitialized_(0),
      compositeMatrix_(Teuchos::null),
      regionMatrices_(0)
{
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  if (comm_->getRank() == 0)
    *out << "Constructing 'RegionMatrix'." << std::endl;

//  // Reset initialization flags
//  regionMatricesInitialized_.clear();
//  for (int i = 0; i < regionManager->getNumRegions(); ++i)
//    regionMatricesInitialized_.push_back(false);

  // Setup matrices in composite and region view
  setupCompositeMatrix(matrixFileName);
  setupRegionMatrices();

  return;
}

template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
void Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::setupCompositeMatrix(const std::string& matrixFileName)
{
  if (comm_->getRank() == 0)
    std::cout << "Starting construction of Composite Map" << std::endl;

  // Create Xpetra map for composite stiffness matrix
  Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > xpetraMap;
  xpetraMap = Xpetra::MapFactory<LO,GO,NO>::Build(lib, regionManager_->getNumNodes(), regionManager_->getCompositeRowMap(), 0, comm_);

  if (comm_->getRank() == 0)
    std::cout << "Finished construction of Composite Map" << std::endl;

  if (comm_->getRank() == 0)
    std::cout << "Started reading composite matrix" << std::endl;

  // Import matrix from an .mm file into an Xpetra wrapper for an Epetra matrix
  compositeMatrix_ = Xpetra::IO<SC,LO,GO,NO>::Read(matrixFileName, regionManager_->getCompositeMap());

  if (comm_->getRank() == 0)
    std::cout << "Finished reading composite matrix" << std::endl;

  return;
}

template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
void Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::printCompositeMatrix(
    Teuchos::FancyOStream& out, Teuchos::EVerbosityLevel verbosity) const
{
  compositeMatrix_->describe(out, verbosity);

  return;
}

template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
void Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::setupRegionMatrices()
{
  // Enlarge composite matrix to duplicate the interface
  Teuchos::RCP<TpetraCrsMatrix> tpetraGlobalMatrix = MueLu::Utilities<SC,LO,GO,NO>::Op2NonConstTpetraCrs(compositeMatrix_);
  Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix =
      Teuchos::rcp(new Ifpack2::OverlappingRowMatrix<TpetraRowMatrix>(tpetraGlobalMatrix, 1));

//  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
//  if (comm_->getRank() == 0)
//    *out << std::endl << "enlargedMatrix:" << std::endl;
//  enlargedMatrix->describe(*out, Teuchos::VERB_EXTREME);

  // initialize region matrices as empty matrices
  regionMatrices_.clear();
  for (int i = 0; i < regionManager_->getNumRegions(); ++i)
  {
    Teuchos::RCP<const Map> xpetraMap = regionManager_->getRegionMap(i);
    const GO numElements = regionManager_->getNumNodesPerRegion(i);

    Teuchos::RCP<CrsMatrix> crsMatrix;
    if (lib == Xpetra::UseEpetra)
      crsMatrix = Teuchos::rcp(new EpetraCrsMatrix(xpetraMap, numElements));
    else if (lib == Xpetra::UseTpetra)
      crsMatrix = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<SC,LO,GO,NO>(xpetraMap, numElements));
    else
      std::cerr << "The library to build matrices must be either Epetra or Tpetra\n";

    Teuchos::RCP<Matrix> matrixPointer = rcp(new CrsMatrixWrap(crsMatrix));
    regionMatrices_.push_back(matrixPointer);
  }

  // extract region matrices from the enlargedMatrix
  for (GO i = 0; i < regionManager_->getNumRegions(); ++i)
  {
    Teuchos::RCP<Matrix> regionMatrix = regionMatrices_[i];
    regionMatrix->resumeFill();

    // Perform collapse or split operation when chopping the composite matrix into regional pieces
    switch (splitMethod)
    {
    case Xpetra::region_collapse:
    {
//      extractRegionByCollapsing(i, regionMatrix, enlargedMatrix);
      break;
    }
    case Xpetra::region_split:
    {
      extractRegionBySplitting(i, regionMatrix, enlargedMatrix);
      break;
    }
    default:
    {
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "The chosen splitting strategy is not available.");
      break;
    }
    }

    regionMatrix->fillComplete();
  }

  return;
}

//template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
//void Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::initializeRegionMatrices(
//    const GO regionID, Teuchos::RCP<Matrix>& regionMatrix,
//    Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix)
//{
////  //Region matrices are initially built to be a chopped version of the composite matrix
////  TEUCHOS_TEST_FOR_EXCEPTION(regionMatricesInitialized_[regionID], Exceptions::RuntimeError,
////      "Surrogate region stiffness matrices are already initialized by chopping the composite stiffness matrix.\n");
//
//  Teuchos::Array<std::tuple<GO,GO> > regionToAll = regionManager_->getRegionToAll(regionID);
//
//  //THIS IS THE CORE OF THE PROBLEM WHERE ONE NEEDS TO POPULATE THE REGIONAL MATRICES BY ACCESSING ENTRIES OF THE GLOBAL MATRIX
//  Teuchos::ArrayView<const GO> myRegionElements = regionMatrix->getRowMap()->getNodeElementList();
//  for (typename ArrayView<const GO>::iterator iter = myRegionElements.begin();
//      iter != myRegionElements.end(); ++iter)
//  {
//    //Nodes are saved in data structures with 1 as base index
//    Xpetra::CheckerRegionToAll<GO> unaryPredicate(*iter+1);
//    typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator;
//    composite_iterator = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//        Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
//
//    TEUCHOS_TEST_FOR_EXCEPTION(composite_iterator == regionToAll.end(),
//        Exceptions::RuntimeError,
//        "Process ID: " << comm_->getRank() << " - Region: " << regionID
//        << " - "<<" node with region index: " << *iter+1
//        << " is not in regionToAll[" << regionID << "]\n");
//
//    GO nodeID = std::get<1>(*composite_iterator);
//    LO nodeLocalID = enlargedMatrix->getRowMap()->getLocalElement(nodeID-1);
//
//    Teuchos::ArrayView<const GO> inds;
//    Teuchos::ArrayView<const SC> vals;
//    enlargedMatrix->getLocalRowView(nodeLocalID, inds, vals);
//
//    std::vector<GO> regionIndsVector(0);
//    std::vector<SC> regionValsVector(0);
//
//    for (LO i = 0; i < inds.size(); ++i)
//    {
//      //Nodes are saved in data structures with 1 as base index
//      GO composite_col_ind = enlargedMatrix->getColMap()->getGlobalElement(inds[i]) + 1;
//      Xpetra::CheckerAllToRegion<GO> unaryPredicate2(composite_col_ind);
//      typename Teuchos::Array<std::tuple<GO,GO> >::iterator region_iterator;
//      region_iterator = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//          Xpetra::CheckerAllToRegion<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicate2);
//      if (region_iterator != regionToAll.end())
//      {
//        regionIndsVector.push_back(std::get<0>(*region_iterator) - 1);
//        regionValsVector.push_back(vals[i]);
//      }
//    }
//
//    Teuchos::ArrayView<GO> regionInds(regionIndsVector);
//    Teuchos::ArrayView<SC> regionVals(regionValsVector);
//    regionMatrix->insertGlobalValues(*iter, regionInds, regionVals);
//  }
//
//  regionMatricesInitialized_[regionID] = true;
//
//  return;
//}

//template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
//void Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::regionCollapse(
//    GO regionID, Teuchos::RCP<Matrix>& regionMatrix,
//    Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix)
//{
//  TEUCHOS_TEST_FOR_EXCEPTION(!regionMatricesInitialized_[regionID], Exceptions::RuntimeError,
//      "The composite stiffness matrix must be chopped into surrogate region matrices before collapsing.\n");
//  TEUCHOS_TEST_FOR_EXCEPTION(regionManager_->getNumNodesPerRegion(regionID) != Teuchos::as<GO>(regionMatrices_[regionID]->getGlobalNumRows()), Exceptions::RuntimeError,
//      "Process ID: " << comm_->getRank() << " - Number of region nodes in region "
//      << regionID+1 << " does not coincide with the value returned by regionMatrices_["
//      << regionID+1 << "]->getGlobalNumRows().\n");
//
//  Teuchos::Array<std::tuple<GO,GO> > regionToAll = regionManager_->getRegionToAll(regionID);
//
//  //This portion of the code assumes that the number of region nodes is the same on each direction of the domain
//  //For a 2D problem we have then nx = ny = sqrt(num_region_nodes_)
//  GO n;
//  GO nx;
//  GO ny;
//
//  n = regionManager_->getNumNodesPerRegion(regionID);
//  nx = std::sqrt(n);
//  ny = nx;
//
//  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<double>(nx - std::floor(static_cast<double>(std::sqrt(static_cast<double>(n))))) != 0.0, Exceptions::RuntimeError,
//      "The code assumes that the regions are 2D and that the number of region nodes is the same on each direction of the domain \n");
//
//  //interfaceNodesToRegions contains nodes on an interface between any regions
//  Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > interfaceNodesToRegions = regionManager_->getInterfaceNodesToRegions();
//
//  Teuchos::ArrayView<const GO> MyRegionElements = regionMatrix->getRowMap()->getNodeElementList();
//  for (typename ArrayView<const GO>::iterator iter = MyRegionElements.begin();
//      iter != MyRegionElements.end(); ++iter)
//  {
//    // Nodes are saved in data structures with 1 as base index
//    GO region_node_idx = *iter+1;
//    Xpetra::CheckerRegionToAll<GO> unaryPredicate(region_node_idx);
//    typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator;
//    composite_iterator = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//        Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
//    TEUCHOS_TEST_FOR_EXCEPTION(composite_iterator == regionToAll.end(), Exceptions::RuntimeError,
//        "Process ID: " << comm_->getRank() << " - Region: " << regionID << " - "
//        << " node with region index: " << region_node_idx << " is not in regionToAll["
//        << regionID << "].\n");
//
//    GO composite_node_idx = std::get<1>(*composite_iterator);
//    Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2(composite_node_idx);
//    typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator;
//    interface_iterator = std::find_if<typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator,
//        Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodesToRegions.begin(), interfaceNodesToRegions.end(), unaryPredicate2);
//
//    // Here we assuming that a specific labeling choice is adopted region wise and we use it to distinguish coarse node from fine nodes
//    bool coarse_point = false;
//    if (region_node_idx % 3 == 1) // ToDo (mayr.mt) Make coarsening factor an input parameter.
//      coarse_point = true;
//
//    GO region_node_idx_neighbor1 = 0;
//    GO region_node_idx_neighbor2 = 0;
//
//    // Horizontal-Vertical Collapse
//    if (interface_iterator != interfaceNodesToRegions.end()
//        and region_node_idx > ny && region_node_idx <= (nx - 1) * ny
//        and !coarse_point)
//    {
//      region_node_idx_neighbor1 = region_node_idx - ny;
//      region_node_idx_neighbor2 = region_node_idx + ny;
//    }
//    else if(interface_iterator != interfaceNodesToRegions.end() && region_node_idx % ny > 1 && !coarse_point )
//    {
//      region_node_idx_neighbor1 = region_node_idx - 1;
//      region_node_idx_neighbor2 = region_node_idx + 1;
//    }
//
//    if (region_node_idx_neighbor1 != 0 && region_node_idx_neighbor2 != 0)
//    {
//      // Computation of composite index for neighbor1 node
//      Xpetra::CheckerRegionToAll<GO> unaryPredicateLeft(region_node_idx_neighbor1);
//      typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator_neighbor1;
//      composite_iterator_neighbor1 = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//          Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateLeft);
//
//      // Computation of composite index for neighbor2 node
//      Xpetra::CheckerRegionToAll<GO> unaryPredicateRight(region_node_idx_neighbor2);
//      typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator_neighbor2;
//      composite_iterator_neighbor2 = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//          Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateRight);
//
//      TEUCHOS_TEST_FOR_EXCEPTION(
//          composite_iterator_neighbor1 == regionToAll.end() || composite_iterator_neighbor2 == regionToAll.end(), Exceptions::RuntimeError,
//          "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<" node with region index: "<<region_node_idx<<" lies on the interface between regions: "<<std::get<1>(*interface_iterator)<<" BUT has compositely mislabeled neighbouring nodes missing from regionToAll \n");
//
//      // Check to see if neighbor1 node lies on a coarse line
//      GO composite_node_idx_neighbor1 = std::get<1>(*composite_iterator_neighbor1);
//      Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2neighbor1(composite_node_idx_neighbor1);
//      typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator_neighbor1;
//      interface_iterator_neighbor1 = std::find_if<typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator,
//          Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodesToRegions.begin(), interfaceNodesToRegions.end(), unaryPredicate2neighbor1);
//
//      // Check to see if neighbor2 node lies on a coarse line
//      GO composite_node_idx_neighbor2 = std::get<1>(*composite_iterator_neighbor2);
//      Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2neighbor2(composite_node_idx_neighbor2);
//      typename Teuchos::Array< std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator_neighbor2;
//      interface_iterator_neighbor2 = std::find_if<typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator,
//          Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodesToRegions.begin(), interfaceNodesToRegions.end(), unaryPredicate2neighbor2);
//
//      //I apply the collapse only if the current node is a fine node which lies on a coarse line
//      //This means that the neighbor1 node and neighbor2 node must both lie on the coarse line as well
//      if( interface_iterator_neighbor1!=interfaceNodesToRegions.end() && interface_iterator_neighbor2!=interfaceNodesToRegions.end() )
//      {
//
//        //For each fine node on a horizontal coarse line on the interface, I extract the rows from the composite matrix
//        LO node_idx = enlargedMatrix->getRowMap()->getLocalElement(composite_node_idx-1);
//        LO node_idx_neighbor1 = enlargedMatrix->getRowMap()->getLocalElement(composite_node_idx_neighbor1-1);
//        LO node_idx_neighbor2 = enlargedMatrix->getRowMap()->getLocalElement(composite_node_idx_neighbor2-1);
//        Teuchos::ArrayView<const LO> inds;
//        Teuchos::ArrayView<const SC> vals;
//        enlargedMatrix->getLocalRowView( node_idx, inds, vals );
//        Teuchos::ArrayView<const LO> inds_neighbor1;
//        Teuchos::ArrayView<const SC> vals_neighbor1;
//        enlargedMatrix->getLocalRowView( node_idx_neighbor1, inds_neighbor1, vals_neighbor1 );
//        Teuchos::ArrayView<const LO> inds_neighbor2;
//        Teuchos::ArrayView<const SC> vals_neighbor2;
//        enlargedMatrix->getLocalRowView( node_idx_neighbor2, inds_neighbor2, vals_neighbor2 );
//
//        std::vector<LO> inds_vector = createVector(inds);
//        std::vector<LO> inds_neighbor1_vector = createVector(inds_neighbor1);
//        std::vector<LO> inds_neighbor2_vector = createVector(inds_neighbor2);
//
//        std::vector<GO> composite_inds_vector(0);
//        std::vector<GO> composite_inds_neighbor1_vector(0);
//        std::vector<GO> composite_inds_neighbor2_vector(0);
//
//        for( typename std::vector<LO>::iterator iter_node = inds_vector.begin(); iter_node!=inds_vector.end(); ++iter_node )
//          composite_inds_vector.push_back( enlargedMatrix->getRowMap()->getGlobalElement(*iter_node) );
//        std::sort( composite_inds_vector.begin(), composite_inds_vector.end() );
//
//        for( typename std::vector<LO>::iterator iter_node = inds_neighbor1_vector.begin(); iter_node!=inds_neighbor1_vector.end(); ++iter_node )
//          composite_inds_neighbor1_vector.push_back( enlargedMatrix->getRowMap()->getGlobalElement(*iter_node) );
//
//        std::sort( composite_inds_neighbor1_vector.begin(), composite_inds_neighbor1_vector.end() );
//
//        for( typename std::vector<LO>::iterator iter_node = inds_neighbor2_vector.begin(); iter_node!=inds_neighbor2_vector.end(); ++iter_node )
//          composite_inds_neighbor2_vector.push_back( enlargedMatrix->getRowMap()->getGlobalElement(*iter_node) );
//
//        std::sort( composite_inds_neighbor2_vector.begin(), composite_inds_neighbor2_vector.end() );
//
//        //IDENTIFICATION OF EXTERNAL NODES THROUGH COMPOSITE INDICES STARTS HERE
//        std::vector<GO> composite_node_idx_neighbor1_extra;
//        std::vector<GO> composite_node_idx_neighbor2_extra;
//        std::vector<GO> composite_node_idx_extra(composite_inds_vector);
//
//        //The following triple of vectors is expected to EVENTUALLY contain only one entry:
//        //the label of the external node with information to collapse close to neighbor1, neighbor2 and central node
//        std::vector<GO> diff_neighbor1;
//        std::vector<GO> diff_neighbor2;
//        std::vector<GO> diff_center;
//
//        //Identification of external node from the side of neighbor1
//        {
//          //Compute the intersection between neoghbourhood of neighbor1 node and neighbourhood of central node
//          std::set_intersection(composite_inds_vector.begin(), composite_inds_vector.end(), composite_inds_neighbor1_vector.begin(), composite_inds_neighbor1_vector.end(), std::back_inserter(composite_node_idx_neighbor1_extra));
//          for( typename std::vector<GO>::iterator iter_node = composite_node_idx_neighbor1_extra.begin(); iter_node!=composite_node_idx_neighbor1_extra.end(); ++iter_node )
//          {
//            Xpetra::CheckerAllToRegion<GO> unaryPredicateExtra(*iter_node+1);
//            typename Array< std::tuple<GO, GO > >::iterator region_iterator_extra;
//            region_iterator_extra = std::find_if<typename Array< std::tuple< GO,GO > >::iterator,
//                Xpetra::CheckerAllToRegion<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateExtra);
//
//            //Invalidation of node indices for nodes belonging to the current region
//            if( region_iterator_extra!=regionToAll.end() )
//              *iter_node = -1;
//          }
//
//          //Removal of invalidated indices associated with nodes belonging to current region: (external nodes do not belong to this region)
//          composite_node_idx_neighbor1_extra.erase(std::remove(composite_node_idx_neighbor1_extra.begin(), composite_node_idx_neighbor1_extra.end(), -1),composite_node_idx_neighbor1_extra.end());
//
//          //External node from neighbor1 side does not belong to the neighborhood of neighbor2
//          std::set_difference(composite_node_idx_neighbor1_extra.begin(), composite_node_idx_neighbor1_extra.end(), composite_inds_neighbor2_vector.begin(), composite_inds_neighbor2_vector.end(), std::inserter(diff_neighbor1, diff_neighbor1.begin()));
//          TEUCHOS_TEST_FOR_EXCEPTION( diff_neighbor1.size()!=1 , Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<"Mislabeling of nodes obstructed the identification of the extra node: region node "<< region_node_idx<<" leads to diff_neighbor1.size()= "<<diff_neighbor1.size()<<" \n");
//        }
//
//        //Identification of external node from the side of neighbor2
//        {
//          //Compute the intersection between neighbourhood of neighbor2 node and neighbourhood of central node
//          std::set_intersection(composite_inds_vector.begin(), composite_inds_vector.end(), composite_inds_neighbor2_vector.begin(), composite_inds_neighbor2_vector.end(), std::back_inserter(composite_node_idx_neighbor2_extra));
//
//          for( typename std::vector<GO>::iterator iter_node = composite_node_idx_neighbor2_extra.begin(); iter_node!=composite_node_idx_neighbor2_extra.end(); ++iter_node )
//          {
//            Xpetra::CheckerAllToRegion<GO> unaryPredicateExtra(*iter_node+1);
//            typename Array< std::tuple<GO, GO > >::iterator region_iterator_extra;
//            region_iterator_extra = std::find_if<typename Array< std::tuple< GO,GO > >::iterator,
//                Xpetra::CheckerAllToRegion<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateExtra);
//
//            //Invalidation of node indices for nodes belonging to the current region
//            if (region_iterator_extra != regionToAll.end())
//              *iter_node = -1;
//          }
//
//          //Removal of invalidated indices associated with nodes belonging to current region: (external nodes do not belong to this region)
//          composite_node_idx_neighbor2_extra.erase(std::remove(composite_node_idx_neighbor2_extra.begin(), composite_node_idx_neighbor2_extra.end(), -1),composite_node_idx_neighbor2_extra.end());
//
//          //External node from neighbor2 side does not belong to the neighborhood of neighbor1
//          std::set_difference(composite_node_idx_neighbor2_extra.begin(), composite_node_idx_neighbor2_extra.end(), composite_inds_neighbor1_vector.begin(), composite_inds_neighbor1_vector.end(), std::inserter(diff_neighbor2, diff_neighbor2.begin()));
//          TEUCHOS_TEST_FOR_EXCEPTION( diff_neighbor2.size()!=1 , Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<"Mislabeling of nodes obstructed the identification of the extra node: region node "<< region_node_idx<<" leads to diff_neighbor2.size()= "<<diff_neighbor2.size()<<" \n");
//        }
//
//        //Identification of external node from the side of central node
//        {
//          for( typename std::vector<GO>::iterator iter_node = composite_node_idx_extra.begin(); iter_node!=composite_node_idx_extra.end(); ++iter_node )
//          {
//            Xpetra::CheckerAllToRegion<GO> unaryPredicateExtra(*iter_node+1);
//            typename Array< std::tuple<GO, GO > >::iterator region_iterator_extra;
//            region_iterator_extra = std::find_if<typename Array< std::tuple< GO,GO > >::iterator,
//                Xpetra::CheckerAllToRegion<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateExtra);
//
//            //Invalidation of node indices for nodes belonging to the current region
//            if( region_iterator_extra!=regionToAll.end() )
//              *iter_node = -1;
//          }
//
//          //Removal of invalidated indices associated with nodes belonging to current region: (external nodes do not belong to this region)
//          composite_node_idx_extra.erase(std::remove(composite_node_idx_extra.begin(), composite_node_idx_extra.end(), -1),composite_node_idx_extra.end());
//          std::vector<GO> diff_center_temp;
//
//          //At thie point composite_node_idx_extra contains indices of all the three external nodes: two of these must be removed since they are already tracked
//          //External nodes from neighbors1's and neighbor2's side must be removed
//          std::set_difference(composite_node_idx_extra.begin(), composite_node_idx_extra.end(), diff_neighbor1.begin(), diff_neighbor1.end(), std::inserter(diff_center_temp, diff_center_temp.begin()));
//          std::set_difference(diff_center_temp.begin(), diff_center_temp.end(), diff_neighbor2.begin(), diff_neighbor2.end(), std::inserter(diff_center, diff_center.begin()));
//          TEUCHOS_TEST_FOR_EXCEPTION( diff_center.size()!=1 , Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<"Mislabeling of nodes obstructed the identification of the extra node: region node "<< region_node_idx<<" leads to diff_center.size()= "<<diff_center.size()<<" \n");
//        }
//
//        //Computation of local indices for central node and its neighbors
//        LO local_region_node_idx = regionMatrix->getRowMap()->getLocalElement( region_node_idx );
//        LO local_region_node_idx_neighbor1 = regionMatrix->getRowMap()->getLocalElement(region_node_idx_neighbor1);
//        LO local_region_node_idx_neighbor2 = regionMatrix->getRowMap()->getLocalElement(region_node_idx_neighbor2);
//
//        //Computation of local indices for external nodes
//        LO local_extra_central = enlargedMatrix->getRowMap()->getLocalElement(diff_center[0]);
//        LO local_extra_neighbor1 = enlargedMatrix->getRowMap()->getLocalElement(diff_neighbor1[0]);
//        LO local_extra_neighbor2 = enlargedMatrix->getRowMap()->getLocalElement(diff_neighbor2[0]);
//
//        Teuchos::ArrayView<const GO> region_row;
//        Teuchos::ArrayView<const GO> region_col;
//        Teuchos::ArrayView<const SC> region_val;
//
//        //Extract Row view of the region matrix before collapsing
//        if (regionMatrix->isLocallyIndexed())
//          regionMatrix->getLocalRowView(local_region_node_idx, region_col, region_val);
//        else
//          regionMatrix->getGlobalRowView(region_node_idx, region_col, region_val);
//
//        //Extract Row of overlapped composite matrix to detect node with information to collapse
//        Teuchos::ArrayView<const GO> external_row;
//        Teuchos::ArrayView<const GO> external_col;
//        Teuchos::ArrayView<const SC> external_val;
//        enlargedMatrix->getLocalRowView(node_idx, external_col, external_val);
//
//        //neighbor1 collapse
//        {
//          SC initial_value = 0.0;
//          for (typename Teuchos::ArrayView<const GO>::iterator iter_view = region_col.begin(); iter_view != region_col.end(); ++iter_view)
//          {
//            if( regionMatrix -> isLocallyIndexed() )
//              if( *iter_view==local_region_node_idx_neighbor1 )
//                initial_value = region_val[iter_view - region_col.begin()];
//            if( regionMatrix -> isGloballyIndexed() )
//              if( *iter_view==region_node_idx_neighbor1 )
//                initial_value = region_val[iter_view - region_col.begin()];
//
//            if (initial_value != 0.0)
//              break;
//          }
//
//          SC external_value = 0.0;
//          for (typename Teuchos::ArrayView<const GO>::iterator iter_view = external_col.begin(); iter_view != external_col.end(); ++iter_view)
//          {
//            if (*iter_view == local_extra_neighbor1)
//              external_value = external_val[iter_view - external_col.begin()];
//
//            if (external_value != 0.0)
//              break;
//          }
//
//          SC new_value = external_value;// new matrix entry generated with the collapsing
//          std::vector<GO> new_entry_ind;
//          std::vector<SC> new_entry_val;
//          new_entry_ind.push_back(region_node_idx_neighbor1-1);
//          new_entry_val.push_back(new_value);
//
//          //If a nonzero value is already stored in the specified position, the new values is SUMMED to the already existing one
//          //See description of insertGlobalValues(...)
//          regionMatrix -> insertGlobalValues( region_node_idx-1, new_entry_ind, new_entry_val );
//        }
//        //neighbor2 collapse
//        {
//          SC initial_value = 0.0;
//          for (typename Teuchos::ArrayView<const GO>::iterator iter_view = region_col.begin(); iter_view != region_col.end(); ++iter_view)
//          {
//            if( regionMatrix -> isLocallyIndexed() )
//              if( *iter_view==local_region_node_idx_neighbor2 )
//                initial_value = region_val[iter_view - region_col.begin()];
//            if( regionMatrix -> isGloballyIndexed() )
//              if( *iter_view==region_node_idx_neighbor2 )
//                initial_value = region_val[iter_view - region_col.begin()];
//
//            if (initial_value != 0.0)
//              break;
//          }
//
//          SC external_value = 0.0;
//          for (typename Teuchos::ArrayView<const GO>::iterator iter_view = external_col.begin(); iter_view != external_col.end(); ++iter_view)
//          {
//            if( *iter_view==local_extra_neighbor2 )
//              external_value = external_val[iter_view - external_col.begin()];
//
//            if (external_value != 0.0)
//              break;
//          }
//
//          SC new_value = external_value;// new matrix entry generated with the collapsing
//          std::vector<GO> new_entry_ind;
//          std::vector<SC> new_entry_val;
//          new_entry_ind.push_back(region_node_idx_neighbor2-1);
//          new_entry_val.push_back(new_value);
//
//          //If a nonzero value is already stored in the specified position, the new values is SUMMED to the already existing one
//          //See description of insertGlobalValues(...)
//          regionMatrix -> insertGlobalValues( region_node_idx-1, new_entry_ind, new_entry_val );
//        }
//        //central node collapse
//        {
//          SC initial_value = 0.0;
//          for( typename ArrayView<const GO>::iterator iter_view = region_col.begin(); iter_view!=region_col.end(); ++iter_view  )
//          {
//            if( regionMatrix -> isLocallyIndexed() )
//              if( *iter_view==local_region_node_idx )
//                initial_value = region_val[iter_view - region_col.begin()];
//            if( regionMatrix -> isGloballyIndexed() )
//              if( *iter_view==region_node_idx )
//                initial_value = region_val[iter_view - region_col.begin()];
//
//            if (initial_value != 0.0)
//              break;
//          }
//
//          SC external_value = 0.0;
//          for( typename ArrayView<const GO>::iterator iter_view = external_col.begin(); iter_view!=external_col.end(); ++iter_view  )
//          {
//            if( *iter_view==local_extra_central )
//              external_value = external_val[iter_view - external_col.begin()];
//
//            if (external_value != 0.0)
//              break;
//          }
//
//          SC new_value = external_value;// new matrix entry generated with the collapsing
//          std::vector<GO> new_entry_ind;
//          std::vector<SC> new_entry_val;
//          new_entry_ind.push_back(region_node_idx-1);
//          new_entry_val.push_back(new_value);
//
//          //If a nonzero value is already stored in the specified position, the new values is SUMMED to the already existing one
//          //See description of insertGlobalValues(...)
//          regionMatrix->insertGlobalValues(region_node_idx - 1, new_entry_ind, new_entry_val);
//        }
//      }
//    }
//  }
//
//  return;
//}

template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
void Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::extractRegionBySplitting(
    GO regionID, Teuchos::RCP<Matrix>& regionMatrix,
    Teuchos::RCP<Ifpack2::OverlappingRowMatrix<TpetraRowMatrix> > enlargedMatrix)
{
//  TEUCHOS_TEST_FOR_EXCEPTION(!regionMatricesInitialized_[regionID], Exceptions::RuntimeError, "The composite stiffness matrix must be chopped into surrogate region matrices before collapsing.\n");
  TEUCHOS_TEST_FOR_EXCEPTION(regionManager_->getNumNodesPerRegion(regionID) != Teuchos::as<GO>(regionMatrices_[regionID]->getGlobalNumRows()), Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - Number of region nodes in region "<<regionID+1<<" does not coincide with the value returned by regionMatrix["<<regionID+1<<"]->getGlobalNumRows().\n");

  Teuchos::Array<std::tuple<GO,GO> > regionToAll = regionManager_->getRegionToAll(regionID);

  //This portion of the code assumes that the number of region nodes is the same on each direction of the domain
  //For a 2D problem we have then nx = ny = sqrt( num_region_ndoes_ )
  GO n;
  GO nx;
  GO ny;

  n = regionManager_->getNumNodesPerRegion(regionID);
  nx = std::sqrt(n);
  ny = nx;

  TEUCHOS_TEST_FOR_EXCEPTION(static_cast<double>(nx - std::floor(static_cast<double>(std::sqrt(static_cast<double>(n))))) != 0.0, Exceptions::RuntimeError,
      "The code assumes that the regions are 2D and that the number of region nodes is the same on each direction of the domain \n");

  //interfaceNodes contains nodes on an interface between any regions
  Teuchos::Array<std::tuple<int,Teuchos::Array<GO> > > interfaceNodes = regionManager_->getInterfaceNodesToRegions();

  Teuchos::ArrayView<const GO> MyRegionElements = regionMatrix->getRowMap()->getNodeElementList();
  for (typename Teuchos::ArrayView<const GO>::iterator iter = MyRegionElements.begin();
      iter != MyRegionElements.end(); ++iter)
  {
    GO region_node_idx = *iter;
    Xpetra::CheckerRegionToAll<GO> unaryPredicate(region_node_idx);
    typename Teuchos::Array< std::tuple<GO, GO > >::iterator composite_iterator;
    composite_iterator = std::find_if<typename Teuchos::Array< std::tuple< GO,GO> >::iterator,
        Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicate);
    TEUCHOS_TEST_FOR_EXCEPTION( composite_iterator==regionToAll.end(), Exceptions::RuntimeError, "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<" node with region index: "<<region_node_idx<<" is not in regionToAll["<<regionID<<"]"<<"\n" );

    GO composite_node_idx = std::get<1>( *composite_iterator );
    Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2( composite_node_idx );
    typename Teuchos::Array< std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator;
    interface_iterator = std::find_if<typename Teuchos::Array< std::tuple<GO,Teuchos::Array<GO> > >::iterator,
        Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2);

    std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;

    GO region_node_idx_neighbor_w = 0;
    GO region_node_idx_neighbor_e = 0;
    GO region_node_idx_neighbor_n = 0;
    GO region_node_idx_neighbor_s = 0;

    int count_neighbours = 0;

    if (interface_iterator != interfaceNodes.end() && region_node_idx + 1 > ny)
    {
      region_node_idx_neighbor_w = region_node_idx + 1 - ny;
      count_neighbours++;
    }

    if (interface_iterator != interfaceNodes.end()
        and region_node_idx + 1 <= (nx - 1) * ny)
    {
      region_node_idx_neighbor_e = region_node_idx + 1 + ny;
      count_neighbours++;
    }

    if (interface_iterator != interfaceNodes.end() && (region_node_idx + 1) % ny != 1)
    {
      region_node_idx_neighbor_s = region_node_idx + 1 - 1;
      count_neighbours++;
    }

    if (interface_iterator != interfaceNodes.end() && (region_node_idx + 1) % ny != 0)
    {
      region_node_idx_neighbor_n = region_node_idx + 1 + 1;
      count_neighbours++;
    }

    bool interface_line = false;
    bool interface_corner = false;

    if (count_neighbours == 3)
      interface_line = true;
    else if (count_neighbours == 2)
      interface_corner = true;

    typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator_neighbor_e;
    typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator_neighbor_e;
    typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator_neighbor_w;
    typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator_neighbor_w;
    typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator_neighbor_s;
    typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator_neighbor_s;
    typename Teuchos::Array<std::tuple<GO,GO> >::iterator composite_iterator_neighbor_n;
    typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator interface_iterator_neighbor_n;


//    if (interface_line || interface_corner)
//    {
//      // Computation of composite index for East node
//      if (region_node_idx_neighbor_e != 0)
//      {
//        Xpetra::CheckerRegionToAll<GO> unaryPredicateEast(region_node_idx_neighbor_e);
//        composite_iterator_neighbor_e = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//            Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateEast);
//
//        //Check to see if neighbor_e node lies on a coarse line
//        GO composite_node_idx_neighbor_e = std::get<1>( *composite_iterator_neighbor_e );
//        Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2neighborEast( composite_node_idx_neighbor_e );
//        interface_iterator_neighbor_e = std::find_if<typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator,
//            Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborEast);
//      }
//      else
//        interface_iterator_neighbor_e = interfaceNodes.end();
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Computation of composite index for West node
//      if (region_node_idx_neighbor_w != 0)
//      {
//        Xpetra::CheckerRegionToAll<GO> unaryPredicateWest(region_node_idx_neighbor_w);
//        composite_iterator_neighbor_w = std::find_if<typename Array<std::tuple<GO,GO> >::iterator, Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateWest);
//
//        //Check to see if neighbor_w node lies on a coarse line
//        GO composite_node_idx_neighbor_w = std::get<1>( *composite_iterator_neighbor_w );
//        Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2neighborWest( composite_node_idx_neighbor_w );
//        interface_iterator_neighbor_w = std::find_if<typename Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator, Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborWest);
//      }
//      else
//        interface_iterator_neighbor_w = interfaceNodes.end();
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Computation of composite index for South node
//      if (region_node_idx_neighbor_s != 0)
//      {
//        Xpetra::CheckerRegionToAll<GO> unaryPredicateSouth(region_node_idx_neighbor_s);
//        composite_iterator_neighbor_s = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//            Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateSouth);
//
//        //Check to see if neighbor_s node lies on a coarse line
//        GO composite_node_idx_neighbor_s = std::get<1>( *composite_iterator_neighbor_s );
//        Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2neighborSouth( composite_node_idx_neighbor_s );
//        interface_iterator_neighbor_s = std::find_if<typename Teuchos::Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator,
//            Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborSouth);
//      }
//      else
//        interface_iterator_neighbor_s = interfaceNodes.end();
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Computation of composite index for North node
//      if (region_node_idx_neighbor_n != 0)
//      {
//        Xpetra::CheckerRegionToAll<GO> unaryPredicateNorth(region_node_idx_neighbor_n);
//        composite_iterator_neighbor_n = std::find_if<typename Teuchos::Array<std::tuple<GO,GO> >::iterator,
//            Xpetra::CheckerRegionToAll<GO> >(regionToAll.begin(), regionToAll.end(), unaryPredicateNorth);
//
//        //Check to see if neighbor_n node lies on a coarse line
//        GO composite_node_idx_neighbor_n = std::get<1>( *composite_iterator_neighbor_n );
//        Xpetra::CheckerInterfaceNodes<GO> unaryPredicate2neighborNorth( composite_node_idx_neighbor_n );
//        interface_iterator_neighbor_n = std::find_if<typename Array<std::tuple<GO,Teuchos::Array<GO> > >::iterator,
//            Xpetra::CheckerInterfaceNodes<GO> >(interfaceNodes.begin(), interfaceNodes.end(), unaryPredicate2neighborNorth);
//      }
//      else
//        interface_iterator_neighbor_n = interfaceNodes.end();
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      int count_neighbours_interface = 0;
//      if (interface_iterator_neighbor_e != interfaceNodes.end())
//        count_neighbours_interface++;
//      if (interface_iterator_neighbor_w != interfaceNodes.end())
//        count_neighbours_interface++;
//      if (interface_iterator_neighbor_s != interfaceNodes.end())
//        count_neighbours_interface++;
//      if (interface_iterator_neighbor_n != interfaceNodes.end())
//        count_neighbours_interface++;
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      TEUCHOS_TEST_FOR_EXCEPTION(count_neighbours_interface > count_neighbours,
//          Exceptions::RuntimeError,
//          "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<" node with region index: "<<region_node_idx<<" has inconsistent information on the number of neighbours: count_neighbours = "<<count_neighbours<<"but count_neighbours_interface ="<<count_neighbours_interface<<"\n");
//
//      //First the splitting is applied on extra-diagonal entries
//
//      //Computation of local indices for central node and its neighbors
//      //Node index base start from 1 in the structures used, but Trilinos maps start from 0, so
//      //indices must be shifted by 1
//      LO local_region_node_idx = regionMatrix->getRowMap()->getLocalElement( region_node_idx-1 );
//      LO local_region_node_idx_neighbor_e = regionMatrix->getRowMap()->getLocalElement(region_node_idx_neighbor_e-1);
//      LO local_region_node_idx_neighbor_w = regionMatrix->getRowMap()->getLocalElement(region_node_idx_neighbor_w-1);
//      LO local_region_node_idx_neighbor_s = regionMatrix->getRowMap()->getLocalElement(region_node_idx_neighbor_s-1);
//      LO local_region_node_idx_neighbor_n = regionMatrix->getRowMap()->getLocalElement(region_node_idx_neighbor_n-1);
//
//      Teuchos::ArrayView<const GO> region_col;
//      Teuchos::ArrayView<const SC> region_val;
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Extract Row view of the region matrix
//      if (regionMatrix->isLocallyIndexed())
//        regionMatrix->getLocalRowView(local_region_node_idx, region_col, region_val);
//      else
//        regionMatrix->getGlobalRowView(*iter, region_col, region_val);
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      std::vector<LO> region_col_vector = createVector(region_col);
//      std::vector<LO> ind_vector(0);
//      std::vector<SC> val_vector(0);
//
//      // Extraction of the info about East neighbour to halve the associated entry in the matrix
//      if (interface_iterator_neighbor_e != interfaceNodes.end())
//      {
//        typename std::vector<GO>::iterator iter_east_vector;
//        GO east_ind;
//        if (regionMatrix->isLocallyIndexed())
//        {
//          iter_east_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_e);
//          east_ind = regionMatrix->getRowMap()->getGlobalElement(*iter_east_vector);
//        }
//        else
//        {
//          iter_east_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_e - 1);
//          east_ind = *iter_east_vector;
//        }
//        SC east_val = -0.5
//            * region_val[iter_east_vector - region_col_vector.begin()];
//        ind_vector.push_back(east_ind);
//        val_vector.push_back(east_val);
//      }
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Extraction of the info about West neighbour to halve the associated entry in the matrix
//      if (interface_iterator_neighbor_w != interfaceNodes.end())
//      {
//        typename std::vector<GO>::iterator iter_west_vector;
//        GO west_ind;
//        if (regionMatrix->isLocallyIndexed())
//        {
//          iter_west_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_w);
//          west_ind = regionMatrix->getRowMap()->getGlobalElement(*iter_west_vector);
//        }
//        else
//        {
//          iter_west_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_w - 1);
//          west_ind = *iter_west_vector;
//        }
//        SC west_val = -0.5 * region_val[iter_west_vector - region_col_vector.begin()];
//        ind_vector.push_back(west_ind);
//        val_vector.push_back(west_val);
//      }
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Extraction of the info about South neighbour to halve the associated entry in the matrix
//      if (interface_iterator_neighbor_s != interfaceNodes.end())
//      {
//        typename std::vector<GO>::iterator iter_south_vector;
//        GO south_ind;
//        if (regionMatrix->isLocallyIndexed())
//        {
//          iter_south_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_s);
//          south_ind = regionMatrix->getRowMap()->getGlobalElement(*iter_south_vector);
//        }
//        else
//        {
//          iter_south_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_s - 1);
//          south_ind = *iter_south_vector;
//        }
//        SC south_val = -0.5 * region_val[iter_south_vector - region_col_vector.begin()];
//        ind_vector.push_back(south_ind);
//        val_vector.push_back(south_val);
//      }
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Extraction of the info about North neighbour to halve the associated entry in the matrix
//      if (interface_iterator_neighbor_n != interfaceNodes.end())
//      {
//        typename std::vector<GO>::iterator iter_north_vector;
//        GO north_ind;
//        if (regionMatrix->isLocallyIndexed())
//        {
//          iter_north_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx_neighbor_n);
//          north_ind = regionMatrix->getRowMap()->getGlobalElement(*iter_north_vector);
//        }
//        else
//        {
//          iter_north_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx_neighbor_n - 1);
//          north_ind = *iter_north_vector;
//        }
//        SC north_val = -0.5 * region_val[iter_north_vector - region_col_vector.begin()];
//        ind_vector.push_back(north_ind);
//        val_vector.push_back(north_val);
//      }
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      // Extraction of the info about my Node ID to split the associated entry in the matrix
//      // The ratio used for the splitting depends on the number of regions this current node
//      // belongs to
//      typename std::vector<LO>::iterator iter_center_vector;
//      GO center_ind;
//      if (regionMatrix->isLocallyIndexed())
//      {
//        iter_center_vector = std::find(region_col_vector.begin(), region_col_vector.end(), local_region_node_idx);
//        center_ind = regionMatrix->getRowMap()->getGlobalElement(*iter_center_vector);
//      }
//      else
//      {
//        iter_center_vector = std::find(region_col_vector.begin(), region_col_vector.end(), region_node_idx - 1);
//        center_ind = *iter_center_vector;
//      }
//
//      //Count of the number of regions the current node belongs to
//      GO region_belonging = std::get<1>(*interface_iterator).size();
//      TEUCHOS_TEST_FOR_EXCEPTION(region_belonging < 2, Exceptions::RuntimeError,
//          "Process ID: "<<comm_->getRank()<<" - Region: "<<regionID<<" - "<<" node with composite index: "<<std::get<0>(*interface_iterator)<<" should lie on an interface between regions but the nubmer of regions it belongs to is only "<<region_belonging<<"\n" );
//
//      //If a node is on a corner between four interfaces, then each the entry A(node_idx,node_idx) must be split in four parts
//      //otherwise the entry must be divided by two, similarly to what done for the neighbours
//
//      SC center_val = -(1.0 - static_cast<SC>(1.0 / static_cast<SC>(region_belonging))) * region_val[iter_center_vector - region_col_vector.begin()];
//      ind_vector.push_back(center_ind);
//      val_vector.push_back(center_val);
//
//      std::cout << "Proc: " << comm_->getRank() << " --- " << __LINE__ << __FILE__ << std::endl;
//
//      regionMatrix->insertGlobalValues(region_node_idx-1, ind_vector, val_vector);
//    }
  }

  return;
}

template<class SC, class LO, class GO, class NO, Xpetra::UnderlyingLib lib, Xpetra::ESplittingMethodHHG splitMethod>
const bool Xpetra::RegionMatrix<SC,LO,GO,NO,lib,splitMethod>::hasRegionMatrix(
    const GO regID) const
{
  return not regionMatrices_[regID].is_null();
}

#endif /* XPETRA_REGION_MATRIX_IMPL_HPP_ */
