// matrix_free_tpetra_operator.hpp

#ifndef DREAM_MATRIX_FREE_TPETRA_OPERATOR_HPP
#define DREAM_MATRIX_FREE_TPETRA_OPERATOR_HPP

#include <iostream>
#include <vector>
#include "Intrepid2_IntegrationTools.hpp"
#include "Intrepid2_CellGeometry.hpp"

// General TODOs:
// 1. Clean up template usages (replace doubles with scalar types, unsigned ints with LO/GO, etc.)
// 2. Check performance on GPU
// 3. Make sure it works on more than one element block

/**
 * This class defines an operator 
 */
class ProjectionOperator : public Tpetra::Operator<> {
public:
  // Tpetra::Operator subclasses should always define typedefs according to Tpetra, although I prefer CamelCase for types in general...
  typedef Tpetra::Operator<>::scalar_type ST;
  typedef Tpetra::Operator<>::local_ordinal_type LO;
  typedef Tpetra::Operator<>::global_ordinal_type GO;
  typedef Tpetra::Operator<>::node_type NO;
  typedef Tpetra::Vector<ST, LO, GO, NO> V;
  typedef Tpetra::MultiVector<ST, LO, GO, NO> MV;
  typedef Tpetra::Map<LO, GO, NO> map_type;
  typedef Tpetra::Import<LO, GO, NO> import_type;
  typedef Tpetra::Export<LO, GO, NO> export_type;
public:
  // Constructor
  ProjectionOperator(const Teuchos::RCP<const panzer_stk::STK_Interface> mesh,
                     const Teuchos::RCP<panzer::DOFManager> dof_manager,
                     const std::pair<Intrepid2::EOperator, Intrepid2::EOperator> physics_operator,
                     const Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space, double, double > > basis,
                     const unsigned int workset_size = 100)
  : mesh_(mesh),
    dof_manager_(dof_manager),
    physics_operator_(physics_operator),
    basis_(basis),
    workset_size_(workset_size)
  {
    auto basis_type = basis->getFunctionSpace();
    if (basis_type == Intrepid2::FUNCTION_SPACE_HGRAD)
      basis_type_ = "HGRAD";
    else if (basis_type == Intrepid2::FUNCTION_SPACE_HCURL)
      basis_type_ = "HCURL";
    else if (basis_type == Intrepid2::FUNCTION_SPACE_HDIV)
      basis_type_ = "HDIV";
    else if (basis_type == Intrepid2::FUNCTION_SPACE_HVOL)
      basis_type_ = "HVOL";
    fe_degree_ = basis->getDegree();

    auto comm = dof_manager->getComm();

    TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::invalid_argument, "ProjectionOperator constructor: The input Comm object must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION(mesh.is_null(), std::invalid_argument, "ProjectionOperator constructor: The input mesh object must be nonnull.");

    const int my_rank = comm->getRank();
    const int num_procs = comm->getSize();
    if (my_rank == 0) {
      std::cout << "ProjectionOperator constructor" << std::endl;
    }

    // -------------------------------------------------------------------
    // Call some setup routines if the user didn't specify a dof_manager
    // -------------------------------------------------------------------

    // if(dof_manager_ == Teuchos::null) {
    //   setup_dof_manager(comm);
    // }
    // this is the only non-const bit of the dof_manager_, but it's safer to do this than assume somebody else did it for us
    dof_manager_->useNeighbors(true); 
    Teuchos::RCP<const panzer::ConnManager> connection_manager = dof_manager_->getConnManager();

    // -------------------------------------------------------------------
    // Call some Intrepid2 setup routines
    // -------------------------------------------------------------------

    setup_intrepid2_values();

    // -----------------------------------------
    // Use the managers to get maps for DOF distribution
    // -----------------------------------------

    // Grab LIDs, get indices from dof_manager_, and compute the global number of DOFs
    LIDs = dof_manager_->getLIDs();
    std::vector<GO> owned, owned_and_shared;
    dof_manager_->getOwnedIndices(owned);
    dof_manager_->getOwnedAndGhostedIndices(owned_and_shared);
    LO num_unknowns = (LO) owned.size();
    GO local_num_unknowns = num_unknowns;
    Teuchos::reduceAll<LO,GO>(*comm, Teuchos::REDUCE_SUM, 1, &local_num_unknowns, &global_num_unknowns);
    std::cout << "Total number of " << global_num_unknowns << " unknowns across all processors" << std::endl;

    // Compute the necessary maps for the object
    const GO index_base = 0;
    operator_map_ = Teuchos::rcp(new map_type(global_num_unknowns, owned, index_base, comm));
    Teuchos::RCP<const map_type> owned_map = Teuchos::rcp(new map_type(global_num_unknowns, owned, 0, comm));
    Teuchos::RCP<const map_type> overlapped_map = Teuchos::rcp(new map_type(global_num_unknowns, owned_and_shared, 0, comm));
    importer_ = Teuchos::rcp(new import_type(owned_map, overlapped_map));

    // PRINT_VAR(operator_map_->getLocalNumElements())

    // column map for handling the redistribution
    Teuchos::ArrayView<const GO> element_list(owned_and_shared);
    redist_map_ = Teuchos::rcp(new map_type(global_num_unknowns, element_list, index_base, comm));

    // --------------------------------------------
    // Grab all elements and neighbors on my block
    // --------------------------------------------

    // TODO: this is not block-independent
    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);
    std::vector<stk::mesh::Entity> elems;
    mesh_->getMyElements(blocknames[0], elems);
    int space_dim = mesh_->getDimension();

    // grab neighbor elements and use a move iterator to move all elements into a single vector. this should be the most performant approach
    {
      std::vector<stk::mesh::Entity> neighbor_elems;
      mesh_->getNeighborElements(blocknames[0], neighbor_elems);
      elems.insert(elems.end(), std::make_move_iterator(neighbor_elems.begin()), std::make_move_iterator(neighbor_elems.end()));
    }
    std::vector<size_t> elem_IDs(elems.size());
    for(size_t i = 0; i < elems.size(); ++i)
      elem_IDs[i] = mesh_->elementLocalId(elems[i]);

    // ------------------------------------------------------
    // Setup the appropriate Intrepid2 information as needed
    // ------------------------------------------------------

    Teuchos::RCP<const shards::CellTopology> cell_topology = mesh_->getCellTopology(blocknames[0]);
    const int num_nodes_per_elem = cell_topology->getNodeCount();
    const unsigned int num_basis = ref_basis_vals.extent(0);
    const unsigned int num_ip = ref_ip.extent(0);

    const size_t num_elems = LIDs.extent(0);

    // Setup for discretization: jacobians, physical points, weights, and values
    Kokkos::DynRankView<double,PHX::Device> nodes("nodes", num_elems, num_nodes_per_elem, space_dim);

    mesh_->getElementVertices(elems, blocknames[0], nodes);

    Kokkos::DynRankView<int,PHX::Device> elem_to_node;
    std::cout << "Allocated?" << elem_to_node.is_allocated() << std::endl;

    if(mesh_.is_null()) {
      // TODO: Intrepid2 CellGeometry is all hardcoded to the 3D case, which is not optimal
      Kokkos::Array<double, 3> origin{0.,0.,0.};
      Kokkos::Array<double, 3> domain_extents{1.,1.,1.};
      Kokkos::Array<int, 3> grid_cell_counts{10,10,10};
      
      using Geometry = Intrepid2::CellGeometry<double, 3, PHX::Device>;
      Geometry::SubdivisionStrategy   subdivision_strategy = Geometry::NO_SUBDIVISION;
      Geometry::HypercubeNodeOrdering node_ordering        = Geometry::HYPERCUBE_NODE_ORDER_CLASSIC_SHARDS;
      geometry_ = Teuchos::rcp(new Geometry(origin, domain_extents, grid_cell_counts, subdivision_strategy, node_ordering));
    } else {
      bool claim_affine = false;
      geometry_ = Teuchos::rcp(new Intrepid2::CellGeometry<double, 3, PHX::Device>(*cell_topology, elem_to_node, nodes, claim_affine));
    }

    // Orientation nonsense
    auto no_connectivity_clone = connection_manager->noConnectivityClone();
    panzer::NodalFieldPattern pattern(*cell_topology);
    no_connectivity_clone->buildConnectivity(pattern);

    orientations = Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>("orientation", num_elems);
    auto orientations_host = create_mirror_view(orientations);
    for (size_t i = 0; i < elem_IDs.size(); ++i) {
      size_t elem_ID = elem_IDs[i];
      const GO* nodes = no_connectivity_clone->getConnectivity(elem_ID);
      Kokkos::View<GO*, Kokkos::DefaultHostExecutionSpace> node_view("nodes", num_nodes_per_elem);
      for (int node=0; node<num_nodes_per_elem; ++node)
        node_view(node) = nodes[node];
      orientations_host(elem_ID) = Intrepid2::Orientation::getOrientation(*cell_topology, node_view);
    }
    deep_copy(orientations,orientations_host);
  };

  // Required since we inherit from Tpetra::Operator
  // Destructor
  virtual ~ProjectionOperator() {}
  // Square matrix, so it's fine to use the same domain and range maps
  Teuchos::RCP<const map_type> getDomainMap() const { return operator_map_; }
  Teuchos::RCP<const map_type> getRangeMap() const { return operator_map_; }
  
  // Compute Y := alpha Op X + beta Y.
  void
  apply(const MV& X,
        MV& Y,
        Teuchos::ETransp mode = Teuchos::NO_TRANS,
        ST alpha = Teuchos::ScalarTraits<ST>::one(),
        ST beta = Teuchos::ScalarTraits<ST>::zero()) const
  {
    Teuchos::RCP<Teuchos::TimeMonitor> tm = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 2.1 - Operator Apply")));
    // ------------------------------
    // Step 1: Compute Y <- beta*Y
    // ------------------------------
    Y.scale(beta);
    
    // ------------------------------
    // Step 2: Compute Y <- Y + alpha*Op(X)
    // ------------------------------

    // Setup: get comms, ranks, procs
    Teuchos::RCP<const Teuchos::Comm<int> > comm = operator_map_->getComm();
    const int my_rank = comm->getRank();
    const int num_procs = comm->getSize();

    // Make a temporary multivector for holding the imported data
    const size_t num_vecs = X.getNumVectors();
    Teuchos::RCP<MV> redist_data_X = Teuchos::rcp(new MV(redist_map_, num_vecs));
    redist_data_X->doImport(X, *importer_, Tpetra::INSERT);

    // Get a view of the multivector
    auto kokkos_view_X = redist_data_X->getLocalView<NO::device_type>(Tpetra::Access::ReadOnly);
    auto kokkos_view_Y = Y.getLocalView<NO::device_type>(Tpetra::Access::ReadWrite);

    // setup bounds for indexing
    const unsigned int space_dim    = mesh_->getDimension();
    const unsigned int value_rank   = (basis_type_ == "HVOL" || basis_type_ == "HGRAD") ? 1 : space_dim;
    const unsigned int num_basis    = ref_basis_vals.extent(0);
    const unsigned int num_ip       = ref_ip.extent(0);
    const size_t       num_elems    = LIDs.extent(0);
    const unsigned int num_worksets = num_elems / workset_size_ + 1;

    // setup the minimum necessary containers for the operator
    Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values_temp;
    // TODO: only allocate what is absolutely necessary
    Intrepid2::Data<double,PHX::Device> ref_data = geometry_->getJacobianRefData(tensor_cubature_points);
    Intrepid2::Data<double,PHX::Device> jacobian = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    Intrepid2::Data<double,PHX::Device> jacobian_det        = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    Intrepid2::TensorData<double,PHX::Device> cell_measures = geometry_->allocateCellMeasure(jacobian_det, tensor_cubature_weights);

    // in the case of HGRAD, HDIV, or HCURL, we need the jacobian inverse
    Intrepid2::Data<double,PHX::Device> jacobian_inv;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD || physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
      jacobian_inv           = Intrepid2::CellTools<PHX::Device>::allocateJacobianInv(jacobian);
    // in the case of HDIV or HCURL, we need the reciprocal of the determinant
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
      jacobian_det_inv       = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    // in the case of HDIV or HCURL, we also need the jacobian divided by its determinant
    Intrepid2::Data<double,PHX::Device> jacobian_divided_by_det;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
      jacobian_divided_by_det = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv_extended; // container with same underlying data as jacobianDet, but extended with CONSTANT type to have same logical shape as Jacobian
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL) {
      auto variation_types = jacobian_det_inv.getVariationTypes(); // defaults to CONSTANT in ranks beyond the rank of the container; this is what we want for our new extents
      auto extents        = jacobian.getExtents();
      jacobian_det_inv_extended = jacobian_det_inv.shallowCopy(jacobian.rank(), extents, variation_types);
      jacobian_divided_by_det.storeInPlaceProduct(jacobian,jacobian_det_inv_extended);
    }
    
    // allocate integral data based on the case we're in
    Intrepid2::Data<double,PHX::Device> integral_data;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_VALUE) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHVOLtransformVALUE(jacobian_det, basis_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_grad_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHDIVtransformDIV(jacobian_det_inv, basis_div_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_inv, basis_grad_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    
    // Setup for the physics: get blocknames, elements, topology, etc.
    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);
    // TODO: finish some logic for grabbing the correct information from element blocks
    for(unsigned int i_block = 0; i_block < blocknames.size(); ++i_block)
    {
      // loop over worksets
      for(unsigned int elem_offset = 0; elem_offset < num_elems; elem_offset+=workset_size_)
      {
        // a bit of logic here; if we're in the last workset, we don't necessarily have the full workset_size_ of elements left
        const unsigned int num_workset_elems = (elem_offset + workset_size_ - 1 < num_elems) ? workset_size_ : num_elems - elem_offset;
        if(num_workset_elems != workset_size_)
        {
          const int CELL_DIM = 0;
          jacobian.setExtent(CELL_DIM, num_workset_elems);
          jacobian_det.setExtent(CELL_DIM, num_workset_elems);
          jacobian_inv.setExtent(CELL_DIM, num_workset_elems);
          jacobian_det_inv.setExtent(CELL_DIM, num_workset_elems);
          jacobian_divided_by_det.setExtent(CELL_DIM, num_workset_elems);
          integral_data.setExtent(CELL_DIM, num_workset_elems);
          cell_measures.setFirstComponentExtentInDimension0(num_workset_elems);
        }
        
        geometry_->setJacobian(jacobian, tensor_cubature_points, ref_data, elem_offset, elem_offset+num_workset_elems);
        Intrepid2::CellTools<PHX::Device>::setJacobianDet(   jacobian_det,        jacobian);
        if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD || physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
          Intrepid2::CellTools<PHX::Device>::setJacobianInv(   jacobian_inv,        jacobian);
        if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
          Intrepid2::CellTools<PHX::Device>::setJacobianDetInv(jacobian_det_inv,    jacobian);

        Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values;
        if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_VALUE)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformVALUE(num_workset_elems, basis_values);
        else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_grad_values);
        else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHDIVtransformDIV(jacobian_det_inv, basis_div_values);
        else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_divided_by_det, basis_curl_values);
        
        geometry_->computeCellMeasure(cell_measures, jacobian_det, tensor_cubature_weights);

        bool sum_into = false;

        // If we're in the case of VALUE*VALUE, or GRAD*GRAD, or DIV*DIV, or CURL*CURL
        if(physics_operator_.first == physics_operator_.second)
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data, transformed_basis_values, cell_measures, transformed_basis_values, sum_into);
        else
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Nonsymmetric physics operators are not yet supported!");

        auto computed_local_matrix = integral_data.getUnderlyingView3();
        LO min_index = operator_map_->getMinLocalIndex();
        LO max_index = operator_map_->getMaxLocalIndex();

        // Main loop: perform the matvec with the locally owned data
        // Some notes: since we're square, X and Y have the same map, so we can use LIDs to index each globally. Otherwise we would have additional work
        // TODO: Optimize loop orders
        //for(unsigned int e = 0; e < num_workset_elems; ++e) { // mesh element
	Kokkos::parallel_for("main assembly loop",Kokkos::RangePolicy<PHX::Device::execution_space>(0,num_workset_elems), KOKKOS_LAMBDA (const int e ) {
          for(unsigned int i = 0; i < num_basis; ++i) {
            int row_index = LIDs(elem_offset+e,i);
            if(row_index >= min_index && row_index <= max_index) { // if we're on a ghosted DOF, don't even bother trying to insert the value
              for(unsigned int j = 0; j < num_basis; ++j) {
                for(size_t c = 0; c < num_vecs; ++c) { // MV column
		  Kokkos::atomic_add(&kokkos_view_Y(row_index,c), alpha*computed_local_matrix(e,i,j)*kokkos_view_X(LIDs(elem_offset+e,j),c));
                }
              }
            }
          }
	});
      }
    }
  }

  bool
  hasDiagonal() const
  {
    return true;
  }

  // Compute Y := alpha Op X + beta Y.
  void
  getLocalDiagCopy(V& diag) const
  {
    Teuchos::RCP<Teuchos::TimeMonitor> tm = Teuchos::rcp(new Teuchos::TimeMonitor(*Teuchos::TimeMonitor::getNewTimer("Driver: 2.1 - Operator diagonal")));

    diag.putScalar(Teuchos::ScalarTraits<ST>::zero());

    // Setup: get comms, ranks, procs
    Teuchos::RCP<const Teuchos::Comm<int> > comm = operator_map_->getComm();
    const int my_rank = comm->getRank();
    const int num_procs = comm->getSize();

    // Get a view of the multivector
    auto kokkos_view_diag = diag.getLocalView<NO::device_type>(Tpetra::Access::ReadWrite);

    // setup bounds for indexing
    const unsigned int space_dim    = mesh_->getDimension();
    const unsigned int value_rank   = (basis_type_ == "HVOL" || basis_type_ == "HGRAD") ? 1 : space_dim;
    const unsigned int num_basis    = ref_basis_vals.extent(0);
    const unsigned int num_ip       = ref_ip.extent(0);
    const size_t       num_elems    = LIDs.extent(0);
    const unsigned int num_worksets = num_elems / workset_size_ + 1;

    // setup the minimum necessary containers for the operator
    Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values_temp;
    // TODO: only allocate what is absolutely necessary
    Intrepid2::Data<double,PHX::Device> ref_data = geometry_->getJacobianRefData(tensor_cubature_points);
    Intrepid2::Data<double,PHX::Device> jacobian = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    Intrepid2::Data<double,PHX::Device> jacobian_det        = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    Intrepid2::TensorData<double,PHX::Device> cell_measures = geometry_->allocateCellMeasure(jacobian_det, tensor_cubature_weights);

    // in the case of HGRAD, HDIV, or HCURL, we need the jacobian inverse
    Intrepid2::Data<double,PHX::Device> jacobian_inv;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD || physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
      jacobian_inv           = Intrepid2::CellTools<PHX::Device>::allocateJacobianInv(jacobian);
    // in the case of HDIV or HCURL, we need the reciprocal of the determinant
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
      jacobian_det_inv       = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    // in the case of HDIV or HCURL, we also need the jacobian divided by its determinant
    Intrepid2::Data<double,PHX::Device> jacobian_divided_by_det;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
      jacobian_divided_by_det = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);

    Intrepid2::Data<double,PHX::Device> jacobian_det_inv_extended; // container with same underlying data as jacobianDet, but extended with CONSTANT type to have same logical shape as Jacobian
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL) {
      auto variation_types = jacobian_det_inv.getVariationTypes(); // defaults to CONSTANT in ranks beyond the rank of the container; this is what we want for our new extents
      auto extents        = jacobian.getExtents();
      jacobian_det_inv_extended = jacobian_det_inv.shallowCopy(jacobian.rank(), extents, variation_types);
      jacobian_divided_by_det.storeInPlaceProduct(jacobian,jacobian_det_inv_extended);
    }

    // allocate integral data based on the case we're in
    Intrepid2::Data<double,PHX::Device> integral_data;
    if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_VALUE) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHVOLtransformVALUE(jacobian_det, basis_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_grad_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHDIVtransformDIV(jacobian_det_inv, basis_div_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }
    else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL) {
      auto transformed_values_temp = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_inv, basis_grad_values);
      integral_data = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp, cell_measures, transformed_values_temp);
    }

    // Setup for the physics: get blocknames, elements, topology, etc.
    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);
    // TODO: finish some logic for grabbing the correct information from element blocks
    for(unsigned int i_block = 0; i_block < blocknames.size(); ++i_block)
    {
      // loop over worksets
      for(unsigned int elem_offset = 0; elem_offset < num_elems; elem_offset+=workset_size_)
      {
        // a bit of logic here; if we're in the last workset, we don't necessarily have the full workset_size_ of elements left
        const unsigned int num_workset_elems = (elem_offset + workset_size_ - 1 < num_elems) ? workset_size_ : num_elems - elem_offset;
        if(num_workset_elems != workset_size_)
        {
          const int CELL_DIM = 0;
          jacobian.setExtent(CELL_DIM, num_workset_elems);
          jacobian_det.setExtent(CELL_DIM, num_workset_elems);
          jacobian_inv.setExtent(CELL_DIM, num_workset_elems);
          jacobian_det_inv.setExtent(CELL_DIM, num_workset_elems);
          jacobian_divided_by_det.setExtent(CELL_DIM, num_workset_elems);
          integral_data.setExtent(CELL_DIM, num_workset_elems);
          cell_measures.setFirstComponentExtentInDimension0(num_workset_elems);
        }

        geometry_->setJacobian(jacobian, tensor_cubature_points, ref_data, elem_offset, elem_offset+num_workset_elems);
        Intrepid2::CellTools<PHX::Device>::setJacobianDet(   jacobian_det,        jacobian);
        if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD || physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
          Intrepid2::CellTools<PHX::Device>::setJacobianInv(   jacobian_inv,        jacobian);
        if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
          Intrepid2::CellTools<PHX::Device>::setJacobianDetInv(jacobian_det_inv,    jacobian);

        Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values;
        if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_VALUE)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformVALUE(num_workset_elems, basis_values);
        else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_grad_values);
        else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHDIVtransformDIV(jacobian_det_inv, basis_div_values);
        else if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL)
          transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_divided_by_det, basis_curl_values);

        geometry_->computeCellMeasure(cell_measures, jacobian_det, tensor_cubature_weights);

        bool sum_into = false;

        // If we're in the case of VALUE*VALUE, or GRAD*GRAD, or DIV*DIV, or CURL*CURL
        if(physics_operator_.first == physics_operator_.second)
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data, transformed_basis_values, cell_measures, transformed_basis_values, sum_into);
        else
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Nonsymmetric physics operators are not yet supported!");

        auto computed_local_matrix = integral_data.getUnderlyingView3();
        LO min_index = operator_map_->getMinLocalIndex();
        LO max_index = operator_map_->getMaxLocalIndex();

        // Main loop: perform the matvec with the locally owned data
        // Some notes: since we're square, X and Y have the same map, so we can use LIDs to index each globally. Otherwise we would have additional work
        // TODO: Optimize loop orders
        //for(unsigned int e = 0; e < num_workset_elems; ++e) { // mesh element
	Kokkos::parallel_for("main assembly loop",Kokkos::RangePolicy<PHX::Device::execution_space>(0,num_workset_elems), KOKKOS_LAMBDA (const int e ) {
          for(unsigned int i = 0; i < num_basis; ++i) {
            int row_index = LIDs(elem_offset+e,i);
            if(row_index >= min_index && row_index <= max_index) { // if we're on a ghosted DOF, don't even bother trying to insert the value
              Kokkos::atomic_add(&kokkos_view_diag(row_index,0), computed_local_matrix(e,i,i));
            }
          }
	});
      }
    }
  }

  Tpetra::MultiVector<> getDofCoordinates() {
    return Tpetra::MultiVector<>();
  }

  GO getGlobalNumUnknowns() {
    return global_num_unknowns;
  }

  Teuchos::RCP<panzer::DOFManager> getDofManager() {
    return dof_manager_;
  }

  Teuchos::RCP<const panzer::ConnManager> getConnManager() {
    return dof_manager_->getConnManager();
  }

private:

  // Mesh
  const Teuchos::RCP<const panzer_stk::STK_Interface> mesh_;

  // Dof manager
  Teuchos::RCP<panzer::DOFManager> dof_manager_; // every dof manager has a conn manager via getConnManager
  
  // Physics
  const std::pair<Intrepid2::EOperator, Intrepid2::EOperator> physics_operator_;

  // Discretization
  std::string basis_type_;
  unsigned int fe_degree_;
  
  Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space, double, double > > basis_;
  Teuchos::RCP<Intrepid2::CellGeometry<double, 3, PHX::Device> > geometry_;
  Intrepid2::BasisValues<ST, PHX::Device> basis_values;
  Intrepid2::BasisValues<ST, PHX::Device> basis_grad_values;
  Intrepid2::BasisValues<ST, PHX::Device> basis_div_values;
  Intrepid2::BasisValues<ST, PHX::Device> basis_curl_values;
  Intrepid2::TensorPoints<double, PHX::Device> tensor_cubature_points;
  Intrepid2::TensorData<double,PHX::Device> tensor_cubature_weights;
  // keep the reference basis values and orientations, but forget everything else
  Kokkos::DynRankView<double,PHX::Device> ref_basis_vals;
  Kokkos::DynRankView<double,PHX::Device> ref_ip;
  Kokkos::DynRankView<double,PHX::Device> ref_weights;
  Kokkos::DynRankView<double,PHX::Device> phys_weights;
  Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> orientations;
  Kokkos::View<const LO**, Kokkos::LayoutRight, PHX::Device> LIDs;

  // Linear algebra
  Teuchos::RCP<const map_type> operator_map_, redist_map_; // every map has a comm via getComm
  Teuchos::RCP<const import_type> importer_;
  
  // Computational details
  GO global_num_unknowns; // requires an MPI communicate; compute this once
  unsigned int workset_size_;

  // Setup the DOF manager if one does not exist.
  // This requires creating a connection manager
  // void setup_dof_manager(const Teuchos::RCP<const Teuchos::MpiComm<int> > comm) {

  //   // -----------------------------------------
  //   // Create a connection manager and dof manager and place variables on the manager
  //   // -----------------------------------------

  //   vector<std::string> blocknames;
  //   mesh_->getElementBlockNames(blocknames);

  //   Teuchos::RCP<panzer::ConnManager> connection_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh_));
  //   dof_manager_ = Teuchos::rcp(new panzer::DOFManager());
  //   dof_manager_->useNeighbors(true); // use neighbors to avoid requiring a reduction across neighboring processes
  //   dof_manager_->setConnManager(connection_manager,*(comm->getRawMpiComm()));

  //   if(fe_degree_ > 2 || basis_type_ != "HGRAD")
  //     dof_manager_->setOrientationsRequired(true);

  //   // TODO GH: this is temporary; I will debug multiple blocks later
  //   if(blocknames.size() > 1) {
  //     std::cout << "Error! This mesh has " << blocknames.size() << " blocks!" << std::endl;
  //     TEUCHOS_TEST_FOR_EXCEPTION(blocknames.size()>1, std::invalid_argument, "The number of blocks in this mesh is greater than 1!");
  //   }

  //   // -----------------------------------------
  //   // Set the discretization and grab the quadrature and weights
  //   // -----------------------------------------

  //   int space_dim = mesh_->getDimension();
  //   std::vector<std::vector<size_t>> myElements;

  //   for (size_t block=0; block<blocknames.size(); ++block) {
  //     std::string block_name = blocknames[block];
  //     Teuchos::RCP<const shards::CellTopology> cell_topology = mesh_->getCellTopology(block_name);
  //     // std::string shape = cell_topology->getName();
  //     // int num_nodes_per_elem = cell_topology->getNodeCount();

  //     // std::vector<stk::mesh::Entity> elems;
  //     // mesh_->getMyElements(block_name, elems);

  //     // // list of all elements on this processor
  //     // std::vector<size_t> blockmyElements = vector<size_t>(elems.size());
  //     // for(size_t e=0; e<elems.size(); e++ ) {
  //     //   blockmyElements[e] = mesh_->elementLocalId(elems[e]);
  //     // }
  //     // myElements.push_back(blockmyElements);

  //     // Get the Intrepid2 basis and assign it to the problem
  //     basis_ = getIntrepid2Basis(space_dim, cell_topology->getName(), basis_type_, fe_degree_);
  //     Teuchos::RCP<const panzer::Intrepid2FieldPattern> pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_));
  //     dof_manager_->addField(blocknames[block], "p", pattern, panzer::FieldType::CG);
  //   }

  //   // this has to come after to tidy everything up
  //   dof_manager_->buildGlobalUnknowns();
  //   if (comm->getRank() == 0) {
  //     dof_manager_->printFieldInformation(std::cout);
  //   }

  // }

  // Compute basis values which are required by the physics
  void setup_intrepid2_values() {

    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);

    // -----------------------------------------
    // Set the discretization and grab the quadrature and weights
    // -----------------------------------------

    int space_dim = mesh_->getDimension();
    std::vector<std::vector<size_t>> myElements;

    for (size_t block=0; block<blocknames.size(); ++block) {
      std::string block_name = blocknames[block];
      Teuchos::RCP<const shards::CellTopology> cell_topology = mesh_->getCellTopology(block_name);

      // Intrepid2 uses the convention that quadrature order means the approximation quality
      // TODO: this depends on the operator, but since we generally use (value*value) or (grad*grad), 2*degree is the most appropriate quadrature for things.
      // TODO: actually, this should move up in priority since we could see runtimes be slower than they should be if this isn't optimized properly
      const Intrepid2::ordinal_type basis_size = basis_->getCardinality();
      const unsigned int quad_order = 2*fe_degree_;

      auto cubature = Intrepid2::DefaultCubatureFactory::create<PHX::Device>(*cell_topology,quad_order);
      tensor_cubature_weights = cubature->allocateCubatureWeights();
      tensor_cubature_points = cubature->allocateCubaturePoints();
      cubature->getCubature(tensor_cubature_points, tensor_cubature_weights);

      // -----------------------------
      // Loop over the entries in the physics operator and check which quantities need to be allocated
      // -----------------------------
      
      if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_VALUE || physics_operator_.second == Intrepid2::EOperator::OPERATOR_VALUE) {
        basis_values = basis_->allocateBasisValues(tensor_cubature_points, Intrepid2::OPERATOR_VALUE);
        basis_->getValues(basis_values, tensor_cubature_points, Intrepid2::OPERATOR_VALUE);
      }
      if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_GRAD || physics_operator_.second == Intrepid2::EOperator::OPERATOR_GRAD) {
        basis_grad_values = basis_->allocateBasisValues(tensor_cubature_points, Intrepid2::OPERATOR_GRAD);
        basis_->getValues(basis_grad_values, tensor_cubature_points, Intrepid2::OPERATOR_GRAD);
      }
      if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_DIV || physics_operator_.second == Intrepid2::EOperator::OPERATOR_DIV) {
        basis_div_values = basis_->allocateBasisValues(tensor_cubature_points, Intrepid2::OPERATOR_DIV);
        basis_->getValues(basis_div_values, tensor_cubature_points, Intrepid2::OPERATOR_DIV);
      }
      if(physics_operator_.first == Intrepid2::EOperator::OPERATOR_CURL || physics_operator_.second == Intrepid2::EOperator::OPERATOR_CURL) {
        basis_curl_values = basis_->allocateBasisValues(tensor_cubature_points, Intrepid2::OPERATOR_CURL);
        basis_->getValues(basis_curl_values, tensor_cubature_points, Intrepid2::OPERATOR_CURL);
      }

      // -----------------------------
      // Generate a cubature and weights for the problem
      // -----------------------------
      Intrepid2::DefaultCubatureFactory cubature_factory;
      Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space, double, double> > basis_cubature  = cubature_factory.create<PHX::Device::execution_space, double, double>(*cell_topology, quad_order);
      
      const int cubature_dim  = basis_cubature->getDimension();
      const int num_cubature_points = basis_cubature->getNumPoints();
      const unsigned int value_rank   = (basis_type_ == "HVOL" || basis_type_ == "HGRAD") ? 1 : space_dim;

      ref_ip = Kokkos::DynRankView<double,PHX::Device>("reference integration points", num_cubature_points, cubature_dim);
      ref_weights = Kokkos::DynRankView<double,PHX::Device>("reference weights", num_cubature_points);
      basis_cubature->getCubature(ref_ip, ref_weights);

      ref_basis_vals = Kokkos::DynRankView<double,PHX::Device>("basis values", basis_size, ref_ip.extent(0));
      basis_->getValues(ref_basis_vals, ref_ip, Intrepid2::OPERATOR_VALUE);
    }
  }

  void setup_orientations() {

    Teuchos::RCP<const panzer::ConnManager> connection_manager = dof_manager_->getConnManager();


  }
};

#endif // DREAM_MATRIX_FREE_TPETRA_OPERATOR_HPP
