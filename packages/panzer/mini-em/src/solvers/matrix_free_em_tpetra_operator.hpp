// matrix_free_em_tpetra_operator.hpp

#ifndef DREAM_MATRIX_FREE_EM_TPETRA_OPERATOR_HPP
#define DREAM_MATRIX_FREE_EM_TPETRA_OPERATOR_HPP

#include <iostream>
#include <vector>
#include "Intrepid2_IntegrationTools.hpp"
#include "Intrepid2_CellGeometry.hpp"
#include <MatrixMarket_Tpetra.hpp>

// print Kokkos::View dimensions
#define PRINT_VIEW_DIMS(view)                            \
std::cout << #view << " (";                              \
for(int r=0; r<view.rank()-1; ++r)                       \
  std::cout << view.extent(r) << ", ";                   \
std::cout << view.extent_int(view.rank()-1) << ")" << std::endl;

#define PRINT_VIEW_DIMS_INT(view)                        \
std::cout << #view << " (";                              \
for(int r=0; r<view.rank()-1; ++r)                       \
  std::cout << view.extent_int(r) << ", ";               \
std::cout << view.extent_int(view.rank()-1) << ")" << std::endl;

// print Kokkos::Views (defined as a macro to avoid template parameter hell) ((use responsibly))
#define PRINT_VIEW1(view)                                \
std::cout << #view << " (" << view.extent(0) << ") = ["; \
for(unsigned int i=0; i<view.extent(0)-1; ++i)           \
  std::cout << view(i) << ", ";                          \
std::cout << view(view.extent(0)-1) << "]" << std::endl;

#define PRINT_VIEW2(view)                                                         \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << ") = " << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i) {                                    \
  for(unsigned int j=0; j<view.extent(1); ++j)                                    \
    std::cout << view(i,j) << " ";                                                \
  std::cout << std::endl << " ";                                                  \
}

// some things like multivectors are "2D views" but only appear as 1D in practice, so the above print isn't as pretty
#define PRINT_VIEW2_LINEAR(view)                                                               \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << ") = [" << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i)                                      \
  for(unsigned int j=0; j<view.extent(1); ++j)                                    \
    std::cout << view(i,j) << " ";                                                \
std::cout << "]" << std::endl;

#define PRINT_VIEW3(view)                                                                                  \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << "," << view.extent(2) << ") = " << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i) {                                                             \
  for(unsigned int j=0; j<view.extent(1); ++j) {                                                           \
    for(unsigned int k=0; k<view.extent(2); ++k)                                                           \
      std::cout << view(i,j,k) << " ";                                                                     \
  std::cout << "]" << std::endl << "    [";                                                                \
  }                                                                                                        \
}                                                                                                          \
std::cout << "]" << std::endl;

#define PRINT_VIEW3_INT(view)                                                                              \
std::cout << #view << " (" << view.extent_int(0) << "," << view.extent_int(1) << "," << view.extent_int(2) << ") = " << std::endl;  \
for(int i=0; i<view.extent_int(0); ++i) {                                                         \
  for(int j=0; j<view.extent_int(1); ++j) {                                                       \
    for(int k=0; k<view.extent_int(2); ++k)                                                       \
      std::cout << view(i,j,k) << " ";                                                                     \
  std::cout << "]" << std::endl << "    [";                                                                \
  }                                                                                                        \
}                                                                                                          \
std::cout << "]" << std::endl;

#define PRINT_VIEW4_LINEAR(view)                                                                           \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << "," << view.extent(2) << "," << view.extent(3) << ") = " << std::endl;  \
for(unsigned int i=0; i<view.extent(0); ++i) {                                                             \
  for(unsigned int j=0; j<view.extent(1); ++j) {                                                           \
    for(unsigned int k=0; k<view.extent(2); ++k)                                                           \
      for(unsigned int l=0; l<view.extent(3); ++l)                                                         \
        std::cout << view(i,j,k,l) << " ";                                                                 \
  std::cout << "]" << std::endl << "    [";                                                                \
  }                                                                                                        \
}                                                                                                          \
std::cout << "]" << std::endl;

#define PRINT_VIEW4_LINEAR_INT(view)                                                                           \
std::cout << #view << " (" << view.extent_int(0) << "," << view.extent_int(1) << "," << view.extent_int(2) << "," << view.extent_int(3) << ") = " << std::endl;  \
for(unsigned int i=0; i<view.extent_int(0); ++i) {                                                             \
  for(unsigned int j=0; j<view.extent_int(1); ++j) {                                                           \
    for(unsigned int k=0; k<view.extent_int(2); ++k)                                                           \
      for(unsigned int l=0; l<view.extent_int(3); ++l)                                                         \
        std::cout << view(i,j,k,l) << " ";                                                                 \
  std::cout << "]" << std::endl << "    [";                                                                \
  }                                                                                                        \
}                                                                                                          \
std::cout << "]" << std::endl;

#define PRINT_DRV(view)                                                                                    \
std::cout << #view << " (" << view.extent(0) << "," << view.extent(1) << "," << view.extent(2) << ") = ["; \
for(unsigned int i=0; i<view.extent(0); ++i)                                                               \
  for(unsigned int j=0; j<view.extent(1); ++j)                                                             \
    for(unsigned int k=0; k<view.extent(2); ++k)                                                           \
      std::cout << view(i,j,k) << " ";                                                                     \
std::cout << "]" << std::endl;

#define PRINT_VAR(var)                         \
std::cout << #var << "=" << var << std::endl;

#define PRINT_VIEW2_MAX(view)                  \
double max = -100000;                          \
for(unsigned int i=0; i<view.extent(0); ++i)   \
  for(unsigned int j=0; j<view.extent(1); ++j) \
    if(view(i,j) > max)                        \
      max = view(i,j);                         \
std::cout << #view << " max=" << max << std::endl;


/**
 * This class defines an operator of the form c_1(\curl u, \curl v) + c_2(u, v)
 */
template<class ST,class LO,class GO,class NO>
class MatrixFreeEMOperator : public Tpetra::Operator<ST,LO,GO,NO> {
public:
  typedef Tpetra::Vector<ST, LO, GO, NO> V;
  typedef Tpetra::MultiVector<ST, LO, GO, NO> MV;
  typedef Tpetra::Map<LO, GO, NO> map_type;
  typedef Tpetra::Import<LO, GO, NO> import_type;
  typedef Tpetra::Export<LO, GO, NO> export_type;
public:
  // Constructor
  MatrixFreeEMOperator(const Teuchos::RCP<const panzer_stk::STK_Interface> mesh,
                       const Teuchos::RCP<panzer::DOFManager> dof_manager,
                       const std::pair<ST,ST> coefficients,
                       Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space, ST, ST > > basis,
                       const unsigned int workset_size = 100)
  : mesh_(mesh),
    dof_manager_(dof_manager),
    c1(coefficients.first),
    c2(coefficients.second),
    fe_degree_(basis->getDegree()),
    workset_size_(workset_size)
  {
    auto comm = dof_manager->getComm();
    basis_ = basis;
    TEUCHOS_TEST_FOR_EXCEPTION(comm.is_null(), std::invalid_argument, "MatrixFreeEMOperator constructor: The input Comm object must be nonnull.");
    
    // -------------------------------------------------------------------
    // Call some setup routines for the DOF manager
    // -------------------------------------------------------------------

    setup_dof_manager(comm);
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

    bool claim_affine = true;
    geometry_ = Teuchos::rcp(new Intrepid2::CellGeometry<double, 3, PHX::Device>(*cell_topology, elem_to_node, nodes, claim_affine));

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
  virtual ~MatrixFreeEMOperator() {}


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

    // Make a temporary multivector for holding the imported data
    const size_t num_vecs = X.getNumVectors();
    Teuchos::RCP<MV> redist_data_X = Teuchos::rcp(new MV(redist_map_, num_vecs));
    redist_data_X->doImport(X, *importer_, Tpetra::INSERT);

    // Get a view of the multivector
    auto kokkos_view_X = redist_data_X->getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto kokkos_view_Y = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);

    // setup bounds for indexing
    const unsigned int space_dim    = mesh_->getDimension();
    // const unsigned int value_rank   = space_dim;
    const unsigned int num_basis    = basis_values.extent(0);
    // const unsigned int num_ip       = ref_ip.extent(0);
    const size_t       num_elems    = LIDs.extent(0);
    // const unsigned int num_worksets = num_elems / workset_size_ + 1;

    // setup the minimum necessary containers for the operator
    Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values_temp;
    // TODO: only allocate what is absolutely necessary
    Intrepid2::Data<double,PHX::Device> ref_data = geometry_->getJacobianRefData(tensor_cubature_points);
    Intrepid2::Data<double,PHX::Device> jacobian = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    Intrepid2::Data<double,PHX::Device> jacobian_det        = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    Intrepid2::TensorData<double,PHX::Device> cell_measures = geometry_->allocateCellMeasure(jacobian_det, tensor_cubature_weights);

    // in the case of HGRAD, HDIV, or HCURL, we need the jacobian inverse
    Intrepid2::Data<double,PHX::Device> jacobian_inv = Intrepid2::CellTools<PHX::Device>::allocateJacobianInv(jacobian);
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    Intrepid2::Data<double,PHX::Device> jacobian_divided_by_det = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv_extended;
    {
      auto variation_types = jacobian_det_inv.getVariationTypes();
      auto extents        = jacobian.getExtents();
      jacobian_det_inv_extended = jacobian_det_inv.shallowCopy(jacobian.rank(), extents, variation_types);
      jacobian_divided_by_det.storeInPlaceProduct(jacobian,jacobian_det_inv_extended);
    }
    
    // allocate integral data based on the case we're in
    auto transformed_values_temp_curl = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_curl_values);
    //auto transformed_values_temp_curl = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_divided_by_det, basis_grad_values);
    Intrepid2::Data<double,PHX::Device> integral_data_curl = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp_curl, cell_measures, transformed_values_temp_curl);
    
    // Setup for the physics: get blocknames, elements, topology, etc.
    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);
    // std::cout << "Workset count start" << std::endl;
    // TODO: finish some logic for grabbing the correct information from element blocks
    for(unsigned int i_block = 0; i_block < blocknames.size(); ++i_block)
    {
      // loop over worksets
      for(unsigned int elem_offset = 0; elem_offset < num_elems; elem_offset+=workset_size_)
      {
        // a bit of logic here; if we're in the last workset, we don't necessarily have the full workset_size_ of elements left
        const unsigned int num_workset_elems = (elem_offset + workset_size_ - 1 < num_elems) ? workset_size_ : num_elems - elem_offset;
        if(num_workset_elems != workset_size_) 
          set_extents_for_workset(num_workset_elems, jacobian, jacobian_det, jacobian_inv, jacobian_det_inv, jacobian_divided_by_det, integral_data_curl, cell_measures);

        geometry_->setJacobian(jacobian, tensor_cubature_points, ref_data, elem_offset, elem_offset+num_workset_elems);
        jacobian.setExtent(2,space_dim);
        jacobian.setExtent(3,space_dim);

        Intrepid2::CellTools<PHX::Device>::setJacobianDet(   jacobian_det,        jacobian);
        Intrepid2::CellTools<PHX::Device>::setJacobianInv(   jacobian_inv,        jacobian);
        Intrepid2::CellTools<PHX::Device>::setJacobianDetInv(jacobian_det_inv,    jacobian);

        Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformVALUE(num_workset_elems, basis_values);
        Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_curls = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_curl_values);
        // Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformVALUE(jacobian_inv, basis_values);
        // Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_curls = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_divided_by_det, basis_curl_values);
        
        geometry_->computeCellMeasure(cell_measures, jacobian_det, tensor_cubature_weights);

        // PRINT_VAR(num_basis)
        // PRINT_VAR(num_workset_elems)
        // PRINT_VAR(elem_offset)
        // PRINT_VIEW2(basis_values)
        // PRINT_VIEW3(basis_curl_values)
        // PRINT_VIEW3_INT(transformed_basis_values)
        // PRINT_VIEW4_LINEAR_INT(transformed_basis_curls)
        // PRINT_VIEW2(cell_measures)
        // PRINT_VIEW_DIMS(basis_values)
        // PRINT_VIEW_DIMS(basis_curl_values)
        // PRINT_VIEW_DIMS_INT(transformed_basis_values)
        // PRINT_VIEW_DIMS_INT(transformed_basis_curls)

        bool sum_into = false;

        // We need two operators: CURLCURL and VALVAL
        // scale by constant c1
        auto cell_measures_data = cell_measures.getTensorComponent(0);
        if (abs(c1)>1e-15) {
          Intrepid2::Data<double,PHX::Device> scaled_c1_measures(c1,cell_measures_data.getExtents());
          cell_measures_data.storeInPlaceProduct(scaled_c1_measures, cell_measures_data);
          // integrate VALVAL
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_values, cell_measures, transformed_basis_values, sum_into);
          sum_into = true;
        }

        if (abs(c1)>1e-15) {
          // scale by constant c2/c1 (there may be a smarter way to do this, but it works currently)
          Intrepid2::Data<double,PHX::Device> scaled_c2_measures(c2/c1,cell_measures_data.getExtents());
          cell_measures_data.storeInPlaceProduct(scaled_c2_measures, cell_measures_data);
          // integrate GRADGRAD
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_curls,  cell_measures,  transformed_basis_curls, sum_into);
        } else {
          Intrepid2::Data<double,PHX::Device> scaled_c2_measures(c2,cell_measures_data.getExtents());
          cell_measures_data.storeInPlaceProduct(scaled_c2_measures, cell_measures_data);
          // integrate GRADGRAD
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_curls,  cell_measures,  transformed_basis_curls, sum_into);
        }

        // We need two operators: CURLCURL and VALVAL
        // scale by constant c1
        // auto cell_measures_data = cell_measures.getTensorComponent(0);
        // if (Teuchos::ScalarTraits<ST>::magnitude(c1)>1e-15) {
        //   Intrepid2::Data<double,PHX::Device> scaled_c1_measures(c1,cell_measures_data.getExtents());
        //   cell_measures_data.storeInPlaceProduct(scaled_c1_measures, cell_measures_data);
        //   // integrate VALVAL
        //   Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_values, cell_measures, transformed_basis_values, sum_into);
        //   sum_into = true;
        // }

        // if (Teuchos::ScalarTraits<ST>::magnitude(c1)>1e-15) {
        //   // scale by constant c2/c1 (there may be a smarter way to do this, but it works currently)
        //   Intrepid2::Data<double,PHX::Device> scaled_c2_measures;
        //   scaled_c2_measures = Intrepid2::Data<double,PHX::Device>(c2/c1,cell_measures_data.getExtents());
        //   cell_measures_data.storeInPlaceProduct(scaled_c2_measures, cell_measures_data);
        // } else {
        //   Intrepid2::Data<double,PHX::Device> scaled_c2_measures;
        //   scaled_c2_measures = Intrepid2::Data<double,PHX::Device>(c2,cell_measures_data.getExtents());
        //   cell_measures_data.storeInPlaceProduct(scaled_c2_measures, cell_measures_data);
        // }
        // // integrate GRADGRAD
        // Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_curls,  cell_measures,  transformed_basis_curls, sum_into);

        auto computed_local_matrix = integral_data_curl.getUnderlyingView3();

        // PRINT_VIEW3(computed_local_matrix)

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

    // Get a view of the multivector
    auto kokkos_view_diag = diag.getLocalViewDevice(Tpetra::Access::ReadWrite);

    // setup bounds for indexing
    const unsigned int space_dim    = mesh_->getDimension();
    // const unsigned int value_rank   = space_dim;
    const unsigned int num_basis    = basis_values.extent(0);
    // const unsigned int num_ip       = ref_ip.extent(0);
    const size_t       num_elems    = LIDs.extent(0);
    // const unsigned int num_worksets = num_elems / workset_size_ + 1;

    // setup the minimum necessary containers for the operator
    Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values_temp;
    // TODO: only allocate what is absolutely necessary
    Intrepid2::Data<double,PHX::Device> ref_data = geometry_->getJacobianRefData(tensor_cubature_points);
    Intrepid2::Data<double,PHX::Device> jacobian = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    Intrepid2::Data<double,PHX::Device> jacobian_det        = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    Intrepid2::TensorData<double,PHX::Device> cell_measures = geometry_->allocateCellMeasure(jacobian_det, tensor_cubature_weights);

    // in the case of HGRAD, HDIV, or HCURL, we need the jacobian inverse
    Intrepid2::Data<double,PHX::Device> jacobian_inv = Intrepid2::CellTools<PHX::Device>::allocateJacobianInv(jacobian);
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv = Intrepid2::CellTools<PHX::Device>::allocateJacobianDet(jacobian);
    Intrepid2::Data<double,PHX::Device> jacobian_divided_by_det = geometry_->allocateJacobianData(tensor_cubature_points, 0, workset_size_);
    Intrepid2::Data<double,PHX::Device> jacobian_det_inv_extended;
    {
      auto variation_types = jacobian_det_inv.getVariationTypes();
      auto extents        = jacobian.getExtents();
      jacobian_det_inv_extended = jacobian_det_inv.shallowCopy(jacobian.rank(), extents, variation_types);
      jacobian_divided_by_det.storeInPlaceProduct(jacobian,jacobian_det_inv_extended);
    }

    // allocate integral data based on the case we're in
    auto transformed_values_temp_curl = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_curl_values);
    //auto transformed_values_temp_curl = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_divided_by_det, basis_grad_values);
    Intrepid2::Data<double,PHX::Device> integral_data_curl = Intrepid2::IntegrationTools<PHX::Device>::allocateIntegralData(transformed_values_temp_curl, cell_measures, transformed_values_temp_curl);

    // Setup for the physics: get blocknames, elements, topology, etc.
    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);
    // std::cout << "Workset count start" << std::endl;
    // TODO: finish some logic for grabbing the correct information from element blocks
    for(unsigned int i_block = 0; i_block < blocknames.size(); ++i_block)
    {
      // loop over worksets
      for(unsigned int elem_offset = 0; elem_offset < num_elems; elem_offset+=workset_size_)
      {
        // a bit of logic here; if we're in the last workset, we don't necessarily have the full workset_size_ of elements left
        const unsigned int num_workset_elems = (elem_offset + workset_size_ - 1 < num_elems) ? workset_size_ : num_elems - elem_offset;
        if(num_workset_elems != workset_size_)
          set_extents_for_workset(num_workset_elems, jacobian, jacobian_det, jacobian_inv, jacobian_det_inv, jacobian_divided_by_det, integral_data_curl, cell_measures);

        geometry_->setJacobian(jacobian, tensor_cubature_points, ref_data, elem_offset, elem_offset+num_workset_elems);
        jacobian.setExtent(2,space_dim);
        jacobian.setExtent(3,space_dim);

        Intrepid2::CellTools<PHX::Device>::setJacobianDet(   jacobian_det,        jacobian);
        Intrepid2::CellTools<PHX::Device>::setJacobianInv(   jacobian_inv,        jacobian);
        Intrepid2::CellTools<PHX::Device>::setJacobianDetInv(jacobian_det_inv,    jacobian);

        Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformVALUE(num_workset_elems, basis_values);
        Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_curls = Intrepid2::FunctionSpaceTools<PHX::Device>::getHGRADtransformGRAD(jacobian_inv, basis_curl_values);
        // Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_values = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformVALUE(jacobian_inv, basis_values);
        // Intrepid2::TransformedBasisValues<double,PHX::Device> transformed_basis_curls = Intrepid2::FunctionSpaceTools<PHX::Device>::getHCURLtransformCURL(jacobian_divided_by_det, basis_curl_values);

        geometry_->computeCellMeasure(cell_measures, jacobian_det, tensor_cubature_weights);

        // PRINT_VAR(num_basis)
        // PRINT_VAR(num_workset_elems)
        // PRINT_VAR(elem_offset)
        // PRINT_VIEW2(basis_values)
        // PRINT_VIEW3(basis_curl_values)
        // PRINT_VIEW3_INT(transformed_basis_values)
        // PRINT_VIEW4_LINEAR_INT(transformed_basis_curls)
        // PRINT_VIEW2(cell_measures)
        // PRINT_VIEW_DIMS(basis_values)
        // PRINT_VIEW_DIMS(basis_curl_values)
        // PRINT_VIEW_DIMS_INT(transformed_basis_values)
        // PRINT_VIEW_DIMS_INT(transformed_basis_curls)

        bool sum_into = false;

        // We need two operators: CURLCURL and VALVAL
        // scale by constant c1
        auto cell_measures_data = cell_measures.getTensorComponent(0);
        if (abs(c1)>1e-15) {
          Intrepid2::Data<double,PHX::Device> scaled_c1_measures(c1,cell_measures_data.getExtents());
          cell_measures_data.storeInPlaceProduct(scaled_c1_measures, cell_measures_data);
          // integrate VALVAL
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_values, cell_measures, transformed_basis_values, sum_into);
          sum_into = true;
        }

        if (abs(c1)>1e-15) {
          // scale by constant c2/c1 (there may be a smarter way to do this, but it works currently)
          Intrepid2::Data<double,PHX::Device> scaled_c2_measures(c2/c1,cell_measures_data.getExtents());
          cell_measures_data.storeInPlaceProduct(scaled_c2_measures, cell_measures_data);
          // integrate GRADGRAD
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_curls,  cell_measures,  transformed_basis_curls, sum_into);
        } else {
          Intrepid2::Data<double,PHX::Device> scaled_c2_measures(c2,cell_measures_data.getExtents());
          cell_measures_data.storeInPlaceProduct(scaled_c2_measures, cell_measures_data);
          // integrate GRADGRAD
          Intrepid2::IntegrationTools<PHX::Device>::integrate(integral_data_curl, transformed_basis_curls,  cell_measures,  transformed_basis_curls, sum_into);
        }

        auto computed_local_matrix = integral_data_curl.getUnderlyingView3();

        // PRINT_VIEW3(computed_local_matrix)

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
              // std::cout << "h " << computed_local_matrix(e,i,i) << std::endl;
            }
          }
        });
      }
    }
    // for (size_t i = 0; i<kokkos_view_diag.extent(0); i++)
    //   std::cout <<kokkos_view_diag(i,0)<< std::endl;
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<ST,LO,GO,NO> >::writeDenseFile("diag", diag);
  }

private:

  // Mesh
  const Teuchos::RCP<const panzer_stk::STK_Interface> mesh_;

  // Dof manager
  Teuchos::RCP<panzer::DOFManager> dof_manager_; // every dof manager has a conn manager via getConnManager

  // Coefficients for physics
  const ST c1;
  const ST c2;

  // Discretization
  const unsigned int fe_degree_;
  
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
  // This requires creating a connection manager, or using an existing one
  void setup_dof_manager(const Teuchos::RCP<const Teuchos::Comm<int> > comm) {

    // -----------------------------------------
    // Create a connection manager and dof manager and place variables on the manager
    // -----------------------------------------
    std::vector<std::string> blocknames;
    mesh_->getElementBlockNames(blocknames);

    // if we have a DOFManager already, grab the connection manager from it and then throw away the user-provided DOFManager
    Teuchos::RCP<panzer::ConnManager> connection_manager;
    if(dof_manager_ == Teuchos::null)
      connection_manager = Teuchos::rcp(new panzer_stk::STKConnManager(mesh_));
    else
      connection_manager = Teuchos::rcp_const_cast<panzer::ConnManager>(dof_manager_->getConnManager());

    dof_manager_ = Teuchos::rcp(new panzer::DOFManager());
    dof_manager_->useNeighbors(true); // use neighbors to avoid requiring a reduction across neighboring processes
    dof_manager_->setConnManager(connection_manager,*(Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm,true)->getRawMpiComm()));
    dof_manager_->setOrientationsRequired(true);

    // TODO GH: this is temporary; I will debug multiple blocks later
    if(blocknames.size() > 1) {
      std::cout << "Error! This mesh has " << blocknames.size() << " blocks!" << std::endl;
      TEUCHOS_TEST_FOR_EXCEPTION(blocknames.size()>1, std::invalid_argument, "The number of blocks in this mesh is greater than 1!");
    }
    
    // -----------------------------------------
    // Set the discretization and grab the quadrature and weights
    // -----------------------------------------

    // int space_dim = mesh_->getDimension();
    std::vector<std::vector<size_t>> myElements;

    for (size_t block=0; block<blocknames.size(); ++block) {
      std::string block_name = blocknames[block];
      Teuchos::RCP<const shards::CellTopology> cell_topology = mesh_->getCellTopology(block_name);
      // std::string shape = cell_topology->getName();
      // int num_nodes_per_elem = cell_topology->getNodeCount();
      
      // std::vector<stk::mesh::Entity> elems;
      // mesh_->getMyElements(block_name, elems);
      
      // // list of all elements on this processor
      // std::vector<size_t> blockmyElements = vector<size_t>(elems.size());
      // for(size_t e=0; e<elems.size(); e++ ) {
      //   blockmyElements[e] = mesh_->elementLocalId(elems[e]);
      // }
      // myElements.push_back(blockmyElements);
      
      // Get the Intrepid2 basis and assign it to the problem
      // basis_ = getIntrepid2Basis(space_dim, cell_topology->getName(), "HGRAD", fe_degree_);
      Teuchos::RCP<const panzer::Intrepid2FieldPattern> pattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(basis_));
      dof_manager_->addField(blocknames[block], "p", pattern, panzer::FieldType::CG);
    }
      
    // this has to come after to tidy everything up
    dof_manager_->buildGlobalUnknowns();
    if (comm->getRank() == 0) {
      dof_manager_->printFieldInformation(std::cout);
    }

  }

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
      // const Intrepid2::ordinal_type basis_size = basis_->getCardinality();
      const unsigned int quad_order = 2*fe_degree_;

      Teuchos::RCP<Intrepid2::Cubature<PHX::Device,double,double> > cubature = Intrepid2::DefaultCubatureFactory::create<PHX::Device>(*cell_topology,quad_order);
      tensor_cubature_weights = cubature->allocateCubatureWeights();
      tensor_cubature_points = cubature->allocateCubaturePoints();
      cubature->getCubature(tensor_cubature_points, tensor_cubature_weights);
      
      // -----------------------------
      // Loop over the entries in the physics operator and check which quantities need to be allocated
      // -----------------------------
      
      basis_values = basis_->allocateBasisValues(tensor_cubature_points, Intrepid2::OPERATOR_VALUE);
      basis_->getValues(basis_values, tensor_cubature_points, Intrepid2::OPERATOR_VALUE);
      basis_curl_values = basis_->allocateBasisValues(tensor_cubature_points, Intrepid2::OPERATOR_GRAD);
      basis_->getValues(basis_curl_values, tensor_cubature_points, Intrepid2::OPERATOR_GRAD);

      // -----------------------------
      // Generate a cubature and weights for the problem
      // -----------------------------
      Intrepid2::DefaultCubatureFactory cubature_factory;
      Teuchos::RCP<Intrepid2::Cubature<PHX::Device::execution_space, double, double> > basis_cubature  = cubature_factory.create<PHX::Device::execution_space, double, double>(*cell_topology, quad_order);
      
      const int cubature_dim  = basis_cubature->getDimension();
      // PRINT_VAR(cubature_dim)
      const int num_cubature_points = basis_cubature->getNumPoints();
      // PRINT_VAR(num_cubature_points)
      // const unsigned int value_rank   = space_dim;

      ref_ip = Kokkos::DynRankView<double,PHX::Device>("reference integration points", num_cubature_points, cubature_dim);
      ref_weights = Kokkos::DynRankView<double,PHX::Device>("reference weights", num_cubature_points);
      basis_cubature->getCubature(ref_ip, ref_weights);
      
      // ref_basis_vals = Kokkos::DynRankView<double,PHX::Device>("basis values", basis_size, ref_ip.extent(0));
      // basis_->getValues(ref_basis_vals, ref_ip, Intrepid2::OPERATOR_VALUE);
      // ref_basis_vals = Kokkos::DynRankView<double,PHX::Device>("basis grads", basis_size, ref_ip.extent(0), ref_ip.);
      // basis_->getValues(ref_basis_vals, ref_ip, Intrepid2::OPERATOR_VALUE);
    }

    // PRINT_VIEW2(tensor_cubature_points)
  }

  KOKKOS_INLINE_FUNCTION
  void set_extents_for_workset(const unsigned int num_workset_elems,
                               Intrepid2::Data<double,PHX::Device> &jacobian, 
                               Intrepid2::Data<double,PHX::Device> &jacobian_det, 
                               Intrepid2::Data<double,PHX::Device> &jacobian_inv, 
                               Intrepid2::Data<double,PHX::Device> &jacobian_det_inv, 
                               Intrepid2::Data<double,PHX::Device> &jacobian_divided_by_det, 
                               Intrepid2::Data<double,PHX::Device> &integral_data, 
                               Intrepid2::TensorData<double,PHX::Device> &cell_measures) const
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

};

#endif // DREAM_MATRIX_FREE_EM_TPETRA_OPERATOR_HPP
