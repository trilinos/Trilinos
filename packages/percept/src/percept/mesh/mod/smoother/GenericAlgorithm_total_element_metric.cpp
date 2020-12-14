#include <percept/mesh/mod/smoother/GenericAlgorithm_total_element_metric.hpp>

#include <percept/PerceptUtils.hpp>

#include <percept/mesh/mod/smoother/SmootherMetric.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>

namespace percept {

template<>
GenericAlgorithm_total_element_metric<STKMesh>::
GenericAlgorithm_total_element_metric(SmootherMetricImpl<STKMesh> *metric, PerceptMesh *eMesh, bool valid_in, size_t *num_invalid_in, Double mtot_in, size_t n_invalid_in)
  : m_metric(metric), m_eMesh(eMesh), valid(valid_in), num_invalid(num_invalid_in), mtot(mtot_in), n_invalid(n_invalid_in)
{
  coord_field = m_eMesh->get_coordinates_field();
  coord_field_current   = coord_field;
  coord_field_lagged  = m_eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

  cg_s_field    = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_s");

  on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
  on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
  spatialDim = m_eMesh->get_spatial_dim();

  {
    elements.resize(0);
    topos.resize(0);

    // element loop
    const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->element_rank() );

    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (MeshSmootherImpl<STKMesh>::select_bucket(**k, m_eMesh) && on_locally_owned_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;

            const unsigned num_elements_in_bucket = bucket.size();

            for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
              {
                stk::mesh::Entity element = bucket[i_element];
                if (m_eMesh->hasFamilyTree(element) && m_eMesh->isParentElement(element, true))
                  continue;
                elements.push_back(element);
                topos.push_back(m_eMesh->get_cell_topology(bucket));
              }
          }
      }
  }
  Kokkos::resize(element_invalid_flags, elements.size());
}

template<>
GenericAlgorithm_total_element_metric<StructuredGrid>::
GenericAlgorithm_total_element_metric(SmootherMetricImpl<StructuredGrid> *metric, PerceptMesh *eMesh, bool valid_in, size_t *num_invalid_in, Double mtot_in, size_t n_invalid_in)
  : m_metric(metric), m_eMesh(eMesh), valid(valid_in), num_invalid(num_invalid_in), mtot(mtot_in), n_invalid(n_invalid_in)
{
  std::shared_ptr<BlockStructuredGrid> bsg = m_eMesh->get_block_structured_grid();
  coord_field                              = bsg->m_fields["coordinates"].get();
  coord_field_current                      = coord_field;
  coord_field_lagged                       = bsg->m_fields["coordinates_lagged"].get();
  cg_s_field                               = bsg->m_fields["cg_s"].get();

  //on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
  //on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
  spatialDim = 3;

  bsg->get_elements(elements);
  topos.resize(elements.size(), static_cast<const typename StructuredGrid::MTCellTopology *>(0));

  Kokkos::resize(element_invalid_flags, elements.size());
}

template<typename MeshType>
void GenericAlgorithm_total_element_metric<MeshType>::
run(unsigned iBlock)
{
#if !defined(KOKKOS_ENABLE_CUDA)
  stk::diag::Timer root_timer("GATM", rootTimerStructured());
  stk::diag::TimeBlock root_block(root_timer);

  {
    stk::diag::Timer reduce1("GATM pr 1", root_timer);
    stk::diag::TimeBlock reduce1_block(reduce1);

    // main loop computes local and global metrics
    Double interim_tot = (mtot);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<SecondaryExecSpace>(0,elements.size()), *this, interim_tot);
    (mtot)=interim_tot;
  }

  {
  	stk::diag::Timer reduce1("GATM pr 2", root_timer);
  	stk::diag::TimeBlock reduce1_block(reduce1);

    // second loop computes n_invalid and valid data
  	size_t interim_n_invalid= (n_invalid);
    Kokkos::parallel_reduce(Kokkos::RangePolicy<SecondaryExecSpace>(0,elements.size()), KOKKOS_LAMBDA (unsigned index, size_t &local_n_invalid) {
        local_n_invalid += element_invalid_flags(index);
      }, interim_n_invalid);
    (n_invalid)=interim_n_invalid;
  }
#else
  for (unsigned index=0; index<elements.size(); index++) {
    operator()(index, mtot);
  }

  for (unsigned index=0; index<elements.size(); index++) {
    n_invalid += element_invalid_flags[index];
  }
#endif
  (valid) = (n_invalid==0);
}


// ETI
template
void GenericAlgorithm_total_element_metric<STKMesh>::
run(unsigned iBlock);

template
void GenericAlgorithm_total_element_metric<StructuredGrid>::
run(unsigned iBlock);

///////////////////////////////////////////////////////////////////////

template<typename MetricType>
SGridGenericAlgorithm_total_element_metric<MetricType>::
SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, Double mtot_in, size_t n_invalid_in)
  : m_metric(eMesh), mtot(mtot_in), n_invalid(n_invalid_in)
{
    if (eMesh->get_block_structured_grid()) {
        valid = false;
        unsigned noBlocks =
                eMesh->get_block_structured_grid()->m_sblocks.size();
        std::shared_ptr<BlockStructuredGrid> bsg =
                eMesh->get_block_structured_grid();

        block_sizes.resize(noBlocks);
        loop_orderings.resize(noBlocks);
        element_invalid_flags_per_block.resize(noBlocks);

        m_coord_field_current = bsg->m_fields["coordinates"].get();
        m_coord_field_original = bsg->m_fields["coordinates_NM1"].get();
        m_coords_current_iterator =
                *m_coord_field_current->m_block_fields[0]; //orient iterators at the beginning of blocks
        m_coords_original_iterator =
                *m_coord_field_original->m_block_fields[0];

        std::vector<StructuredCellIndex> elems_from_block;

        for (unsigned iBlock = 0; iBlock < noBlocks; iBlock++) { //populate loop orderings and block sizes for all blocks in mesh
            block_sizes[iBlock] = bsg->m_sblocks[iBlock]->m_sizes;
            loop_orderings[iBlock] = bsg->m_sblocks[iBlock]->m_loop_ordering;
            unsigned total_elems_this_block = (1
                    + block_sizes[iBlock].cell_max[0]
                    - block_sizes[iBlock].cell_min[0])
                    * (1 + block_sizes[iBlock].cell_max[1]
                            - block_sizes[iBlock].cell_min[1])
                    * (1 + block_sizes[iBlock].cell_max[2]
                            - block_sizes[iBlock].cell_min[2]);
            Kokkos::resize(
                    element_invalid_flags_per_block[iBlock],
                    total_elems_this_block);
        }
        m_iBlock= 0;;
    }
}

// ETI
template
SGridGenericAlgorithm_total_element_metric< SmootherMetricUntangleImpl<StructuredGrid> >::
SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, Double mtot_in, size_t n_invalid_in);

template
SGridGenericAlgorithm_total_element_metric<HexMeshSmootherMetric>::
SGridGenericAlgorithm_total_element_metric(PerceptMesh *eMesh, Double mtot_in, size_t n_invalid_in);

template<typename MetricType>
void SGridGenericAlgorithm_total_element_metric<MetricType>::run(
        unsigned iBlock) {
    m_iBlock= iBlock;
    unsigned total_elems_this_block = (1 + block_sizes[iBlock].cell_max[0]
            - block_sizes[iBlock].cell_min[0])
            * (1 + block_sizes[iBlock].cell_max[1]
                    - block_sizes[iBlock].cell_min[1])
            * (1 + block_sizes[iBlock].cell_max[2]
                    - block_sizes[iBlock].cell_min[2]);
    stk::diag::Timer root_timer("GATM"+std::to_string(total_elems_this_block), rootTimerStructured());
    stk::diag::TimeBlock root_block(root_timer);

    {
        flags_iterator = element_invalid_flags_per_block[iBlock];
        block_sizes_iterator = block_sizes[iBlock];
        loop_orderings_iterator[0] = loop_orderings[iBlock][0];
        loop_orderings_iterator[1] = loop_orderings[iBlock][1];
        loop_orderings_iterator[2] = loop_orderings[iBlock][2];

        m_coords_current_iterator =
                *m_coord_field_current->m_block_fields[iBlock]; //orient iterators at the beginning of blocks
        m_coords_original_iterator =
                *m_coord_field_original->m_block_fields[iBlock];


        stk::diag::Timer reduce1(
                std::string("GATM pr 1 ")
                        + std::to_string(
                                (int) std::cbrt((double) total_elems_this_block)),
                root_timer);
        stk::diag::TimeBlock reduce1_block(reduce1);
        // main loop computes local and global metrics
        mtot = 0;
        Double interim_mtot = mtot;
        Kokkos::parallel_reduce(
                Kokkos::RangePolicy<ExecSpace>(0, total_elems_this_block),
                *this, interim_mtot);

        mtot = interim_mtot;
    }

    {
        stk::diag::Timer reduce2(
                std::string("GATM pr 2 ")
                        + std::to_string(
                                (int) std::cbrt((double) total_elems_this_block)),
                                root_timer);
        stk::diag::TimeBlock reduce2_block(reduce2);
        // second loop computes n_invalid and valid data
        element_invalid_flags_reducer eifr;
        eifr.flags = flags_iterator;
        eifr.n_invalid = n_invalid;
        eifr.reduce();
        n_invalid = eifr.n_invalid;
    }
    valid = (n_invalid==0);
}

// ETI
template
void SGridGenericAlgorithm_total_element_metric<SmootherMetricUntangleImpl<StructuredGrid> >::run(unsigned iBlock);

template/*<>*/
void SGridGenericAlgorithm_total_element_metric<HexMeshSmootherMetric>::run(unsigned iBlock);


template<typename MetricType>
KOKKOS_INLINE_FUNCTION
void SGridGenericAlgorithm_total_element_metric<MetricType>::
operator()(const unsigned& index, Double& mtot_loc) const
{
  bool local_valid = false;
  StructuredGrid::MTElement elem;
  Kokkos::Array<unsigned,3> cell_ijk;
  sgrid_multi_dim_indices_from_index_cell(block_sizes_iterator, loop_orderings_iterator,index,cell_ijk);

  elem[0]=cell_ijk[0];
  elem[1]=cell_ijk[1];
  elem[2]=cell_ijk[2];
  Double mm = m_metric.metric(elem, local_valid);
  flags_iterator(index) = !local_valid;

  mtot_loc += mm;
}

// ETI
template
void SGridGenericAlgorithm_total_element_metric<SmootherMetricUntangleImpl<StructuredGrid> >::
operator()(const unsigned& index, Double& mtot_loc) const;

template<>
KOKKOS_INLINE_FUNCTION
void SGridGenericAlgorithm_total_element_metric<HexMeshSmootherMetric>::
operator()(const unsigned& index,Double& mtot_loc) const
{
  bool local_valid = false;
  double v_i_current[8][3];
  double v_i_org[8][3];
  unsigned indx[3] = { 0, 0, 0 };
  unsigned II[3] = { 0, 0, 0 };
  Kokkos::Array<unsigned,3> cell_ijk;
  const int A0 = 0, A1 = 1, A2 = 2;

  sgrid_multi_dim_indices_from_index_cell(block_sizes_iterator, loop_orderings_iterator,index,cell_ijk);

  unsigned cnt = 0;
  for (indx[2] = 0; indx[2] < 2; ++indx[2]) {
      II[2] = indx[2] + cell_ijk[2];
      for (indx[1] = 0; indx[1] < 2; ++indx[1]) {
          II[1] = indx[1] + cell_ijk[1];
          for (indx[0] = 0; indx[0] < 2; ++indx[0]) {
              II[0] = indx[0] + cell_ijk[0];
              for (unsigned ic = 0; ic < 3; ++ic) {
                  v_i_current[cnt][ic] = m_coords_current_iterator(II[A0],
                          II[A1], II[A2], ic);
                  v_i_org[cnt][ic] = m_coords_original_iterator(II[A0],
                          II[A1], II[A2], ic);
              }
              ++cnt;
          }
      }
  }

  Double mm = m_metric.metric(v_i_current, v_i_org, local_valid);
  flags_iterator(index) = !local_valid;
  mtot_loc += mm;
}

///////////////////////////////////////////////////////////////////////

template<>
KOKKOS_INLINE_FUNCTION
void GenericAlgorithm_total_element_metric<STKMesh>::
operator()(const unsigned& index, Double& mtot_loc)
{
  typename STKMesh::MTElement element = elements[index];
  // FIXME
  m_metric->m_topology_data = topos[index];

  //gather coords

  bool local_valid = false;
  Double mm = m_metric->metric(element, local_valid);
  element_invalid_flags(index) = !local_valid;

  mtot_loc += mm;
}

template<>
KOKKOS_INLINE_FUNCTION
void GenericAlgorithm_total_element_metric<StructuredGrid>::

operator()(const unsigned& index, Double& mtot_loc)
{

  typename StructuredGrid::MTElement element = elements[index];//{elements[index][0],elements[index][1],elements[index][2],elements[index][3]};
  //gather coords

  bool local_valid = false;
  Double mm = m_metric->metric(element, local_valid);
  element_invalid_flags(index) = !local_valid;

  mtot_loc += mm;

}

} // namespace percept

