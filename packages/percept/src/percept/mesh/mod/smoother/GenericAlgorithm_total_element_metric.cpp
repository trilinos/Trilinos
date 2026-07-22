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
  element_invalid_flags.resize(elements.size());
}

// ETI
template
void GenericAlgorithm_total_element_metric<STKMesh>::
run(unsigned iBlock);

///////////////////////////////////////////////////////////////////////
// explicit specialization of 'operator()' after instantiation
template<>
inline
void GenericAlgorithm_total_element_metric<STKMesh>::
operator()(const unsigned& index, Double& mtot_loc)
{
  typename STKMesh::MTElement element = elements[index];
  // FIXME
  m_metric->m_topology_data = topos[index];

  //gather coords

  bool local_valid = false;
  Double mm = m_metric->metric(element, local_valid);
  element_invalid_flags[index] = !local_valid;

  mtot_loc += mm;
}

template<typename MeshType>
void GenericAlgorithm_total_element_metric<MeshType>::
run(unsigned /*iBlock*/)
{
  for (unsigned index=0; index<elements.size(); index++) {
    operator()(index, mtot); // implicit instantiation first required here
  }

  for (unsigned index=0; index<elements.size(); index++) {
    n_invalid += element_invalid_flags[index];
  }

  (valid) = (n_invalid==0);
}

} // namespace percept

