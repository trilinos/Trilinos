#include <percept/mesh/mod/smoother/GenericAlgorithm_update_coordinates.hpp>

#include <percept/PerceptUtils.hpp>
#include <percept/PerceptMesh.hpp>

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherConjugateGradient.hpp>

namespace percept {

template<>
GenericAlgorithm_update_coordinates<STKMesh>::
GenericAlgorithm_update_coordinates(RefMeshSmoother *rms, PerceptMesh *eMesh, Double alpha_in)
  : m_rms(rms), m_eMesh(eMesh), alpha(alpha_in)
{
	 stk::diag::Timer     cumulative_timer(eMesh->getProperty("in_filename"), rootTimerStructured());
	 stk::diag::Timer constrGATM("GATM <STKMesh> constructor",cumulative_timer);
	 stk::diag::TimeBlock my_time_block(constrGATM);

  dmax=m_rms->m_dmax;
  dmax_relative=m_rms->m_dmax_relative;
  m_calc_deltas=false;

  coord_field = m_eMesh->get_coordinates_field();
  coord_field_current   = coord_field;
  coord_field_lagged  = m_eMesh->get_field(stk::topology::NODE_RANK, "coordinates_lagged");

  cg_s_field = m_eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
  cg_edge_length_field = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

  on_locally_owned_part =  ( m_eMesh->get_fem_meta_data()->locally_owned_part() );
  on_globally_shared_part =  ( m_eMesh->get_fem_meta_data()->globally_shared_part() );
  spatialDim = m_eMesh->get_spatial_dim();

  {
    std::vector< const stk::mesh::FieldBase *> fields;
    fields.push_back(cg_s_field);
    fields.push_back(coord_field);
    // stk::mesh::copy_owned_to_shared(*m_eMesh->get_bulk_data(), fields);
    // stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
    MTcommFields<STKMesh>(fields, m_eMesh);
  }

  // cache coordinates
  m_eMesh->copy_field("coordinates_lagged", "coordinates");

  nodes.resize(0);

  // node loop - update positions
  {
    const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
    for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        // update local and globally shared
        if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                nodes.push_back(node);
              }
          }
      }
  }

  //filter out fixed nodes
  std::vector<typename STKMesh::MTNode> unFixedNodes;
  unFixedNodes.resize(nodes.size());
  std::pair<bool,int> fixed;
  int64_t numUnFixed = 0;

  for (int64_t iNode = 0; iNode < (int64_t)nodes.size(); ++iNode)
	{
		fixed = m_rms->get_fixed_flag(nodes[iNode]);
		if(!fixed.first)
		{
			unFixedNodes[numUnFixed]=nodes[iNode];
			numUnFixed++;
		}
	}
 nodes.resize(numUnFixed);

 for (int64_t iNode = 0; iNode < numUnFixed; ++iNode)
	{ nodes[iNode]=unFixedNodes[iNode]; } //only operate on unFixedNodes

}//GATM constrcr <STKMesh>


template<>
KOKKOS_INLINE_FUNCTION
void GenericAlgorithm_update_coordinates<STKMesh>::
operator()(const int64_t& index) const
{
	  STKMesh::MTNode node = nodes[index];
	  double coord_current[3];
	  double cg_s[3];
	  get_field<STKMesh>(coord_current, spatialDim, m_eMesh, coord_field_current, node);
	  get_field<STKMesh>(cg_s, spatialDim, m_eMesh, cg_s_field, node);
	  double cg_edge_length[1];

	  //double coord_project[3] = {0,0,0};
	  for (int i=0; i < spatialDim; i++)
	    {
	      Double dt = alpha * cg_s[i];
	      coord_current[i] += dt;
	      if(m_calc_deltas)
	      {
	          get_field<STKMesh>(cg_edge_length, 1, m_eMesh, cg_edge_length_field, node);
	          m_rms->m_dmax = std::max(std::fabs(dt), m_rms->m_dmax);
	          m_rms->m_dmax_relative = std::max(std::fabs(dt)/cg_edge_length[0], m_rms->m_dmax_relative);
	      }
	    }
	  set_field<STKMesh>(coord_current, spatialDim, m_eMesh, coord_field_current, node);
}


template<>
void GenericAlgorithm_update_coordinates<STKMesh>::
run(bool calc_deltas)
{
      m_calc_deltas=calc_deltas;
      if(m_calc_deltas)
      {//ensure deltas are updated before calculating them
          m_rms->m_dmax=0;
          m_rms->m_dmax_relative=0;
          dmax=m_rms->m_dmax;
          dmax_relative=m_rms->m_dmax_relative;
      }

	  {
	   	  std::vector< const typename STKMesh::MTField *> fields;
	   	  fields.push_back(cg_s_field);
	   	  fields.push_back(coord_field);
	      MTcommFields<STKMesh>(fields, m_eMesh);
	  }

	  // cache coordinates
	  m_eMesh->copy_field("coordinates_lagged", "coordinates");

	  for (int64_t index = 0; index < (int64_t)nodes.size(); ++index)
	   {
	     (*this)(index);
	   }


}

template struct GenericAlgorithm_update_coordinates<STKMesh>;

} // namespace percept

