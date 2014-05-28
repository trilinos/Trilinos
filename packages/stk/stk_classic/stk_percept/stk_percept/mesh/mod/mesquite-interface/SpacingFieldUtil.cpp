#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/mesh/mod/mesquite-interface/PerceptMesquiteMesh.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_percept/mesh/mod/mesquite-interface/PMMMsqMatrix.hpp>

#include "SpacingFieldUtil.hpp"

#include "mpi.h"

namespace MESQUITE_NS {

  extern int get_parallel_rank();
}

namespace stk {
  namespace percept {

    using namespace Mesquite;
    
    void SpacingFieldUtil::compute_spacing_field()
    {
      stk::mesh::FieldBase *spacing_field    = m_eMesh.get_field("ref_spacing_field");
      stk::mesh::FieldBase *spacing_field_counter    = m_eMesh.get_field("ref_spacing_field_counter");

      m_eMesh.nodal_field_set_value(spacing_field, 0.0);
      m_eMesh.nodal_field_set_value(spacing_field_counter, 0.0);

      MsqMatrix<3,3> AI;

      int spatial_dim = m_eMesh.get_spatial_dim();

      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh.get_fem_meta_data()->globally_shared_part() );

      {
        // element loop: compute deltas
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (PerceptMesquiteMesh::select_bucket(**k, &m_eMesh) && on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                const CellTopologyData *topology_data = m_eMesh.get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity& element = bucket[i_element];
                    JacobianUtil jacA;

                    double A_ = 0.0;
                    jacA(A_, m_eMesh, element, m_eMesh.get_coordinates_field(), topology_data);

                    const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::fem::FEMMetaData::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity & node = *elem_nodes[ inode ].entity();
                        double *spacing = PerceptMesh::field_data(spacing_field, node);
                        double *spacing_counter = PerceptMesh::field_data(spacing_field_counter, node);
                        if (m_type == SPACING_AVE)
                          spacing_counter[0] += 1.0;
                        else
                          spacing_counter[0] = 1.0;

                        inverse(jacA.m_J[inode], AI);
                        for (int jdim=0; jdim < spatial_dim; jdim++)
                          {
                            double sum=0.0;
                            for (int idim=0; idim < spatial_dim; idim++)
                              {
                                sum += AI(idim,jdim)*AI(idim,jdim);
                              }
                            sum = std::sqrt(sum);
                            if (m_type == SPACING_AVE)
                              spacing[jdim] += 1.0/sum;
                            else
                              spacing[jdim] = std::max(spacing[jdim], 1.0/sum);
                          }
                      }
                  }
              }
          }

        VectorFieldType *spacing_field_v = static_cast<VectorFieldType *>(spacing_field);
        stk::mesh::parallel_reduce(*m_eMesh.get_bulk_data(), stk::mesh::sum(*spacing_field_v));

        VectorFieldType *spacing_field_counter_v = static_cast<VectorFieldType *>(spacing_field_counter);
        stk::mesh::parallel_reduce(*m_eMesh.get_bulk_data(), stk::mesh::sum(*spacing_field_counter_v));

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(spacing_field);
          fields.push_back(spacing_field_counter);

          // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields); 

          // the shared part (just the shared boundary)
          //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
        }

      }

      {
        // nodal loop
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity& node = bucket[i_node];
                    double *spacing = PerceptMesh::field_data(spacing_field, node);
                    double *spacing_counter = PerceptMesh::field_data(spacing_field_counter, node);

                    for (int idim=0; idim < spatial_dim; idim++)
                      {
                        spacing[idim] /= spacing_counter[0];
                      }
                  }
              }
          }
      }


    }

  }
}

#endif
