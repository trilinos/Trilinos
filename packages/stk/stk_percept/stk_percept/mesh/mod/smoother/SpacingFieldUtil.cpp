#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) 

#include <stk_percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

#include <stk_percept/math/DenseMatrix.hpp>

#include "SpacingFieldUtil.hpp"

#include "mpi.h"

namespace stk {
  namespace percept {

    void SpacingFieldUtil::compute_spacing_field()
    {
      stk::mesh::FieldBase *spacing_field    = m_eMesh.get_field("ref_spacing_field");
      stk::mesh::FieldBase *spacing_field_counter    = m_eMesh.get_field("ref_spacing_field_counter");

      m_eMesh.nodal_field_set_value(spacing_field, 0.0);
      m_eMesh.nodal_field_set_value(spacing_field_counter, 0.0);

      DenseMatrix<3,3> AI;

      int spatial_dim = m_eMesh.get_spatial_dim();

      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh.get_fem_meta_data()->globally_shared_part() );

      {
        // element loop: compute deltas
        const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( stk::mesh::MetaData::ELEMENT_RANK );

        for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (MeshSmoother::select_bucket(**k, &m_eMesh) && on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_elements_in_bucket = bucket.size();
                const CellTopologyData *topology_data = m_eMesh.get_cell_topology(bucket);

                for (unsigned i_element = 0; i_element < num_elements_in_bucket; i_element++)
                  {
                    stk::mesh::Entity element = bucket[i_element];
                    JacobianUtil jacA;

                    double A_ = 0.0;
                    jacA(A_, m_eMesh, element, m_eMesh.get_coordinates_field(), topology_data);

                    const mesh::PairIterRelation elem_nodes = element.relations( stk::mesh::MetaData::NODE_RANK );
                    unsigned num_node = elem_nodes.size();

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        mesh::Entity node = elem_nodes[ inode ].entity();
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
        stk::mesh::Selector sel(m_eMesh.get_fem_meta_data()->universal_part());
        stk::mesh::parallel_reduce(*m_eMesh.get_bulk_data(), stk::mesh::sum(*spacing_field_v, &sel));

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
                    stk::mesh::Entity node = bucket[i_node];
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

    double SpacingFieldUtil::spacing_at_node_in_direction(double unit_vector_dir[3], stk::mesh::Entity node, stk::mesh::Selector *element_selector)
    {
      int spatial_dim = m_eMesh.get_spatial_dim();

      double spacing_ave=0.0;
      double *coord = stk::mesh::field_data( *static_cast<const VectorFieldType *>(m_eMesh.get_coordinates_field()) , node );

      const mesh::PairIterRelation node_elems = node.relations( m_eMesh.element_rank() );
      double num_elem = 0;

      for (unsigned i_element = 0; i_element < node_elems.size(); i_element++)
        {
          stk::mesh::Entity element = node_elems[i_element].entity();
          if (!MeshSmoother::select_bucket(element.bucket(), &m_eMesh))
            continue;

          if (!element_selector || (*element_selector)(element))
            {
              ++num_elem;
              //if (element.entity_rank() == m_eMesh.element_rank() );
              const CellTopologyData *topology_data = m_eMesh.get_cell_topology(element);
              shards::CellTopology topo(topology_data);
              JacobianUtil jacA;

              DenseMatrix<3,3> AI;
              double A_ = 0.0;
              jacA(A_, m_eMesh, element, m_eMesh.get_coordinates_field(), topology_data);

              unsigned inode = node_elems[i_element].relation_ordinal();
              const mesh::PairIterRelation elem_nodes = element.relations( m_eMesh.node_rank() );
              VERIFY_OP_ON(inode, <, elem_nodes.size(), "elem_nodes 2");
              VERIFY_OP_ON(elem_nodes[inode].entity().identifier(), ==, node.identifier(), "elem_nodes 3");

              DenseMatrix<3,3>& A = jacA.m_J[inode];
              double detJ = jacA.m_detJ[inode];
              if (detJ < 1.e-12)
                {
                  std::cout << "ERROR: bad Jacobian in SpacingFieldUtil, detJ = " << detJ << std::endl;
                  for (unsigned jn=0; jn < elem_nodes.size(); jn++)
                    {
                      double *cr = stk::mesh::field_data( *static_cast<const VectorFieldType *>(m_eMesh.get_coordinates_field()) , elem_nodes[jn].entity() );
                      for (int idim=0; idim < spatial_dim; idim++)
                        {
                          std::cout << "coord[node-" << jn << ", " << idim << "]= " << cr[idim] << std::endl;
                        }
                    }
                  for (int idim=0; idim < spatial_dim; idim++)
                    {
                      for (int jdim=0; jdim < spatial_dim; jdim++)
                        {
                          std::cout << "(i,j)= " << idim << " " << jdim << " coord= " << coord[0] << " " << coord[1]
                                    << " AI= " << AI(idim,jdim) <<  " num_elem= " << num_elem << " i_element= " << i_element
                                    << " A_ = " << A_ << " A= " << A(idim,jdim) << " uv= " << unit_vector_dir[jdim] << " topo= " << topo.getName() << "\n";
                        }
                    }
                  VERIFY_OP_ON(detJ, >=, 1.e-12, "ERROR: bad Jacobian in SpacingFieldUtil.");
                }

              inverse(A, AI);
              double spacing[3] = {0,0,0};
              for (int idim=0; idim < spatial_dim; idim++)
                {
                  for (int jdim=0; jdim < spatial_dim; jdim++)
                    {
                      if (0)
                        std::cout << "(i,j)= " << idim << " " << jdim << " coord= " << coord[0] << " " << coord[1] 
                                  << " AI= " << AI(idim,jdim) <<  " num_elem= " << num_elem << " i_element= " << i_element 
                                  << " A_ = " << A_ << " A= " << A(idim,jdim) << " uv= " << unit_vector_dir[jdim] << " topo= " << topo.getName() << "\n";
                      spacing[idim] += AI(idim,jdim)*unit_vector_dir[jdim];
                    }
                }
              double sum=0.0;
              for (int jdim=0; jdim < spatial_dim; jdim++)
                {
                  if ( 0 && element.identifier() == 6659) 
                    {
                      std::cout << " nodeid= " << node.identifier() << " spacing[" << jdim << "]= " << spacing[jdim] << " unit_vector_dir= " << unit_vector_dir[jdim] << std::endl;
                      PerceptMesh::get_static_instance()->print(element);
                    }
                  sum += spacing[jdim]*spacing[jdim];
                }
              sum = std::sqrt(sum);
              VERIFY_OP_ON(sum, >, 1.e-12, "bad sum");
              spacing_ave += 1.0/sum;

              if (0 && element.identifier() == 6659) 
                {
                  PerceptMesh::get_static_instance()->print(element, false);
                  std::cout << " nodeid= " << node.identifier() << " spacing= " << 1.0/sum << std::endl;
                }

            }
        }
      //std::cout << "spacing_ave= " << spacing_ave << std::endl;
      if (spacing_ave < 1.e-10) {
        std::cout << "spacing_ave= " << spacing_ave << std::endl;
        throw std::runtime_error("spacing too small before invert in SpacingFieldUtil::spacing_at_node_in_direction");
      }
      if (num_elem) spacing_ave /= num_elem;
      return spacing_ave;
    }

  }
}

#endif
