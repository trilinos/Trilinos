// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#if defined(STK_PERCEPT_HAS_GEOMETRY)
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#endif

#include <stk_mesh/base/FieldParallel.hpp>

#include <percept/math/DenseMatrix.hpp>

#include <percept/mesh/mod/smoother/SpacingFieldUtil.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/FieldTypes.hpp>

#include "mpi.h"

  namespace percept {

    void SpacingFieldUtil::compute_spacing_field()
    {
      stk::mesh::FieldBase *spacing_field    = m_eMesh.get_field(stk::topology::NODE_RANK, "ref_spacing_field");
      stk::mesh::FieldBase *spacing_field_counter    = m_eMesh.get_field(stk::topology::NODE_RANK, "ref_spacing_field_counter");

      m_eMesh.nodal_field_set_value("ref_spacing_field", 0.0);
      m_eMesh.nodal_field_set_value("ref_spacing_field_counter", 0.0);

      DenseMatrix<3,3> AI;

      int spatial_dim = m_eMesh.get_spatial_dim();

      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
      stk::mesh::Selector on_globally_shared_part =  ( m_eMesh.get_fem_meta_data()->globally_shared_part() );

      {
        // element loop: compute deltas
        const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k))
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

                    const MyPairIterRelation elem_nodes(m_eMesh, element, stk::topology::NODE_RANK);
                    unsigned num_node = elem_nodes.size();

                    for (unsigned inode=0; inode < num_node; inode++)
                      {
                        stk::mesh::Entity node = elem_nodes[ inode ].entity();
                        double *spacing = m_eMesh.field_data(spacing_field, node);
                        double *spacing_counter = m_eMesh.field_data(spacing_field_counter, node);
                        if (m_type == SPACING_AVE)
                          spacing_counter[0] += 1.0;
                        else
                          spacing_counter[0] = 1.0;

                        DenseMatrix<3,3>& A = jacA.m_J[inode];
                        inverse(A, AI);
                        for (int jdim=0; jdim < spatial_dim; jdim++)
                          {
                            double norm_AI_jdir=0.0;
                            for (int idim=0; idim < spatial_dim; idim++)
                              {
                                norm_AI_jdir += AI(idim,jdim)*AI(idim,jdim);
                              }
                            norm_AI_jdir = std::sqrt(norm_AI_jdir);
                            if (m_type == SPACING_AVE)
                              {
                                spacing[jdim] += 1.0/norm_AI_jdir;
                              }
                            else
                              spacing[jdim] = std::max(spacing[jdim], 1.0/norm_AI_jdir);
                          }
                      }
                  }
              }
          }

        std::vector<const stk::mesh::FieldBase*> spacing_field_vec(1,spacing_field);
        stk::mesh::parallel_sum(*m_eMesh.get_bulk_data(), spacing_field_vec);

        std::vector<const stk::mesh::FieldBase*> spacing_field_counter_vec(1,spacing_field_counter);
        stk::mesh::parallel_sum(*m_eMesh.get_bulk_data(), spacing_field_counter_vec);

        {
          std::vector< const stk::mesh::FieldBase *> fields;
          fields.push_back(spacing_field);
          fields.push_back(spacing_field_counter);

          // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
          stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->aura_ghosting(), fields);
        }

      }

      {
        // nodal loop
        const stk::mesh::BucketVector & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.node_rank() );

        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity node = bucket[i_node];
                    double *spacing = m_eMesh.field_data(spacing_field, node);
                    double *spacing_counter = m_eMesh.field_data(spacing_field_counter, node);

                    for (int idim=0; idim < spatial_dim; idim++)
                      {
                        spacing[idim] /= spacing_counter[0];
                      }
                  }
              }
          }
      }
    }

    double SpacingFieldUtil::spacing_at_node_in_direction(
            const double unit_vector_dir[3], stk::mesh::Entity node, stk::mesh::Selector *element_selector)
    {
      int spatial_dim = m_eMesh.get_spatial_dim();

      double spacing_ave=0.0;

      const MyPairIterRelation node_elems(m_eMesh, node, m_eMesh.element_rank() );
      double num_elem = 0;

      for (unsigned i_element = 0; i_element < node_elems.size(); i_element++)
        {
          stk::mesh::Entity element = node_elems[i_element].entity();

          if (!element_selector || (*element_selector)(m_eMesh.bucket(element)))
            {
              ++num_elem;
              const CellTopologyData *topology_data = m_eMesh.get_cell_topology(element);
              shards::CellTopology topo(topology_data);
              JacobianUtil jacA;

              double A_ = 0.0;
              jacA(A_, m_eMesh, element, m_eMesh.get_coordinates_field(), topology_data);

              unsigned inode = node_elems[i_element].relation_ordinal();
              const MyPairIterRelation elem_nodes(m_eMesh, element, m_eMesh.node_rank() );
              VERIFY_OP_ON(inode, <, elem_nodes.size(), "elem_nodes 2");
              VERIFY_OP_ON(m_eMesh.identifier(elem_nodes[inode].entity()), ==, m_eMesh.identifier(node), "elem_nodes 3");

              DenseMatrix<3,3>& A = jacA.m_J[inode];
              DenseMatrix<3,3> AI;
              inverse(A, AI);
              const double normA = Frobenius(A);

              VERIFY_OP_ON(normA, >=, 0.0, "ERROR: bad norm(Jacobian) in SpacingFieldUtil.");
              VERIFY_OP_ON(jacA.m_detJ[inode], >=, 1.e-12*normA, "ERROR: bad det(Jacobian) in SpacingFieldUtil.");

              double spacing[3] = {0,0,0};
              for (int idim=0; idim < spatial_dim; idim++)
                {
                  for (int jdim=0; jdim < spatial_dim; jdim++)
                    {
                      spacing[idim] += AI(idim,jdim)*unit_vector_dir[jdim];
                    }
                }
              double sum=0.0;
              for (int jdim=0; jdim < spatial_dim; jdim++)
                {
                  sum += spacing[jdim]*spacing[jdim];
                }
              sum = std::sqrt(sum);
              VERIFY_OP_ON(sum, >, 1.e-12*normA, "bad sum");
              spacing_ave += 1.0/sum;
            }
        }
      if (spacing_ave == 0.0) {
        std::cout << "spacing_ave= " << spacing_ave << std::endl;
        throw std::runtime_error("spacing too small before invert in SpacingFieldUtil::spacing_at_node_in_direction");
      }
      if (num_elem) spacing_ave /= num_elem;
      return spacing_ave;
    }
  }

#endif
