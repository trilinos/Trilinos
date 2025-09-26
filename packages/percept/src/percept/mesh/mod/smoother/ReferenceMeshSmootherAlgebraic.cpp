// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#if !defined(NO_GEOM_SUPPORT)

#include <percept/mesh/mod/smoother/ReferenceMeshSmootherAlgebraic.hpp>
#include <percept/mesh/mod/smoother/MeshSmoother.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>
#include <percept/math/DenseMatrix.hpp>

#include <stk_mesh/base/FieldParallel.hpp>
#include <stdio.h>
#include <queue>
#include <limits>

#include "mpi.h"
#include <cstdio>

#define DEBUG_PRINT 0
#define PRINT(a) do { if (DEBUG_PRINT && !m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_1(a) do { if (!m_eMesh->get_rank()) std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)
#define PRINT_2(a) do {  std::cout << "P[" << m_eMesh->get_rank() <<"] " << a << std::endl; } while(0)

namespace percept {

  int ReferenceMeshSmootherAlgebraic::find_new_value(stk::mesh::Entity node, int valOld, WallDistanceFieldType *wall_distance_field, stk::mesh::FieldBase */*coord_field_orig*/)
  {
    int valNew = valOld;
    typedef std::set<stk::mesh::Entity> EntitySet;
    EntitySet neighbors;
    //double dmin = std::numeric_limits<double>::max();
    m_eMesh->get_node_node_neighbors(node, neighbors);
    for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
      {
        stk::mesh::Entity nnode = *it;
        WallDistanceFieldType::value_type *valn = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, nnode);
        if (valn[0] != 0)
          {
            int valTmp = valn[0] + 1;
            if (valTmp == 0) valTmp = 1;
            if (valNew != 0)
              valTmp = std::min(valNew, valTmp);
            if (valNew != valTmp)
              {
                valNew = valTmp;
              }
          }
      }
    return valNew;
  }

  /// find (integer) distance to wall
  void ReferenceMeshSmootherAlgebraic::get_wall_distances()
  {
    PerceptMesh *eMesh = m_eMesh;
    stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    WallDistanceFieldType *wall_distance_field = eMesh->m_wall_distance_field;
    stk::mesh::FieldBase *coord_field_orig   = m_coord_field_original;
    stk::mesh::FieldBase *coord_field   = m_eMesh->get_coordinates_field();

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );

    typedef std::queue<stk::mesh::Entity> EntityQ;
    EntityQ nodeQ;

    {
      // node loop
      std::vector<stk::mesh::Entity> nodes;
      const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
      for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          //if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          //if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket = bucket.size();

            for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
              {
                stk::mesh::Entity node = bucket[i_node];
                VERIFY_OP_ON(m_eMesh->is_valid(node), ==, true, "bad node");
                WallDistanceFieldType::value_type *val = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, node);
                val[0] = 0;

                double *current_coord = static_cast<double*>(stk::mesh::field_data(*coord_field, node));
                double *orig_coord = static_cast<double*>(stk::mesh::field_data(*coord_field_orig, node));
                double d = 0.0;
                for (int jc = 0; jc < m_eMesh->get_spatial_dim(); ++jc)
                  {
                    d += (current_coord[jc] - orig_coord[jc])*(current_coord[jc] - orig_coord[jc]);
                  }
                std::pair<bool,int> fixed = this->get_fixed_flag(node);
                //if (fixed.first && d != 0.0)
                double *cg_edge_length = m_eMesh->field_data(cg_edge_length_field, node);
                bool moved = d/cg_edge_length[0] > 1.e-6;
                if (moved && fixed.first && (fixed.second == MS_CURVE || fixed.second == MS_ON_BOUNDARY))
                  {
                    //std::cout << "found fixed c/s = " << fixed.second << " curve = " << MS_CURVE << " surface= " << MS_SURFACE << " node= " << m_eMesh->id(node) << std::endl;
                    val[0] = -1;
                  }
                // if (fixed.second == MS_SURFACE)
                //   {
                //     project_delta_to_tangent_plane(node, cg_g);
                //   }
                if (m_eMesh->owned(node))
                  nodes.push_back(node);
              }
          }
        }

      bool continue_outer_loop = true;
      int outer_iter = -1;
      while (continue_outer_loop)
        {
          {
            std::vector< const stk::mesh::FieldBase *> fields;
            fields.push_back(wall_distance_field);
            stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);
            stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
          }

          ++outer_iter;
          for (size_t ii = 0; ii < nodes.size(); ++ii)
            {
              stk::mesh::Entity node = nodes[ii];
              WallDistanceFieldType::value_type *val = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, node);
              int valNew = find_new_value(node, val[0], wall_distance_field, coord_field_orig);

              if (val[0] == 0 || valNew != val[0])
                {
                  if (valNew <= m_nlayers_drop_off && valNew)
                    nodeQ.push(node);
                }
            }
          size_t nqs = nodeQ.size();
          bool did_change = nqs != 0;
          while (nodeQ.size())
            {
              stk::mesh::Entity node = nodeQ.front();
              nodeQ.pop();
              WallDistanceFieldType::value_type *val = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, node);
              int valNew = find_new_value(node, val[0], wall_distance_field, coord_field_orig);
              if (valNew != val[0])
                {
                  val[0] = valNew;
                  //nodeQ.push(node);
                }
            }

          continue_outer_loop = did_change;
          stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & continue_outer_loop ) );
          stk::all_reduce( m_eMesh->get_bulk_data()->parallel() , stk::ReduceMax<1>( & nqs ) );
          if (m_eMesh->get_rank() == 0)
            std::cout << "P[" << m_eMesh->get_rank() << "] outer_iter= " << outer_iter << " nodeQ= " << nqs << " continue_outer_loop= " << continue_outer_loop << std::endl;
        }
    }
  }

  // fills cg_s
  void ReferenceMeshSmootherAlgebraic::get_step()
  {
    PerceptMesh *eMesh = m_eMesh;
    CoordinatesFieldType *coord_field   = static_cast<CoordinatesFieldType*>(eMesh->get_coordinates_field());
    CoordinatesFieldType *coord_field_orig   = static_cast<CoordinatesFieldType*>(m_coord_field_original);
    WallDistanceFieldType *wall_distance_field = static_cast<WallDistanceFieldType *> (eMesh->m_wall_distance_field);
    stk::mesh::FieldBase *cg_s_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_s");
    //stk::mesh::FieldBase *cg_edge_length_field    = eMesh->get_field(stk::topology::NODE_RANK, "cg_edge_length");

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    int spatialDim = eMesh->get_spatial_dim();

    m_scale = 1.e-10;

    // r=0
    eMesh->nodal_field_set_value("cg_s", 0.0);

    std::vector<stk::mesh::Entity> nodes;

    static int anim_step = 0;
    if (m_do_animation)
      {
        std::ostringstream fileid_ss;
        fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

        std::string oname = "step.e";
        if (anim_step > 0) oname += "-s" + fileid_ss.str();
        eMesh->save_as(oname);
        ++anim_step;
      }

    //std::cout << "m_nlayers_drop_off= " << m_nlayers_drop_off << std::endl;
    for (int layers = 1; layers <= m_nlayers_drop_off; ++layers)
      {
        // node loop
        const stk::mesh::BucketVector & buckets = m_eMesh->get_bulk_data()->buckets( m_eMesh->node_rank() );
        for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
          {
            //if (on_locally_owned_part(**k) || on_globally_shared_part(**k))
            if (on_locally_owned_part(**k))
              {
                stk::mesh::Bucket & bucket = **k ;
                const unsigned num_nodes_in_bucket = bucket.size();

                for (unsigned i_node = 0; i_node < num_nodes_in_bucket; i_node++)
                  {
                    stk::mesh::Entity node = bucket[i_node];
                    double *current_coord = stk::mesh::field_data<CoordinatesFieldType>(*coord_field, node);
                    double *orig_coord = stk::mesh::field_data<CoordinatesFieldType>(*coord_field_orig, node);
                    WallDistanceFieldType::value_type *val = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, node);
                    if (val[0] > m_nlayers_drop_off)
                      continue;
                    if (val[0] == -1)
                      continue;
                    std::pair<bool,int> fixed = this->get_fixed_flag(node);
                    if (fixed.first)
                      continue;

                    bool ldebug = false;
                    if (ldebug) std::cout << "node= " << m_eMesh->id(node) << " val[0]= " << val[0] << std::endl;
                    if (val[0] != layers)
                      continue;

                    typedef std::set<stk::mesh::Entity> EntitySet;
                    EntitySet neighbors;
                    m_eMesh->get_node_node_neighbors(node, neighbors);
                    stk::mesh::Entity nodeMin = stk::mesh::Entity();
                    int valMin = std::numeric_limits<int>::max();
                    double dmin = std::numeric_limits<double>::max();
                    for (EntitySet::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
                      {
                        stk::mesh::Entity nnode = *it;
                        if (nnode == node)
                          continue;
                        CoordinatesFieldType::value_type *ncoord = stk::mesh::field_data<CoordinatesFieldType>(*coord_field_orig, nnode);
                        WallDistanceFieldType::value_type *valn = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, nnode);
                        double d = 0.0;
                        for (int jc=0; jc < m_eMesh->get_spatial_dim(); ++jc)
                          {
                            d += (ncoord[jc] - orig_coord[jc])*(ncoord[jc] - orig_coord[jc]);
                          }
                        if (valn[0] <= valMin && d < dmin)
                          {
                            dmin = d;
                            nodeMin = nnode;
                            valMin = valn[0];
                            if (ldebug) std::cout << "000 valMin= " << valMin << " nodeMin= " << m_eMesh->id(nodeMin) << std::endl;
                         }
                        if (ldebug) std::cout << "000 valMin= " << valMin << " nnode= " << m_eMesh->id(nnode) << std::endl;

                        //valMin = std::min(valMin, valn[0]);
                      }
                    if (ldebug) std::cout << "valMin= " << valMin << " nodeMin= " << m_eMesh->id(nodeMin) << std::endl;

                    if (valMin < val[0])
                      {
                        double *n_current_coord = stk::mesh::field_data<CoordinatesFieldType>(*coord_field, nodeMin);
                        double *n_orig_coord = stk::mesh::field_data<CoordinatesFieldType>(*coord_field_orig, nodeMin);
                        //WallDistanceFieldType::value_type *n_val = stk::mesh::field_data<WallDistanceFieldType>(*wall_distance_field, nodeMin);
                        double fac = double(m_nlayers_drop_off + 1 - val[0])/double(m_nlayers_drop_off);
                        double *cg_s = m_eMesh->field_data(*cg_s_field, node);
                        fac = 1.0; // FIXME
                        if (ldebug) std::cout << "fac= " << fac << " val[0]= " << val[0] << std::endl;

                        // here's where to add a drop_off function
                        fac = std::sqrt(std::fabs(fac));
                        for (int jc=0; jc < spatialDim; ++jc)
                          {
                            if (0)
                              {
                                double delta_jc = fac*(n_orig_coord[jc] - orig_coord[jc]);
                                //current_coord[jc] = current_coord[jc] + delta_jc;
                                cg_s[jc] = delta_jc;
                              }
                            else
                              {
                                double delta_jc = fac*(n_current_coord[jc] - n_orig_coord[jc]);
                                //current_coord[jc] = current_coord[jc] + delta_jc;
                                cg_s[jc] = delta_jc;
                              }
                          }


                        // project deltas to surface
                        if (m_eMesh->get_smooth_surfaces())
                          {
                            if (m_eMesh->owned(node))
                              {
                                std::pair<bool,int> fixed2 = this->get_fixed_flag(node);
                                if (!fixed2.first)
                                  {
                                    if (fixed2.second == MS_SURFACE)
                                      {
                                        project_delta_to_tangent_plane(node, cg_s);
                                      }
                                  }
                              }
                          }
                        for (int jc=0; jc < spatialDim; ++jc)
                          {
                            current_coord[jc] = current_coord[jc] + cg_s[jc];
                          }

                      }
                  }
              }
          }

        if (m_do_animation)
          {
            std::ostringstream fileid_ss;
            fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

            std::string oname = "step.e";
            if (anim_step > 0) oname += "-s" + fileid_ss.str();
            eMesh->save_as(oname);
            ++anim_step;
          }

      }


    // project deltas to surface
    if (m_eMesh->get_smooth_surfaces())
      {
        snap_nodes();

        // node loop
        for (size_t ii = 0; ii < nodes.size(); ++ii)
          {
            stk::mesh::Entity node = nodes[ii];
            if (m_eMesh->owned(node))
              {
                std::pair<bool,int> fixed = this->get_fixed_flag(node);
                if (fixed.first)
                  {
                    continue;
                  }

                double *cg_s = m_eMesh->field_data(cg_s_field, node);

                if (fixed.second == MS_SURFACE)
                  {
                    project_delta_to_tangent_plane(node, cg_s);
                  }
              }
          }
      }

    if (m_do_animation)
      {
        std::ostringstream fileid_ss;
        fileid_ss << std::setfill('0') << std::setw(4) << (anim_step);

        std::string oname = "step.e";
        if (anim_step > 0) oname += "-s" + fileid_ss.str();
        eMesh->save_as(oname);
        ++anim_step;
      }

    {
      std::vector< const stk::mesh::FieldBase *> fields;
      fields.push_back(cg_s_field);

      // only the aura = !locally_owned_part && !globally_shared_part (outer layer)
      stk::mesh::communicate_field_data(m_eMesh->get_bulk_data()->aura_ghosting(), fields);
      stk::mesh::copy_owned_to_shared(*eMesh->get_bulk_data(), fields);

      // the shared part (just the shared boundary)
      //stk::mesh::communicate_field_data(*m_eMesh->get_bulk_data()->ghostings()[0], fields);
    }
  }

  double ReferenceMeshSmootherAlgebraic::run_one_iteration()
  {
    PerceptMesh *eMesh = m_eMesh;

    stk::mesh::Selector on_locally_owned_part =  ( eMesh->get_fem_meta_data()->locally_owned_part() );
    stk::mesh::Selector on_globally_shared_part =  ( eMesh->get_fem_meta_data()->globally_shared_part() );
    bool total_valid=false;

    if (m_iter == 0)
      {
        get_edge_lengths(m_eMesh);
        eMesh->nodal_field_set_value("cg_g", 0.0);
        eMesh->nodal_field_set_value("cg_r", 0.0);
        eMesh->nodal_field_set_value("cg_d", 0.0);
        eMesh->nodal_field_set_value("cg_s", 0.0);
      }

    if (m_eMesh->get_rank() == 0)
      std::cout << "ReferenceMeshSmootherAlgebraic: get_wall_distances" << std::endl;
    get_wall_distances();
    if (m_stage == 0 && m_iter== 0)
      eMesh->save_as("wall.e");

    if (m_eMesh->get_rank() == 0)
      std::cout << "ReferenceMeshSmootherAlgebraic: get_step" << std::endl;
    get_step();

    if (m_eMesh->get_rank() == 0)
      std::cout << "ReferenceMeshSmootherAlgebraic: get_step done" << std::endl;

    eMesh->copy_field("cg_g", "cg_s");
    eMesh->nodal_field_axpby(-1.0, "cg_g", 0.0, "cg_r");

    m_dnew = eMesh->nodal_field_dot("cg_s", "cg_s");
    if (m_iter == 0)
      {
        m_d0 = m_dnew;
      }

    /// x = x + alpha*d
    double alpha = 1.0;
    m_alpha = alpha;
    update_node_positions(alpha);

    // check if re-snapped geometry is acceptable
    if (m_eMesh->get_smooth_surfaces())
      {
        snap_nodes();
        if (m_stage != 0)
          {
            bool total_valid_0=true;
            total_metric( 0.0, 1.0, total_valid_0);
            VERIFY_OP_ON(total_valid_0, ==, true, "bad mesh after snap_node_positions...");
          }
      }

    double metric_check = total_metric( 0.0, 1.0, total_valid);
    m_total_metric = metric_check;
    //PRINT_1( "INFO: tmp srk m_iter= " << m_iter << " m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " metric_check= " << metric_check );
    if (check_convergence() || metric_check == 0.0)
      {
        PRINT_1( "INFO: tmp srk already converged m_dnew= " << m_dnew << " gradNorm= " << gradNorm << " metric_check= " << metric_check );
        //update_node_positions
        return total_metric(0.0,1.0, total_valid);
      }

    double tm = total_metric(0.0,1.0, total_valid);
    //PRINT_1( "tmp srk iter= "<< m_iter << " dmax= " << m_dmax << " alpha= " << alpha << " global metric= " << tm);

    return tm;
  }


}


#endif
