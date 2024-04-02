// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/DihedralAngleCheck.hpp>
#include <percept/mesh/mod/smoother/JacobianUtil.hpp>

#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/Refiner.hpp>
#include <adapt/IEdgeBasedAdapterPredicate.hpp>
#include <adapt/PredicateBasedEdgeAdapter.hpp>

#include <utility>
#include <vector>

namespace percept {

  void
  DihedralAngleCheck::find_simplex_elements_with_obtuse_angles(std::ostream& out)
  {
    PerceptMesh& mesh = *m_eMesh;
    m_map.clear();

    //std::cerr << "find_simplex_elements_with_obtuse_angles..." << std::endl;
    double max_angle = 0.;

    // Skip 2D for now
    if (mesh.get_spatial_dim() != 3) {
      return;
    }

    constexpr int spatial_dim = 3;
    const stk::topology simplex_topology = spatial_dim == 2 ?
      stk::topology::TRIANGLE_3_2D :
      stk::topology::TETRAHEDRON_4 ;

    stk::mesh::Selector active_owned_selector =
      PerceptMesh::select_active_elements(*mesh.get_bulk_data(), std::vector<stk::mesh::EntityRank>(1, mesh.element_rank()))
      & mesh.get_fem_meta_data()->locally_owned_part()
      & mesh.get_fem_meta_data()->get_topology_root_part(simplex_topology);

    size_t num_bad = 0;
    JacobianUtil jac;
    DenseMatrix<3,3> mat;
    double detJ=0;
    int side_node_ordinals[4][3];
    for (auto && b_ptr : mesh.get_bulk_data()->get_buckets(stk::topology::ELEMENT_RANK, active_owned_selector))
      {
        const auto & topology = b_ptr->topology();
        const int num_nodes = topology.num_nodes();
        const int num_sides = topology.num_sides();
        const int num_side_nodes = topology.side_topology().num_nodes();

        const double *nodal_coords[num_nodes];
        const double *side_nodal_coords[3];
        double side_normals[num_sides][spatial_dim];

        //auto surface_topology = topology.side_topology();
        for (int s = 0; s < num_sides; ++s)
          {
            topology.side_node_ordinals(s, &side_node_ordinals[s][0]);
            if (0) std::cout << "sn[" << s << "] = "
                      << side_node_ordinals[s][0] << " "
                      << side_node_ordinals[s][1] << " "
                      << side_node_ordinals[s][2] << std::endl;
          }

        for (auto && elem : *b_ptr)
          {
            const auto * nodes = mesh.get_bulk_data()->begin(elem, stk::topology::NODE_RANK);
            for (int i = 0; i < num_nodes; ++i)
              {
                nodal_coords[i] = static_cast<double*>(stk::mesh::field_data(*mesh.get_coordinates_field(), nodes[i]));
              }

            // Determine normal vectors of each side
            for (int s = 0; s < num_sides; ++s)
              {
                for (int sn = 0; sn < num_side_nodes; ++sn)
                  {
                    const int node_ordinal = side_node_ordinals[s][sn];
                    side_nodal_coords[sn] = nodal_coords[node_ordinal];
                  }
                jac.jacobian_matrix_2D_in_3D(detJ, mat, side_nodal_coords, &side_normals[s][0]);
                double vl = Math::norm_3d(&side_normals[s][0]);
                if (vl == 0.0) vl = 1.0;
                for (int d = 0; d < spatial_dim; ++d)
                  side_normals[s][d] /= vl;

                if (0) std::cout << "side_normals[" << s << "] = " << Math::print_3d(&side_normals[s][0])
                          << " side_nodal_coords= \n"
                          << Math::print_3d(side_nodal_coords[0]) << "\n"
                          << Math::print_3d(side_nodal_coords[1]) << "\n"
                          << Math::print_3d(side_nodal_coords[2]) << std::endl;
              }

            // Need to iterate adjacent pairs of sides and determine dihedral angles
            EntityObtuseAngleSidePairs obtuse_angle_pairs(elem, {});
            for (int s1 = 0; s1 < num_sides - 1; ++s1)
              {
                for (int s2 = s1 + 1; s2 < num_sides; ++s2)
                  {
                    double dot = 0.;
                    for (int d = 0; d < spatial_dim; ++d)
                      {
                        dot += side_normals[s1][d] * side_normals[s2][d];
                      }
                    const double angle = std::acos(-dot) * 180. / M_PI;
                    if (0) out << mesh.rank() << " DihedralAngleCheck: Elem " << mesh.identifier(elem) << " dihedral angle " << angle << " dot= " << dot << " between faces " << s1 << " and " << s2 << std::endl;;
                    if (angle >= 90.)
                      {
                        obtuse_angle_pairs.second.emplace_back(s1, s2);
                        ++num_bad;
                        if (m_print_level > 2) out << mesh.rank() << " DihedralAngleCheck: Elem " << mesh.identifier(elem) << " has obtuse dihedral angle " << angle << " between faces " << s1 << " and " << s2 << std::endl;;
                      }
                    max_angle = std::max(max_angle, angle);
                  }
              }
            if (!obtuse_angle_pairs.second.empty())
              {
                m_map[elem] = obtuse_angle_pairs.second;
              }
          }
      }

    stk::ParallelMachine pm = mesh.get_bulk_data()->parallel();
    double max_angle_global = max_angle;
    stk::all_reduce( pm, stk::ReduceMax<1>( &max_angle_global ) );
    stk::all_reduce( pm, stk::ReduceSum<1>( &num_bad ) );

    if (m_print_level > 0 && mesh.get_rank() == 0)
      out << mesh.rank() << " DihedralAngleCheck: max_angle= " << max_angle_global << " number of obtuse angles= " << num_bad << std::endl;


  }



}
