// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_IC_Alg.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_DistanceSweeper.hpp>

#include <Akri_CDFEM_Support.hpp>
#include <Akri_AuxMetaData.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_MasterElement.hpp>
#include <Akri_Vec.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_io/IossBridge.hpp>
#include <Akri_MasterElementDeterminer.hpp>

namespace krino{

namespace {
double relative_crossing_position(const double ls0, const double ls1)
{
  return LevelSet::sign_change(ls0, ls1) ? ls0 / ( ls0 - ls1 ) : -10.;
}
}

//----------------------------------------------------------------
void IC_Alg::execute(const double time, const bool requires_additional_initialization)
{ /* %TRACE[ON]% */ Trace trace__("krino::IC_Analytic_Alg::execute()"); /* %TRACE% */

  if (surface_list.size() == 0 && my_calculators.empty())
  {
    return;
  }

  const stk::mesh::BulkData& mesh = levelSet.mesh();
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  const FieldRef xField = levelSet.get_coordinates_field();
  const FieldRef dField = levelSet.get_distance_field();
  const int spatial_dim = meta.spatial_dimension();

  BoundingBox node_bbox;
  levelSet.compute_nodal_bbox( mesh.mesh_meta_data().universal_part(), node_bbox );
  surface_list.prepare_to_compute(time, node_bbox, levelSet.narrow_band_size());

  stk::mesh::BucketVector const& buckets = mesh.get_buckets( stk::topology::NODE_RANK, stk::mesh::selectField(dField) );

  for ( auto && bucket_ptr : buckets )
  {
    const stk::mesh::Bucket & b = *bucket_ptr;
    const int length = b.size();
    double *dist = field_data<double>(dField, b);
    double * coord = field_data<double>(xField, b);

    for (int n = 0; n < length; ++n)
    {
      ThrowAssert(&(dist[n]) != NULL);

      const Vector3d x(&coord[spatial_dim*n], spatial_dim);

      dist[n] = surface_list.point_signed_distance_with_narrow_band(x, levelSet.narrow_band_size());
    }
  }

  stk::mesh::communicate_field_data(mesh, {&dField.field()});

  if (levelSet.narrow_band_size() > 0. && surface_list.truncated_distance_may_have_wrong_sign())
  {
    DistanceSweeper::fix_sign_by_sweeping(mesh, dField, surface_list.get_signed_narrow_band_size(levelSet.narrow_band_size()));
  }

  if(CDFEM_Support::get(meta).get_nonconformal_adapt_target_count() > 0)
  {
    compute_IC_error_indicator();
  }

  for (auto && calc : my_calculators)
  {
    calc->compute_signed_distance(levelSet);
  }

  if (!levelSet.get_keep_IC_surfaces() && !requires_additional_initialization)
  {
    clear();
  }
}

void IC_Alg::compute_IC_error_indicator()
{
  const stk::mesh::BulkData& mesh = levelSet.mesh();
  const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
  const FieldRef xField = levelSet.get_coordinates_field();
  const FieldRef dField = levelSet.get_distance_field();
  stk::mesh::Selector fieldSelector(dField);
  const int spatial_dim = meta.spatial_dimension();
  const auto & cdfem_support = CDFEM_Support::get(meta);

  const auto & aux_meta = AuxMetaData::get(meta);
  const auto & indicator_field_name = cdfem_support.get_nonconformal_adapt_indicator_name();
  auto indicator_field =
      aux_meta.get_field(stk::topology::ELEMENT_RANK, indicator_field_name);
  const auto & elem_buckets =
      mesh.get_buckets(stk::topology::ELEMENT_RANK, meta.locally_owned_part() & fieldSelector);

  std::vector<double> nodal_signed_distances;
  std::vector<Vector3d> nodal_coordinates;
  std::vector<Vector3d> edge_midpoints;
  std::vector<double> midpoint_signed_distances, midpoint_interp_signed_distances;
  int edge_nodes[] = {0, 0};

  for(auto && b_ptr : elem_buckets)
  {
    const auto & b = *b_ptr;

    const int size = b.size();
    double * indicator_data = field_data<double>(indicator_field, b);

    const stk::topology & topo = b.topology();

    const auto & master_elem = MasterElementDeterminer::getMasterElement(*b_ptr, xField);
    const int num_nodes = master_elem.num_nodes();

    nodal_signed_distances.resize(num_nodes);
    nodal_coordinates.resize(num_nodes);

    for(int el=0; el < size; ++el)
    {
      const auto * nodes = mesh.begin_nodes(b[el]);
      for(int n=0; n < num_nodes; ++n)
      {
        auto node = nodes[n];
        nodal_signed_distances[n] = *field_data<double>(dField, node);
        const double * coords_data = field_data<double>(xField, node);
        nodal_coordinates[n] = Vector3d(coords_data, spatial_dim);
      }

      // Iterate edges, find location of crossing on the edge and compare to location of crossing
      // if we add a node at the midpoint with the exact signed distance there.
      const int num_edges = topo.num_edges();
      edge_midpoints.resize(num_edges);
      midpoint_signed_distances.resize(num_edges);
      midpoint_interp_signed_distances.resize(num_edges);
      double err = 0.;
      for (int e=0; e < num_edges; ++e)
      {
        ThrowAssert(topo.edge_topology(e).num_nodes() == 2);
        topo.edge_node_ordinals(e, edge_nodes);

        const Vector3d & x0 = nodal_coordinates[edge_nodes[0]];
        const Vector3d & x1 = nodal_coordinates[edge_nodes[1]];
        edge_midpoints[e] = 0.5*(x0+x1);
        const double edge_length_sqr = (x1-x0).length_squared();

        double midpoint_signed_distance =
            surface_list.point_signed_distance_with_narrow_band(edge_midpoints[e], levelSet.narrow_band_size());
        midpoint_signed_distances[e] = midpoint_signed_distance;

        const double ls0 = nodal_signed_distances[edge_nodes[0]];
        const double ls1 = nodal_signed_distances[edge_nodes[1]];
        midpoint_interp_signed_distances[e] = 0.5 * (ls0 + ls1);
        const double orig_crossing_pos = relative_crossing_position(ls0, ls1);
        const double left_half_crossing_pos =
            0.5 * relative_crossing_position(ls0, midpoint_signed_distance);
        const double right_half_crossing_pos =
            0.5 * (1. + relative_crossing_position(midpoint_signed_distance, ls1));

        if(orig_crossing_pos >= 0.)
        {
          if(left_half_crossing_pos >= 0.)
          {
            const double delta = orig_crossing_pos - left_half_crossing_pos;
            err += edge_length_sqr * delta * delta;
          }
          if(right_half_crossing_pos >= 0.)
          {
            const double delta = orig_crossing_pos - right_half_crossing_pos;
            err += edge_length_sqr * delta * delta;
          }
        }
        else
        {
          if(left_half_crossing_pos >= 0. || right_half_crossing_pos >= 0.)
          {
            err += edge_length_sqr;
          }
        }
      }

      for(int ei=0; ei < num_edges; ++ei)
      {
        for(int ej=ei+1; ej < num_edges; ++ej)
        {
          const double edge_length_sqr = (edge_midpoints[ei] - edge_midpoints[ej]).length_squared();
          const double interp_crossing_pos =
              relative_crossing_position(midpoint_interp_signed_distances[ei],
                  midpoint_interp_signed_distances[ej]);
          const double crossing_pos =
              relative_crossing_position(midpoint_signed_distances[ei],
                  midpoint_signed_distances[ej]);

          const bool has_interp_crossing = interp_crossing_pos >= 0.;
          const bool has_crossing = crossing_pos >= 0.;
          if(has_interp_crossing != has_crossing) err += edge_length_sqr;
          if(has_crossing && has_interp_crossing)
          {
            const double delta = crossing_pos - interp_crossing_pos;
            err += edge_length_sqr * delta * delta;
          }
        }
      }
      indicator_data[el] += err;
    }
  }
}

//----------------------------------------------------------------
} // namespace krino
