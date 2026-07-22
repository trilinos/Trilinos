// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshDiagnostics.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_Element.hpp>
#include <Akri_MasterElement.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <Akri_MasterElementDeterminer.hpp>

namespace krino{

void print_volume_or_surface_area(const stk::mesh::BulkData& mesh, const stk::mesh::EntityRank entity_rank, const stk::mesh::Selector & active_selector, const stk::mesh::PartVector & parts)
{
  double overall_sum = 0.;

  const std::string label = (entity_rank == stk::topology::ELEMENT_RANK) ? "Volume" : "Surface Area";

  krinolog << stk::diag::dendl;

  for (auto && part : parts)
  {
    const stk::mesh::Selector selector = mesh.mesh_meta_data().locally_owned_part() & active_selector & *part;
    const double local_part_sum = compute_volume_or_surface_area(mesh, entity_rank, selector);
    double global_part_sum = 0.;

    stk::all_reduce_sum(mesh.parallel(), &local_part_sum, &global_part_sum, 1);
    overall_sum += global_part_sum;

    krinolog << label << " for part " << part->name() << " = " << global_part_sum << stk::diag::dendl;
  }

  krinolog << "Sum of " << label << " for all parts (which may overlap) = " << overall_sum << stk::diag::dendl;
}

double
compute_volume_or_surface_area(const stk::mesh::BulkData& mesh, const stk::mesh::EntityRank entity_rank, const stk::mesh::Selector & selector)
{
  const FieldRef coords_field(mesh.mesh_meta_data().coordinate_field());
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();

  std::vector<double> coords;
  std::vector<double> intg_weights;

  double volume_or_area = 0.;
  for ( auto && bucket : mesh.get_buckets( entity_rank, selector ) )
  {
    const krino::MasterElement& master_elem = MasterElementDeterminer::getMasterElement(bucket->topology());

    for ( auto && entity : *bucket )
    {
      ElementObj::gather_nodal_field(mesh, entity, coords_field, coords, dim);
      ElementObj::integration_weights( intg_weights, dim, coords, master_elem );

      for ( double intg_wt : intg_weights ) volume_or_area += intg_wt;
    }
  }
  return volume_or_area;
}

} // namespace krino
