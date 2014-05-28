/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <assert.h>
#include <iostream>
#include <transfer/UseCaseIsInElement.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/base/FieldData.hpp>

namespace stk {
namespace usecase {

size_t is_in_element(VectorField *domain_coordinates,
                     VectorField *range_coordinates,
                     const std::vector<std::pair<stk::mesh::Entity*, stk::mesh::Entity*> > &entity_map,
                     std::vector<std::size_t> &not_in_element)
{
  // For each pair of mesh entities, determine if the second is contained in the first.
  const std::size_t num_entities = entity_map.size();

  for (std::size_t i = 0; i < num_entities; ++i) {
    stk::mesh::Entity* domain_entity = entity_map[i].first;
    assert(domain_entity->entity_rank() > stk::mesh::fem::FEMMetaData::NODE_RANK);

    stk::mesh::Entity* range_entity  = entity_map[i].second;
    assert(range_entity->entity_rank() == stk::mesh::fem::FEMMetaData::NODE_RANK);

    const CellTopologyData * const cell_topo = stk::mesh::fem::get_cell_topology(*domain_entity).getCellTopologyData();
    const int nodes_per_entity = cell_topo->node_count;

    const stk::mesh::PairIterRelation entity_nodes = domain_entity->relations(stk::mesh::fem::FEMMetaData::NODE_RANK);
    double *domain_fld_data = static_cast<double*>(stk::mesh::field_data(*domain_coordinates, *entity_nodes[0].entity()));
    assert(domain_fld_data != NULL);
    double xmin = domain_fld_data[0];
    double ymin = domain_fld_data[1];
    double zmin = domain_fld_data[2];

    double xmax = domain_fld_data[0];
    double ymax = domain_fld_data[1];
    double zmax = domain_fld_data[2];

    for (int j = 1; j < nodes_per_entity; ++j) {
      domain_fld_data = static_cast<double*>(stk::mesh::field_data(*domain_coordinates, *entity_nodes[j].entity()));
      assert(domain_fld_data != NULL);
      xmin = domain_fld_data[0] < xmin ? domain_fld_data[0] : xmin;
      ymin = domain_fld_data[1] < ymin ? domain_fld_data[1] : ymin;
      zmin = domain_fld_data[2] < zmin ? domain_fld_data[2] : zmin;

      xmax = domain_fld_data[0] > xmax ? domain_fld_data[0] : xmax;
      ymax = domain_fld_data[1] > ymax ? domain_fld_data[1] : ymax;
      zmax = domain_fld_data[2] > zmax ? domain_fld_data[2] : zmax;
    }
    double *range_fld_data = static_cast<double*>(stk::mesh::field_data(*range_coordinates, *range_entity));
    assert(range_fld_data != NULL);
    const double x=range_fld_data[0];
    const double y=range_fld_data[1];
    const double z=range_fld_data[2];
    if (x<xmin || xmax<x || y<ymin || ymax<y || z<zmin || zmax<z) {
      not_in_element.push_back(i);
    }
  }
  return not_in_element.size();
}

}
}
