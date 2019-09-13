// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#if HAVE_CUBIT

#include <percept/structured/StructuredGridSnap.hpp>
//#include <percept/PerceptUtils.hpp>
//#include <adapt/UniformRefinerPattern.hpp>
//#include <adapt/UniformRefiner.hpp>
#include <sys/time.h>
#include <PGeom.hpp>
#include <percept/structured/PGeomAssocStructured.hpp>
#include <CubitVector.hpp>

namespace percept {

  void StructuredGridSnap::snap_to_geometry(std::shared_ptr<PGeom> pgeom)
  {
    for (unsigned iblock=0; iblock < m_grid->m_sblocks.size(); ++iblock)
    {
        std::shared_ptr<StructuredBlock> sgi = m_grid->m_sblocks[iblock];

        unsigned i_min = sgi->m_sizes.node_min[0];
        unsigned j_min = sgi->m_sizes.node_min[1];
        unsigned k_min = sgi->m_sizes.node_min[2];
        unsigned i_max = sgi->m_sizes.node_max[0];
        unsigned j_max = sgi->m_sizes.node_max[1];
        unsigned k_max = sgi->m_sizes.node_max[2];

        // snap minimum j block face to the outer cylinder surface
        std::array<unsigned, 3> range_min = {{i_min, j_min, k_min}};
        std::array<unsigned, 3> range_max = {{i_max, j_min, k_max}};
        int outer_cyl_id = 5;
        snap_to_surface(pgeom, sgi, range_min, range_max, outer_cyl_id);

        // snap maximum j block face to the inner cylinder surface
        range_min[1] = j_max;
        range_max[1] = j_max;
        int inner_cyl_id = 4;
        snap_to_surface(pgeom, sgi, range_min, range_max, inner_cyl_id);


    }
  }

  void StructuredGridSnap::snap_to_geometry(std::shared_ptr<PGeom> pgeom,
                                            PGeomAssocStructured &passoc)
  {
    for (unsigned iblock=0; iblock < m_grid->m_sblocks.size(); ++iblock)
    {
        std::shared_ptr<StructuredBlock> sgi = m_grid->m_sblocks[iblock];

        std::map<uint64_t, NodeAssociatedGeom> &block_node_map = passoc.get_refined_node_map(iblock);


        for (std::pair<const int, NodeAssociatedGeom> node_data : block_node_map)
        {
          uint64_t local_node_id = node_data.first;
          NodeAssociatedGeom &node_geom = node_data.second;

          std::array<unsigned, 3> node_ijk;
          sgi->multi_dim_indices_from_local_offset(local_node_id, node_ijk);

          if (node_geom.dimension == 1)
          {
            // TODO - handle multiple associations
            int curve_id = node_geom.geom_ids[0];
            snap_to_curve(pgeom, sgi, node_ijk, node_ijk, curve_id);
          }
          else if (node_geom.dimension == 2)
          {
            int surf_id = node_geom.geom_ids[0];
            snap_to_surface(pgeom, sgi, node_ijk, node_ijk, surf_id);
          }
        }
    }
  }

  void StructuredGridSnap::snap_to_surface(std::shared_ptr<PGeom> pgeom,
                       std::shared_ptr<StructuredBlock> sgi,
                       const std::array<unsigned, 3> &range_min,
                       const std::array<unsigned, 3> &range_max,
                       const int surf_id)
  {
    std::array<unsigned,3> indx{{0,0,0}};

    for (indx[2] = range_min[2]; indx[2] <= range_max[2]; ++indx[2]) // K
    {
        for (indx[1] = range_min[1]; indx[1] <= range_max[1]; ++indx[1]) // J
        {
            for (indx[0] = range_min[0]; indx[0] <= range_max[0]; ++indx[0]) // I
            {
                double coords[3];
                for (unsigned ic = 0; ic < 3; ++ic)
                {
                    coords[ic] = sgi->m_sgrid_coords(indx[0], indx[1], indx[2], ic);
                }

                CubitVector original_pt(coords);
                CubitVector projected_pt;

                // snap coordinates to geometry
                projected_pt = pgeom->get_surf_closest_point(surf_id, original_pt);

                for (unsigned ic = 0; ic < 3; ++ic)
                {
                    sgi->m_sgrid_coords(indx[0], indx[1], indx[2], ic) = projected_pt[ic];
                }
            }
        }
    }
  }

  void StructuredGridSnap::snap_to_curve(std::shared_ptr<PGeom> pgeom,
                       std::shared_ptr<StructuredBlock> sgi,
                       const std::array<unsigned, 3> &range_min,
                       const std::array<unsigned, 3> &range_max,
                       const int curve_id)
  {
    std::array<unsigned,3> indx{{0,0,0}};

    for (indx[2] = range_min[2]; indx[2] <= range_max[2]; ++indx[2]) // K
      {
        for (indx[1] = range_min[1]; indx[1] <= range_max[1]; ++indx[1]) // J
          {
            for (indx[0] = range_min[0]; indx[0] <= range_max[0]; ++indx[0]) // I
              {
                double coords[3];
                for (unsigned ic = 0; ic < 3; ++ic)
                  {
                    coords[ic] = sgi->m_sgrid_coords(indx[0], indx[1], indx[2], ic);
                  }

                CubitVector original_pt(coords);
                CubitVector projected_pt;

                // snap coordinates to geometry
                projected_pt = pgeom->get_curv_closest_point(curve_id, original_pt);

                for (unsigned ic = 0; ic < 3; ++ic)
                {
                    sgi->m_sgrid_coords(indx[0], indx[1], indx[2], ic) = projected_pt[ic];
                }
              }
          }
      }
  }


}

#endif
