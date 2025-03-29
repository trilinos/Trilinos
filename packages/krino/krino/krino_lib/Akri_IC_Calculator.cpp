// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_IC_Calculator.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_FieldRef.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_Surface_Manager.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <limits>

namespace krino{

void IC_Binder::compute_signed_distance(const LevelSet &ls) const
{
  const stk::mesh::MetaData & meta = ls.meta();
  const stk::mesh::BulkData & mesh = ls.mesh();
  const FieldRef distRef = ls.get_isovar_field();

  const stk::mesh::Selector selector = stk::mesh::selectField(distRef);

  const krino::Surface_Manager & surfaceManager = krino::Surface_Manager::get(meta);

  std::vector<FieldRef> otherDistRefs;
  for (auto && otherLS : surfaceManager.get_levelsets())
  {
    if(ls.get_identifier() == otherLS->get_identifier())
    {
      break;
    }
    otherDistRefs.push_back(otherLS->get_isovar_field());
  }

  const unsigned num_other_dist = otherDistRefs.size();
  std::vector<double *> otherDist(num_other_dist, nullptr);
  for ( auto && bucket : mesh.get_buckets( stk::topology::NODE_RANK, selector) )
  {
    const stk::mesh::Bucket & b = *bucket;
    const size_t length = b.size();

    double * dist = field_data<double>(distRef, b);
    for (unsigned i=0; i<num_other_dist; ++i)
    {
      otherDist[i] = field_data<double>(otherDistRefs[i], b);
    }

    for ( size_t n = 0; n < length; ++n )
    {
      dist[n] = std::numeric_limits<double>::max();
      for (unsigned i=0; i<(num_other_dist-my_binder_type); ++i)
      {
        for (unsigned j=i+1; j<(num_other_dist-my_binder_type); ++j)
        {
          //ensure binder separates particles by a distance of my_interface_size
          if (my_binder_type==0)
          {
            dist[n] = std::min(dist[n], std::max(otherDist[i][n]-my_interface_size, otherDist[j][n]-my_interface_size));
          }

          //create smooth connection between particles if bridge size specified
          //take the product of two offset particle level sets and move the zero crossing by bridge_size
          if (my_binder_type==1)
          {
            if (my_smooth_bridge_size > 0.0)
            {
              double bridge = (otherDist[i][n]+my_smooth_bridge_offset)*(otherDist[j][n]+my_smooth_bridge_offset)-my_smooth_bridge_size;
              if (my_root_smooth_bridge)
              {
                //root the level set product to bring it back to distance function magnitudes
                const int sign = (bridge >= 0.0) ? 1 : -1;
                bridge = sign*std::sqrt(std::abs(bridge));
              }
              dist[n] = std::min(dist[n], bridge);
            }
          }
        }
      }
      for (unsigned i=0; i<num_other_dist; ++i)
      {
        if (my_binder_type==1)
        {
          otherDist[i][n] *= my_other_ls_scale_factor;
        }
      }
    }
  }
}

} // namespace krino
