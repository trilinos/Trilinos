// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ExcludeWedgesNotConnectedToPyramids_hpp
#define percept_ExcludeWedgesNotConnectedToPyramids_hpp

#include <percept/PerceptMesh.hpp>

namespace percept
{

void excludeWedgesNotConnectedToPyramids(PerceptMesh& eMesh, stk::mesh::PartVector& exclude_part_vector)
{
  std::set<stk::mesh::Entity> exclude_wedges_set;
  
  // init set of all wedges
  const stk::mesh::BucketVector & buckets = eMesh.get_bulk_data()->buckets( stk::topology::ELEMENT_RANK );
  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      stk::mesh::Bucket & bucket = **k ;
      if (bucket.owned() && bucket.topology() == stk::topology::WEDGE_6)
        {
          const unsigned num_elems_in_bucket = bucket.size();
          for (unsigned iWedge = 0; iWedge < num_elems_in_bucket; iWedge++)
            {
              stk::mesh::Entity wedge = bucket[iWedge];
              exclude_wedges_set.insert(wedge);
            }
        }
    }
  
  std::cout << "init wedges excluded: " << exclude_wedges_set.size() << std::endl;
  
  typedef std::pair<stk::mesh::EntityId,stk::mesh::EntityId> Edge;
  std::set<Edge> all_edges_to_remove;
  
  // remove wedges adjacent to a pyramid - put them in old_wedges_front
  for ( stk::mesh::BucketVector::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
    {
      stk::mesh::Bucket & bucket = **k ;
      if (bucket.owned() && bucket.topology() == stk::topology::PYRAMID_5)
        {
          const unsigned num_elems_in_bucket = bucket.size();
          for (unsigned iPyr = 0; iPyr < num_elems_in_bucket; iPyr++)
            {
              stk::mesh::Entity pyramid = bucket[iPyr];
              stk::mesh::Entity neigh = eMesh.get_face_neighbor(pyramid, 4); // guess that face 4 is the quad face
              if (eMesh.is_valid(neigh) &&
                  eMesh.topology(neigh) == stk::topology::WEDGE_6) {
                exclude_wedges_set.erase(neigh);
                stk::mesh::EntityId pyramid_id = eMesh.identifier(pyramid);
                stk::mesh::EntityId neigh_id = eMesh.identifier(neigh);
                Edge edge(pyramid_id,neigh_id);
                all_edges_to_remove.insert(edge);
              }
            }
        }
    }
  
  std::cout << "found num wedges adjacent to pyramids: " << all_edges_to_remove.size() << std::endl;
  std::cout << "current wedges excluded: " << exclude_wedges_set.size() << std::endl;
  
  std::set<Edge> old_edges_to_remove(all_edges_to_remove.begin(), all_edges_to_remove.end());
  
  unsigned num_changed=old_edges_to_remove.size();
  int counter=0;
  while (num_changed && counter++<2000) {
    
    std::set<Edge> new_edges_to_remove;
    
    for (std::set<Edge>::iterator iter = old_edges_to_remove.begin();
         iter != old_edges_to_remove.end(); ++iter)
      {
        stk::mesh::EntityId wedge1_id = iter->second;
        stk::mesh::Entity wedge1 = eMesh.get_element(wedge1_id);
        for (unsigned face=0; face<3; face++) { // assuming faces 0,1,2 are quad faces
          stk::mesh::Entity neigh = eMesh.get_face_neighbor(wedge1, face);
          stk::mesh::EntityId neigh_id = eMesh.identifier(neigh);
          
          Edge new_edge = Edge(wedge1_id,neigh_id); // point from me to neigh
          Edge old_edge = Edge(neigh_id,wedge1_id); // reversed edge from neigh
          
          if (eMesh.is_valid(neigh) &&
              eMesh.topology(neigh) == stk::topology::WEDGE_6 &&
              all_edges_to_remove.find(old_edge) == all_edges_to_remove.end()) {
            
            new_edges_to_remove.insert(new_edge);
            all_edges_to_remove.insert(new_edge);
            exclude_wedges_set.erase(neigh);
          }
        }
      }
    num_changed=new_edges_to_remove.size();
    
    old_edges_to_remove.clear();
    old_edges_to_remove.insert(new_edges_to_remove.begin(), new_edges_to_remove.end());
    
    std::cout << "counter = " << counter
              << " num_changed = " << num_changed << std::endl;
  }
  
  std::cout << "final wedges excluded: " << exclude_wedges_set.size() << std::endl;
  
  std::vector<stk::mesh::Entity> exclude_wedges(exclude_wedges_set.begin(), exclude_wedges_set.end());
  
  std::vector<stk::mesh::PartVector> add_parts(exclude_wedges.size(), exclude_part_vector),
            remove_parts(exclude_wedges.size());
  
  eMesh.get_bulk_data()->batch_change_entity_parts(exclude_wedges, add_parts, remove_parts);          
}
  
}

#endif
