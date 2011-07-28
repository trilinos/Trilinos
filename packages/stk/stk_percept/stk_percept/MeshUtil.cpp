/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_percept/MeshUtil.hpp>

namespace stk {
  namespace percept {

    bool MeshUtil::m_debug = false;

    // ================================================================================================================================================================
    // ================================================================================================================================================================
    // ================================================================================================================================================================

    void MeshUtil::fillSideNodes(stk::mesh::Entity& element, unsigned iside, std::vector<stk::mesh::EntityId>& side_nodes)
    {
      CellTopology cell_topo(stk::percept::PerceptMesh::get_cell_topology(element));
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      int nfn = cell_topo.getCellTopologyData()->side[iside].topology->vertex_count;
      side_nodes.resize(nfn);
      for (int inode = 0; inode < nfn; inode++)
        {
          int jnode = cell_topo.getCellTopologyData()->side[iside].node[inode];
          side_nodes[inode] = elem_nodes[jnode].entity()->identifier();
        }
    }

    void MeshUtil::fillSideNodes(stk::mesh::Entity& element, unsigned iside, std::vector<stk::mesh::Entity *>& side_nodes)
    {
      CellTopology cell_topo(stk::percept::PerceptMesh::get_cell_topology(element));
      const mesh::PairIterRelation elem_nodes = element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK);

      int nfn = cell_topo.getCellTopologyData()->side[iside].topology->vertex_count;
      side_nodes.resize(nfn);
      for (int inode = 0; inode < nfn; inode++)
        {
          int jnode = cell_topo.getCellTopologyData()->side[iside].node[inode];
          side_nodes[inode] = elem_nodes[jnode].entity();
        }
    }

    double MeshUtil::triFaceArea(percept::PerceptMesh& eMesh, stk::mesh::Entity& element, unsigned iside)
    {
      std::vector<stk::mesh::Entity *> side_nodes;
      fillSideNodes(element, iside, side_nodes);

      //int spatialDim = eMesh.getSpatialDim();
      double a[3]={0,0,0};
      double b[3]={0,0,0};
      double c[3]={0,0,0};

      double *fdata0 = stk::mesh::field_data( *eMesh.getCoordinatesField() , *side_nodes[0]);
      double *fdata1 = stk::mesh::field_data( *eMesh.getCoordinatesField() , *side_nodes[1]);
      double *fdata2 = stk::mesh::field_data( *eMesh.getCoordinatesField() , *side_nodes[2]);
      for (int idim=0; idim < 3; idim++)
        {
          a[idim] = fdata1[idim] - fdata0[idim];
          b[idim] = fdata2[idim] - fdata0[idim];
        }
      c[0] = a[1]*b[2] - a[2]*b[1];
      c[1] = -(a[0]*b[2] - a[2]*b[0]);
      c[2] = a[0]*b[1] - a[1]*b[0];

      return 0.5*std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
    }

    bool MeshUtil::nodesMatch(  std::vector<stk::mesh::EntityId>& side1, std::vector<stk::mesh::EntityId>& side2, bool reverse)
    {
      if (side1.size() != side2.size()) return false;
      //std::cout << "tmp side1= " << side1 << " side2= " << side2 << std::endl;
      int nfn = side1.size();
      for (int i = 0; i < nfn; i++)
        {
          bool found = true;
          if (reverse)
            {
              int k=0;
              for (int j = nfn-1; j >= 0; --j)
                {
                  if (side1[k++] != side2[(i+j) % nfn])
                    {
                      found = false;
                      break;
                    }
                }
            }
          else
            {
              for (int j = 0; j < nfn; j++)
                {
                  if (side1[j] != side2[(i+j) % nfn])
                    {
                      found = false;
                      break;
                    }
                }
            }

          if (found)
            {
              return true;
            }
        }
      return false;
    }

    bool MeshUtil::sharesFace(stk::mesh::Entity& element1, stk::mesh::Entity& element2, unsigned& iside1, unsigned& iside2)
    {
      CellTopology cell_topo1(stk::percept::PerceptMesh::get_cell_topology(element1));
      CellTopology cell_topo2(stk::percept::PerceptMesh::get_cell_topology(element2));
      if (cell_topo1.getKey() != cell_topo2.getKey())
        return false;
        
      std::vector<stk::mesh::EntityId> side1;
      std::vector<stk::mesh::EntityId> side2;
      unsigned nsides = (unsigned)cell_topo1.getSideCount();
      for (iside1 = 0; iside1 < nsides; iside1++)
        {
          fillSideNodes(element1, iside1, side1);
          for (iside2 = 0; iside2 < nsides; iside2++)
            {
              fillSideNodes(element2, iside2, side2);

              // true = reverse the nodes since we are looking at a face from two sides
              if (nodesMatch(side1, side2, true))
                {
                  return true;
                }
                
            }
        }
      return false;
    }

    bool MeshUtil::facesConsistent1(percept::PerceptMesh& eMesh, stk::mesh::Entity& element1, stk::mesh::Entity& element2)
    {
      //int spatialDim = eMesh.getSpatialDim();
      //unsigned side_rank = (spatialDim == 3 ? m_eMesh.face_rank() : m_eMesh.edge_rank());
      double tol = 1.e-5;

      CellTopology cell_topo1(stk::percept::PerceptMesh::get_cell_topology(element1));
      CellTopology cell_topo2(stk::percept::PerceptMesh::get_cell_topology(element2));

      if (cell_topo1.getKey() != cell_topo2.getKey())
        return false;

      unsigned ichild_side1=0, ichild_side2=0;
      unsigned iside1=0, iside2=0;
      if (eMesh.isParentElement(element1, false) && eMesh.isParentElement(element2, false))
        {
          if (sharesFace(element1, element2, iside1, iside2))
            {
              double ptriFaceArea1 = triFaceArea(eMesh, element1, iside1);
              double ptriFaceArea2 = triFaceArea(eMesh, element2, iside2);
              if (std::abs(ptriFaceArea1 - ptriFaceArea2) > tol*(ptriFaceArea1 + ptriFaceArea2)/2.)
                {
                  throw std::runtime_error("triFaceArea inconsistent parent");
                }
              std::vector<stk::mesh::Entity *> children1, children2;
              bool hasChildren1 = eMesh.getChildren(element1, children1, true, true);  // check_for_family_tree, only_if_element_is_parent_leaf
              bool hasChildren2 = eMesh.getChildren(element2, children2, true, true);
              VERIFY_OP_ON(hasChildren1, &&, hasChildren2, "no children");
              double cfaceTot=0.0;

              for (unsigned ichild1 = 0; ichild1 < children1.size(); ichild1++)
                {
                  for (unsigned ichild2 = 0; ichild2 < children2.size(); ichild2++)
                    {
                      stk::mesh::Entity& child1 = *children1[ichild1];
                      stk::mesh::Entity& child2 = *children2[ichild2];
                      if (sharesFace(child1, child2, ichild_side1, ichild_side2))
                        {
                          double ctriFaceArea1 = triFaceArea(eMesh, child1, ichild_side1);
                          double ctriFaceArea2 = triFaceArea(eMesh, child2, ichild_side2);
                          if (std::abs(ctriFaceArea1 - ctriFaceArea2) > tol*(ctriFaceArea1 + ctriFaceArea2)/2.)
                            {
                              throw std::runtime_error("triFaceArea inconsistent child");
                            }
                          cfaceTot += ctriFaceArea1;
                        }
                    }
                }
              if (m_debug)
                {
                  std::cout << "E1= " << element1.identifier() << " E2= " << element2.identifier() << " cfaceTot= " << cfaceTot 
                            << " ptriFaceArea1= " << ptriFaceArea1 << std::endl;
                }
              if (std::abs(cfaceTot - ptriFaceArea1) > tol*(cfaceTot+ptriFaceArea1)/2.)
                {
                  throw std::runtime_error("triFaceArea inconsistent child/parent");
                }
            }
        }
      return true;
    }


    bool MeshUtil::facesConsistent(percept::PerceptMesh& eMesh)
    {
      const std::vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( eMesh.element_rank() );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& element = bucket[iElement];

              for ( std::vector<stk::mesh::Bucket*>::const_iterator k1 = buckets.begin() ; k1 != buckets.end() ; ++k1 ) 
                {
                  stk::mesh::Bucket & bucket1 = **k1 ;

                  const unsigned num_elements_in_bucket1 = bucket1.size();
                  for (unsigned iElement1 = 0; iElement1 < num_elements_in_bucket1; iElement1++)
                    {
                      stk::mesh::Entity& element1 = bucket1[iElement1];
                      if (&element1 != &element)
                        {
                          bool isConsistent = facesConsistent1(eMesh, element, element1);
                          if (!isConsistent)
                            return false;
                        }
                    }
                }
            }
        }
      return true;
    }

    void MeshUtil::checkTopology(percept::PerceptMesh& eMesh)
    {
      const std::vector<stk::mesh::Bucket*> & buckets = eMesh.getBulkData()->buckets( eMesh.element_rank() );

      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) 
        {
          stk::mesh::Bucket & bucket = **k ;

          const unsigned num_elements_in_bucket = bucket.size();
          for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
            {
              stk::mesh::Entity& element = bucket[iElement];
              if (!stk::percept::PerceptMesh::get_cell_topology(element))
                {
                  // an empty element
                  if (element.relations(stk::mesh::fem::FEMMetaData::NODE_RANK).size() == 0)
                    continue;

                  std::cout << "Error MeshUtil::checkTopology null" << std::endl;
                  eMesh.printEntity(std::cout, element);
                  throw std::logic_error("Error MeshUtil::checkTopology null" );
                }
            }
        }
    }
  }
}

