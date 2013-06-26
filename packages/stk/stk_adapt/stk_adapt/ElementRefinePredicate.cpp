
#include <stk_adapt/ElementRefinePredicate.hpp>

namespace stk {
  namespace adapt {

    /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
    int ElementRefinePredicate::operator()(const stk::mesh::Entity entity)
    {
      double *fdata = 0;
      if (m_field)
        fdata = m_eMesh.field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
      bool selected = (m_eb_selector==0 || (*m_eb_selector)(m_eMesh.bucket(entity)));
      bool ref_field_criterion = (fdata  && fdata[0] > 0);
      bool unref_field_criterion = (fdata && fdata[0] < 0);
      int mark = 0;
      if (selected && ref_field_criterion) mark |= DO_REFINE;
      if (selected && unref_field_criterion) mark |= DO_UNREFINE;
      return mark;
    }

      // void check_two_to_one(PerceptMesh& eMesh);
      // void enforce_two_to_one_refine(PerceptMesh& eMesh);
      // void ok_to_unrefine(PerceptMesh& eMesh);


    bool ElementRefinePredicate::is_face_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level)
    {
      stk::mesh::Entity parent_element = element;
      stk::mesh::Entity parent_neighbor = neighbor;
      if (element_level < neighbor_level)
        {
          parent_element = m_eMesh.getParent(element, true);
        }
      else if (element_level > neighbor_level)
        {
          parent_neighbor = m_eMesh.getParent(neighbor, true);
        }
      bool fn = m_eMesh.is_face_neighbor(parent_element, parent_neighbor);
      return fn;
    }

    bool ElementRefinePredicate::check_two_to_one()
    {
      ScalarIntFieldType *refine_level = m_eMesh.get_fem_meta_data()->get_field<ScalarIntFieldType>("refine_level");
      if (!refine_level)
        {
          throw std::logic_error("must have refine_level field for hanging-node refinement");
        }
      bool valid = true;
      std::set<stk::mesh::Entity> neighbors;
      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
      const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  neighbors.clear();
                  m_eMesh.get_node_neighbors(element, neighbors);
                  for (std::set<stk::mesh::Entity>::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); neighbor++)
                    {
                      int *ref_lev_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                      int *ref_lev_neigh = m_eMesh.get_bulk_data()->field_data( *refine_level , *neighbor );
                      if (std::abs(ref_lev_neigh[0] - ref_lev_elem[0]) > 1)
                        {
                          bool fn = is_face_neighbor(element, ref_lev_elem[0], *neighbor, ref_lev_neigh[0]);
                          if (fn)
                            {
                              valid = false;
                            }
                        }

                    }
                }
            }
        }

      return valid;
    }


  }
}

