
#include <stk_adapt/ElementRefinePredicate.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace stk {
  namespace adapt {

    static bool debug_print = false;

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

    void ElementRefinePredicate::refine(IAdapter& breaker, bool enforce_what[3])
    {
      int max_iter=100;
      int iter=0;
      bool did_change=false;
      while ((iter++ < max_iter) && (did_change = this->enforce_two_to_one_refine(enforce_what)) )
        {
          if (debug_print)
            std::cout << "P[" << m_eMesh.get_rank() << " iter= " << iter << " did_change= " << did_change
                      << std::endl;
        }
      breaker.doBreak();
    }

    void ElementRefinePredicate::unrefine(IAdapter& breaker, bool enforce_what[3])
    {
      int max_iter=100;
      int iter=0;
      bool did_change=false;
      while ((iter++ < max_iter) && (did_change = this->enforce_two_to_one_unrefine(enforce_what)) )
        {
          if (debug_print)
            std::cout << "P[" << m_eMesh.get_rank() << "]  iter= " << iter
                      << " did_change= " << did_change
                      << std::endl;
        }

      ElementUnrefineCollection elements_to_unref = breaker.buildUnrefineList();
      breaker.unrefineTheseElements(elements_to_unref);
    }

    bool ElementRefinePredicate::enforce_two_to_one_refine(bool enforce_what[3])
    {
      bool did_change = false;
      ScalarIntFieldType *refine_level = m_eMesh.get_fem_meta_data()->get_field<ScalarIntFieldType>("refine_level");
      if (!refine_level)
        {
          throw std::logic_error("must have refine_level field for hanging-node refinement");
        }
      ScalarFieldType *refine_field = static_cast<ScalarFieldType *>(m_field);

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_level);
        fields.push_back(refine_field);
        //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);

      }

      if (0 && debug_print)
        {
          std::ostringstream ostr;
          stk::mesh::EntityId id[]={12,13};
          for (unsigned i=0; i < 2; i++)
            {
              stk::mesh::Entity element = m_eMesh.get_entity(m_eMesh.element_rank(), id[i]);
              if (m_eMesh.is_valid(element))
                {
                  int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                  double *refine_field_elem = m_eMesh.get_bulk_data()->field_data( *refine_field , element );
                  ostr << "P[" << m_eMesh.get_rank() << "] elem " << id[i] << " lev= " << refine_level_elem[0]
                       << " ref_field= " << refine_field_elem[0] << "\n";
                }
            }
          std::cout << ostr.str() << std::endl;
        }

      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
      const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (1 || on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (m_eMesh.numChildren(element))
                    continue;

                  int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                  double *refine_field_elem = m_eMesh.get_bulk_data()->field_data( *refine_field , element );
                  std::set<stk::mesh::Entity> selected_neighbors;
                  get_neighbors(element, refine_level, enforce_what, selected_neighbors);
                  //std::cout << "selected_neighbors.size= " << selected_neighbors.size() << std::endl;
                  for (std::set<stk::mesh::Entity>::iterator neighbor = selected_neighbors.begin();
                       neighbor != selected_neighbors.end(); ++neighbor)
                    {
                      if (m_eMesh.numChildren(*neighbor))
                        continue;

                      int *refine_level_neigh = m_eMesh.get_bulk_data()->field_data( *refine_level , *neighbor );
                      double *refine_field_neigh = m_eMesh.get_bulk_data()->field_data( *refine_field , *neighbor );

                      // if any neighbor is my level + 1 and marked for refine,
                      //   and I am not marked for refine, mark me for refine
                      if ( (refine_level_neigh[0] > refine_level_elem[0])
                           && refine_field_neigh[0] > 0.0
                           && refine_field_elem[0] <= 0.0)
                        {
                          refine_field_elem[0] = 1.0;
                          if (debug_print)
                            std::cout << "P[" << m_eMesh.get_rank() << "] enforce_two_to_one_refine: upgrading element (" 
                                      << m_eMesh.identifier(element) << "," << m_eMesh.isGhostElement(element) << ")"
                                      << " due to neighbor: " 
                                      << m_eMesh.identifier(*neighbor) << "," << m_eMesh.isGhostElement(*neighbor) << ")"
                                      << std::endl;
                          did_change = true;
                        }
                    }
                }
            }
        }

      if (0 && debug_print)
        {
          std::ostringstream ostr;
          stk::mesh::EntityId id[]={12,13};
          unsigned nid=1;
          for (unsigned i=0; i < nid; i++)
            {
              stk::mesh::Entity element = m_eMesh.get_entity(m_eMesh.element_rank(), id[i]);
              if (m_eMesh.is_valid(element))
                {
                  int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                  double *refine_field_elem = m_eMesh.get_bulk_data()->field_data( *refine_field , element );
                  ostr << "P[" << m_eMesh.get_rank() << "] t4 b4 copy_owned_to_shared elem " << id[i] << " lev= " << refine_level_elem[0]
                       << " ref_field= " << refine_field_elem[0] << "\n";
                }
            }
          std::cout << ostr.str() << std::endl;
        }

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_level);
        fields.push_back(refine_field);
        //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
        stk::mesh::communicate_field_data( m_eMesh.get_bulk_data()->shared_aura(), fields);

      }

      if (0 && debug_print)
        {
          MPI_Barrier( MPI_COMM_WORLD );
          std::ostringstream ostr;
          stk::mesh::EntityId id[]={12,13};
          unsigned nid=1;
          for (unsigned i=0; i < nid; i++)
            {
              stk::mesh::Entity element = m_eMesh.get_entity(m_eMesh.element_rank(), id[i]);
              if (m_eMesh.is_valid(element))
                {
                  int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                  double *refine_field_elem = m_eMesh.get_bulk_data()->field_data( *refine_field , element );
                  ostr << "P[" << m_eMesh.get_rank() << "] t4 after copy_owned_to_shared elem " << id[i] << " lev= " << refine_level_elem[0]
                       << " ref_field= " << refine_field_elem[0] << "\n";
                }
            }
          std::cout << ostr.str() << std::endl;
        }

      stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceMax<1>( &did_change ) );

      return did_change;
    }

    bool ElementRefinePredicate::enforce_two_to_one_unrefine(bool enforce_what[3])
    {
      bool did_change = false;
      ScalarIntFieldType *refine_level = m_eMesh.get_fem_meta_data()->get_field<ScalarIntFieldType>("refine_level");
      if (!refine_level)
        {
          throw std::logic_error("must have refine_level field for hanging-node refinement");
        }
      ScalarFieldType *refine_field = static_cast<ScalarFieldType *>(m_field);

      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_level);
        fields.push_back(refine_field);
        //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);
      }

      stk::mesh::Selector on_locally_owned_part =  ( m_eMesh.get_fem_meta_data()->locally_owned_part() );
      const std::vector<stk::mesh::Bucket*> & buckets = m_eMesh.get_bulk_data()->buckets( m_eMesh.element_rank() );
      for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
        {
          if (1 || on_locally_owned_part(**k))
            {
              stk::mesh::Bucket & bucket = **k ;
              const unsigned num_elements_in_bucket = bucket.size();
              for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
                {
                  stk::mesh::Entity element = bucket[iElement];
                  if (m_eMesh.numChildren(element))
                    continue;

                  int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                  double *refine_field_elem = m_eMesh.get_bulk_data()->field_data( *refine_field , element );
                  std::set<stk::mesh::Entity> selected_neighbors;
                  get_neighbors(element, refine_level, enforce_what, selected_neighbors);
                  //std::cout << "selected_neighbors.size= " << selected_neighbors.size() << std::endl;
                  for (std::set<stk::mesh::Entity>::iterator neighbor = selected_neighbors.begin();
                       neighbor != selected_neighbors.end(); ++neighbor)
                    {
                      if (m_eMesh.numChildren(*neighbor))
                        continue;

                      int *refine_level_neigh = m_eMesh.get_bulk_data()->field_data( *refine_level , *neighbor );
                      double *refine_field_neigh = m_eMesh.get_bulk_data()->field_data( *refine_field , *neighbor );

                      // if any neighbor is more refined (level is higher)
                      //   and I am marked for unrefine, unmark me for unrefine
                      if ( (refine_level_neigh[0] > refine_level_elem[0])
                           && refine_field_elem[0] < 0.0)
                        {
                          refine_field_elem[0] = 0.0;
                          if (debug_print)
                            std::cout << "P[" << m_eMesh.get_rank() << "] enforce_two_to_one_unrefine: downgrading element " << m_eMesh.identifier(element)
                                      << " due to neighbor: " << m_eMesh.identifier(*neighbor)
                                      << " with refine_field_neigh= " << refine_field_neigh[0]
                                      << std::endl;
                          did_change = true;
                        }
                    }
                }
            }
        }
      {
        std::vector< const stk::mesh::FieldBase *> fields;
        fields.push_back(refine_level);
        fields.push_back(refine_field);
        //stk::mesh::copy_owned_to_shared( *m_eMesh.get_bulk_data(), fields);
        stk::mesh::communicate_field_data(m_eMesh.get_bulk_data()->shared_aura(), fields);
      }

      stk::ParallelMachine pm = m_eMesh.get_bulk_data()->parallel();
      stk::all_reduce( pm, stk::ReduceMax<1>( &did_change ) );

      return did_change;
    }

    bool ElementRefinePredicate::is_face_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level)
    {

      stk::mesh::Entity parent_element = element;
      stk::mesh::Entity parent_neighbor = neighbor;
      if (std::abs(element_level - neighbor_level) > 1) return false;
      if (element_level > neighbor_level)
        {
          parent_element = m_eMesh.getParent(element, true);
        }
      else if (element_level < neighbor_level)
        {
          parent_neighbor = m_eMesh.getParent(neighbor, true);
        }
      bool fn = m_eMesh.is_face_neighbor(parent_element, parent_neighbor);
      return fn;
    }

    bool ElementRefinePredicate::is_edge_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level)
    {
      stk::mesh::Entity parent_element = element;
      stk::mesh::Entity parent_neighbor = neighbor;
      if (std::abs(element_level - neighbor_level) > 1) return false;
      if (element_level > neighbor_level)
        {
          parent_element = m_eMesh.getParent(element, true);
        }
      else if (element_level < neighbor_level)
        {
          parent_neighbor = m_eMesh.getParent(neighbor, true);
        }
      bool fn = m_eMesh.is_edge_neighbor(parent_element, parent_neighbor);
      return fn;
    }

    bool ElementRefinePredicate::is_node_neighbor(stk::mesh::Entity element, int element_level, stk::mesh::Entity neighbor, int neighbor_level)
    {
      stk::mesh::Entity parent_element = element;
      stk::mesh::Entity parent_neighbor = neighbor;
      if (std::abs(element_level - neighbor_level) > 1) return false;
      if (element_level > neighbor_level)
        {
          parent_element = m_eMesh.getParent(element, true);
        }
      else if (element_level < neighbor_level)
        {
          parent_neighbor = m_eMesh.getParent(neighbor, true);
        }
      bool fn = m_eMesh.is_node_neighbor(parent_element, parent_neighbor);
      return fn;
    }

    bool ElementRefinePredicate::min_max_neighbors_level(stk::mesh::Entity element, int min_max[2], ScalarIntFieldType *refine_level, bool check_what[3])
    {
      min_max[0] = std::numeric_limits<int>::max();
      min_max[1] = 0;
      bool found_a_face_neighbor = false;

      std::set<stk::mesh::Entity> neighbors;
      m_eMesh.get_node_neighbors(element, neighbors);
      int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
      for (std::set<stk::mesh::Entity>::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor)
        {
          int *refine_level_neigh = m_eMesh.get_bulk_data()->field_data( *refine_level , *neighbor );
          bool fn = false;
          if (check_what[2])
            fn = fn || is_face_neighbor(element, refine_level_elem[0], *neighbor, refine_level_neigh[0]);
          if (check_what[1] && m_eMesh.get_spatial_dim() == 3)
            fn = fn || is_edge_neighbor(element, refine_level_elem[0], *neighbor, refine_level_neigh[0]);
          if (check_what[0])
            fn = fn || is_node_neighbor(element, refine_level_elem[0], *neighbor, refine_level_neigh[0]);

          if (fn)
            {
              min_max[0] = std::min(min_max[0], refine_level_neigh[0]);
              min_max[1] = std::max(min_max[1], refine_level_neigh[0]);
              found_a_face_neighbor = true;
            }
        }
      return found_a_face_neighbor;
    }

    void ElementRefinePredicate::get_neighbors(stk::mesh::Entity element, ScalarIntFieldType *refine_level, bool get_what[3],
                                               std::set<stk::mesh::Entity>& selected_neighbors)
    {
      std::set<stk::mesh::Entity> neighbors;
      m_eMesh.get_node_neighbors(element, neighbors);
      //std::cout << "node_neighbors.size= " << neighbors.size() << std::endl;
      int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
      for (std::set<stk::mesh::Entity>::iterator neighbor = neighbors.begin(); neighbor != neighbors.end(); ++neighbor)
        {
          int *refine_level_neigh = m_eMesh.get_bulk_data()->field_data( *refine_level , *neighbor );
          bool fn = false;
          if (get_what[2])
            {
              bool isfn = is_face_neighbor(element, refine_level_elem[0], *neighbor, refine_level_neigh[0]);
              fn = fn || isfn;
            }
          if (get_what[1] && m_eMesh.get_spatial_dim() == 3)
            fn = fn || is_edge_neighbor(element, refine_level_elem[0], *neighbor, refine_level_neigh[0]);
          if (get_what[0])
            fn = fn || is_node_neighbor(element, refine_level_elem[0], *neighbor, refine_level_neigh[0]);

          if (fn)
            {
              selected_neighbors.insert(*neighbor);
            }
        }
    }

    bool ElementRefinePredicate::check_two_to_one(bool check_what[3])
    {
      bool valid = true;
      ScalarIntFieldType *refine_level = m_eMesh.get_fem_meta_data()->get_field<ScalarIntFieldType>("refine_level");
      if (!refine_level)
        {
          throw std::logic_error("must have refine_level field for hanging-node refinement");
        }

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
                  if (m_eMesh.numChildren(element))
                    continue;

                  int *refine_level_elem = m_eMesh.get_bulk_data()->field_data( *refine_level , element );
                  std::set<stk::mesh::Entity> selected_neighbors;
                  get_neighbors(element, refine_level, check_what, selected_neighbors);
                  //std::cout << "selected_neighbors.size= " << selected_neighbors.size() << std::endl;
                  for (std::set<stk::mesh::Entity>::iterator neighbor = selected_neighbors.begin();
                       neighbor != selected_neighbors.end(); ++neighbor)
                    {
                      if (m_eMesh.numChildren(*neighbor))
                        continue;

                      int *refine_level_neigh = m_eMesh.get_bulk_data()->field_data( *refine_level , *neighbor );

                      if ( std::abs(refine_level_neigh[0] - refine_level_elem[0]) > 1)
                        {
                          if (debug_print)
                            std::cout << "check_two_to_one: invalid element (id,level)= ("
                                      << m_eMesh.identifier(element) << ", " << refine_level_elem[0]
                                      << ") due to neighbor: (" << m_eMesh.identifier(*neighbor)
                                      << ", " << refine_level_neigh[0] << ")"
                                      << std::endl;
                          valid = false;
                        }
                    }
                }
            }
        }
      return valid;

    }


  }
}

