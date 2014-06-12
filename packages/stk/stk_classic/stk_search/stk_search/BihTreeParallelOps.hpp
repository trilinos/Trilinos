/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_BihTreeParallelOps_hpp
#define stk_search_BihTreeParallelOps_hpp

#include <stk_search/BihTree.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/OctTreeOps.hpp>

namespace stk_classic {
namespace search {

template <class DomainBoundingBox, class RangeBoundingBox>
bool oct_tree_bih_tree_proximity_search(
  ParallelMachine            arg_comm ,
  const float        * const arg_global_box ,
  const size_t               arg_domain_boxes_number ,
  const DomainBoundingBox * const arg_domain_boxes ,
  const size_t               arg_range_boxes_number ,
  const RangeBoundingBox * const arg_range_boxes ,
  const OctTreeKey   * const arg_cuts ,
  std::vector< std::pair< typename DomainBoundingBox::Key,  typename RangeBoundingBox::Key > > & arg_relation ,
  unsigned * const arg_search_tree_stats = NULL )
{
  typedef typename DomainBoundingBox::Key DomainKey;
  typedef typename RangeBoundingBox::Key RangeKey;
  typedef std::map< stk_classic::OctTreeKey, std::pair< std::list< DomainBoundingBox >, std::list< RangeBoundingBox > > > SearchTree ;

  enum { Dim = 3 };

  const unsigned p_rank = parallel_machine_rank( arg_comm );

  //----------------------------------------------------------------------
  // Search tree defined by oct-tree covering for boxes

  bool local_violations = false ;
  bool global_violations = false ;

  SearchTree search_tree ;

  {
    stk_classic::OctTreeKey covering[8] ;
    unsigned   number = 0 ;

    double scale = arg_global_box[0+Dim] - arg_global_box[0];
    for ( unsigned i = 1 ; i < Dim ; ++i ) {
      double tst_scale = arg_global_box[i+Dim] - arg_global_box[i];
      if (tst_scale > scale) scale = tst_scale;
    }
    if (scale > 0.0) // Not an error. Could arise with a point bounding box in the range/domain...
      scale = 1.0 / scale;
    else
      scale = 1.0;

    for ( size_t i = 0 ; i < arg_domain_boxes_number ; ++i ) {

      DomainBoundingBox tmp( arg_domain_boxes[i] );

      tmp.key.proc = p_rank ;

      float box[6];

      for (int x=0; x<Dim; ++x) {
        box[x] = tmp.lower(x);
        box[x+Dim] = tmp.upper(x);
      }

      const bool valid =
        hsfc_box_covering( arg_global_box, box, covering, number, scale );

      if ( ! valid ) { local_violations = true ; }

      for ( unsigned k = 0 ; k < number ; ++k ) {
        const stk_classic::OctTreeKey key = covering[k] ;
        search_tree[key].first.push_back(tmp);
      }
    }

    for ( size_t i = 0 ; i < arg_range_boxes_number ; ++i ) {

      RangeBoundingBox tmp( arg_range_boxes[i] );

      tmp.key. proc = p_rank ;

      float box[6];

      for (int x=0; x<Dim; ++x) {
        box[x] = tmp.lower(x);
        box[x+Dim] = tmp.upper(x);
      }
      
      const bool valid =
        hsfc_box_covering( arg_global_box, box, covering, number, scale );

      if ( ! valid ) { local_violations = true ; }

      for ( unsigned k = 0 ; k < number ; ++k ) {
        const stk_classic::OctTreeKey key = covering[k] ;
        search_tree[key].second.push_back(tmp);
      }
    }
  }

  //----------------------------------------------------------------------
  // Use a set to provide a unique and sorted result.

  std::set< std::pair<DomainKey, RangeKey> > tmp_relation ;

  {
    // Communicate search_tree members

    SearchTree local_tree ;

    std::set< std::pair<DomainKey, RangeKey> > local_relation ;

    if ( arg_cuts ) {
      global_violations =
        communicate<DomainBoundingBox, RangeBoundingBox>( arg_comm , arg_cuts , search_tree , local_tree ,
                           local_violations );
    }
    else {
      const double tolerance = 0.001 ;

      std::vector< stk_classic::OctTreeKey > cuts ;

      oct_tree_partition( arg_comm , search_tree , tolerance , cuts );

      global_violations =
        communicate<DomainBoundingBox, RangeBoundingBox>(arg_comm , & cuts[0] , search_tree , local_tree ,
                          local_violations );
    }

    // Local proximity search with received members

    if ( arg_search_tree_stats ) {
      search_tree_statistics( arg_comm , local_tree ,
          arg_search_tree_stats );
    }

    std::set<DomainBoundingBox, stk_classic::search::box::compare::Compare< DomainBoundingBox, box::compare::KEY> > domain;
    std::set<RangeBoundingBox,  stk_classic::search::box::compare::Compare< RangeBoundingBox , box::compare::KEY> > range;

    for (typename SearchTree::iterator i = local_tree.begin(); i != local_tree.end(); ++i) {
      //insert domain list
      std::list<DomainBoundingBox> & domain_list = i->second.first;
      for (typename std::list<DomainBoundingBox>::iterator j = domain_list.begin(); j != domain_list.end(); ++j)
        domain.insert(*j);

      //insert range list
      std::list<RangeBoundingBox> & range_list = i->second.second;
      for (typename std::list<RangeBoundingBox>::iterator j = range_list.begin(); j != range_list.end(); ++j)
        range.insert(*j);
    }

    stk_classic::search::bih::BihTree<RangeBoundingBox> tree(range.begin(),range.end());
    for (typename  std::set<DomainBoundingBox, box::compare::Compare<DomainBoundingBox,box::compare::KEY> >::iterator i = domain.begin();
        i != domain.end() ; ++i)
    {
      tree.intersect(*i,local_relation);
    }




    // Communicate relations back to domain and range processors

    communicate<DomainBoundingBox,RangeBoundingBox>( arg_comm , local_relation , tmp_relation );
  }

  arg_relation.clear();
  arg_relation.reserve( tmp_relation.size() );

  typename std::set< std::pair<DomainKey,RangeKey> >::iterator ir ;
  for ( ir = tmp_relation.begin() ; ir != tmp_relation.end() ; ++ir ) {
    arg_relation.push_back( *ir );
  }
  return global_violations ;
}

} //namespace search
} //namespace stk_classic
#endif //stk_search_BihTreeParallelOps_hpp
