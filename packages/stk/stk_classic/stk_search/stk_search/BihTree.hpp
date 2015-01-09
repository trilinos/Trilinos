/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_search_BihTree_hpp
#define stk_search_BihTree_hpp

#include <float.h>

#include <algorithm>
#include <list>
#include <vector>
#include <set>

#include <stk_search/BoundingBox.hpp>

namespace stk_classic {
namespace search {
namespace bih {

class UnitTestBihTree;

/** Bih tree */

template <class RangeBoundingVolume, int COMPARE_MODE = stk_classic::search::box::compare::MIDDLE>
class BihTree
{
  public:
    typedef typename RangeBoundingVolume::Key RangeKey;


    template <class InputIterator>
    BihTree(
        const InputIterator & first,
        const InputIterator & last,
        unsigned max_depth = 32,
        unsigned items_per_leaf = 1 ) :
      MAX_DEPTH(max_depth), ITEMS_PER_LEAF(items_per_leaf), m_range(), m_root(), m_depth(0), m_num_nodes(0)
    {
      m_range.insert(m_range.begin(),first,last);
      build_tree(m_root, m_range.begin(), m_range.end(), 0);
    }

    /**
     * RelationshipContainer can be
     *  std::vector< RangeBoundingVoulme >
     *  std::vector< std::pair< DomainBoundingVolume::Key, RangeBoundingVolume::Key> > >
     *  std::set<    std::pair< DomainBoundingVolume::Key, RangeBoundingVolume::Key> > >
     */
    template<class DomainBoundingVolume, class RelationshipContainer>
    void intersect(
        const DomainBoundingVolume & domain,
        RelationshipContainer  & relations) const
    {
      intersect_helper(m_root,domain,relations);
    }

    /**
     * RelationshipContainer can be
     *  std::vector< RangeBoundingVoulme >
     *  std::vector< std::pair< DomainBoundingVolume::Key, RangeBoundingVolume::Key> > >
     *  std::set<    std::pair< DomainBoundingVolume::Key, RangeBoundingVolume::Key> > >
     */
    template<class DomainIterator, class RelationshipContainer>
    void intersect(
        const DomainIterator & first,
        const DomainIterator & last,
        RelationshipContainer  & relations) const
    {
      for (DomainIterator i = first; i != last; ++i)
        intersect_helper(m_root,*i,relations);
    }

    /**
     * RelationshipContainer can be
     *  std::vector< RangeBoundingVoulme >
     *  std::vector< std::pair< RangeBoundingVolume::Key, RangeBoundingVolume::Key> > >
     *  std::set<    std::pair< RangeBoundingVolume::Key, RangeBoundingVolume::Key> > >
     */
    template<class RelationshipContainer>
    void self_intersect (RelationshipContainer & relations) const {
      intersect(m_range.begin(),m_range.end(),relations);
    }

    unsigned depth() const { return m_depth; }
    unsigned num_nodes() const { return m_num_nodes; }

    const unsigned MAX_DEPTH;
    const unsigned ITEMS_PER_LEAF;

  private:

    /** Node for a BihTree.  Each Node is 16 bytes.  Leaf Nodes store an array of type K.
    */
    class Node
    {
      private: // functions to use inside of BihTree
        static const uint64_t LEAF_NODE = 3;
        static const uint64_t BIT_MASK = 3;

        Node()
          : m_children_ptr(LEAF_NODE)
        {
          m_data.volumes = 0;
        }

        ~Node() {
          if (!is_leaf()) // remove volumes
            delete [] children();
        }

        /** return what axis the node is split along */
        inline int axis() const { return static_cast<int>(m_children_ptr & BIT_MASK); }

        /** true if the node is a leaf node */
        inline bool is_leaf() const { return (m_children_ptr & BIT_MASK) == LEAF_NODE; }

        /** number of items in the leaf node
         *  returns -1 if the node is not a leaf node  */
        inline int size() const { return is_leaf() ? static_cast<int>(m_children_ptr >> 2) : -1; }

        /** begin of key list
         * it is an unchecked error to call begin or end on an interior node */
        inline const RangeBoundingVolume * begin() const { return m_data.volumes; }

        /** end of key list
         * it is an unchecked error to call begin or end on an interior node */
        inline const RangeBoundingVolume * end() const { return m_data.volumes + size(); }

        /** if !is_leaf() returns a pointer to a pair of children, otherwise returns NULL */
        inline Node * children() const { return is_leaf() ? 0 : reinterpret_cast<Node *>(m_children_ptr & ~BIT_MASK); }

        /** left child */
        inline Node & left() const { return children()[0]; }

        /** right child */
        inline Node & right() const { return children()[1]; }

        /** left_max clip plane
         * it is an unchecked error to call left_max or right_min on a leaf node */
        inline float left_max() const { return m_data.planes[0]; }

        /** right_min clip plane
         * it is an unchecked error to call left_max or right_min on a leaf node */
        inline float right_min() const { return m_data.planes[1]; }

        /** make the current node and interior node */
        inline void make_interior_node(uint64_t axis, const float clipplanes[]) {
          Node * child = new Node[2];
          m_children_ptr = reinterpret_cast<uint64_t>(child) | axis;
          m_data.planes[0] = clipplanes[0];
          m_data.planes[1] = clipplanes[1];
        }

        /** make the current node a leaf node */
        void make_leaf_node(const RangeBoundingVolume * volumes, size_t num) {
          // the size is stored in the upper bytes of the children_ptr
          m_children_ptr = (num << 2) | LEAF_NODE;
          m_data.volumes = volumes;
        }


        /** stores either a pointer to the children or the number of volumes in the leaf in the upper bits
         * the lowest 2 bits store the axis */
        uint64_t m_children_ptr;

        union {
          const RangeBoundingVolume * volumes  ;
          float                  planes[2] ;
        }m_data;

        //disallow copy constructor and assignment operator
        Node(const Node &);
        void operator = (const Node &);

        friend class BihTree<RangeBoundingVolume,COMPARE_MODE>;
        friend class ::stk_classic::search::bih::UnitTestBihTree;
    }; //end class Node


  private:

    /** build a tree using the middle_longest_side heristic */
    void build_tree(
        Node & node,
        typename std::vector<RangeBoundingVolume>::iterator begin,
        typename std::vector<RangeBoundingVolume>::iterator end,
        unsigned current_depth);

    template<class DomainBoundingVolume>
      void intersect_helper(
          const Node & node,
          const DomainBoundingVolume & domain,
          std::vector<std::pair<typename DomainBoundingVolume::Key, RangeKey> >  & relations) const;

    template<class DomainBoundingVolume>
      void intersect_helper(
          const Node & node,
          const DomainBoundingVolume & domain,
          std::set<std::pair<typename DomainBoundingVolume::Key, RangeKey> >  & relations) const;

    template<class DomainBoundingVolume>
      void intersect_helper(
          const Node & node,
          const DomainBoundingVolume & domain,
          std::vector<RangeBoundingVolume>  & range) const;
  private:


    /** find max value of the bounding boxes along the given axis */
    float find_max(
        typename std::vector<RangeBoundingVolume>::iterator begin,
        typename std::vector<RangeBoundingVolume>::iterator end,
        int axis) const
    {
      float value = begin->upper(axis);

      for (typename std::vector<RangeBoundingVolume>::const_iterator i = begin+1; i!=end; ++i) {
        if ( static_cast<float>(i->upper(axis)) > value) {
          value = static_cast<float>(i->upper(axis));
        }
      }
      return value;
    }

    /** find ratio of left to right child using the given split plane and axis */
    float find_ratio(
        typename std::vector<RangeBoundingVolume>::iterator begin,
        typename std::vector<RangeBoundingVolume>::iterator end,
        int axis,
        float split_plane) const
    {
      int count = 0;
      for (typename std::vector<RangeBoundingVolume>::const_iterator i = begin; i!=end; ++i) {
        if ( static_cast<float>(i->center(axis)) < split_plane) {
          ++count;
        }
      }
      float value = static_cast<float>(count)/(end-begin) - 0.5f;
      return value >=0 ? value : -1 * value;
    }

    /** find min value of the bounding boxes along the given axis */
    float find_min(
        typename std::vector<RangeBoundingVolume>::iterator begin,
        typename std::vector<RangeBoundingVolume>::iterator end,
        int axis) const
    {
      float value = begin->lower(axis);

      for (typename std::vector<RangeBoundingVolume>::const_iterator i = begin+1; i!=end; ++i) {
        if ( static_cast<float>(i->lower(axis)) < value) {
          value = static_cast<float>(i->lower(axis));
        }
      }
      return value;
    }


    bool find_split_plane(
        const typename std::vector<RangeBoundingVolume>::iterator begin,
        const typename std::vector<RangeBoundingVolume>::iterator end,
        int & axis,
        typename RangeBoundingVolume::Data & split_plane )
    {
      bool split_by_space = true;
      axis = 0;

      int size = end - begin;

      typename RangeBoundingVolume::Data length = 0, mean_length = 0, max_length = 0;


      for (int j = 0; j < RangeBoundingVolume::DIMENSION; ++j) {
        typename RangeBoundingVolume::Data min = FLT_MAX, max = -FLT_MAX, temp_mean_length = 0, temp_max_length = 0;
        for (typename std::vector<RangeBoundingVolume>::const_iterator i = begin; i!=end; ++i) {
          temp_mean_length += i->length(j);
          if (i->length(j) > temp_max_length) temp_max_length = i->length(j);
          if (i->lower(j) < min) min = i->lower(j);
          if (i->upper(j) > max) max = i->upper(j);
        }
        temp_mean_length /= size;
        if (max-min > length) {
          length = max-min;
          axis = j;
          split_plane = min + length/2;
          mean_length = temp_mean_length;
          max_length = temp_max_length;
        }
      }

      if (mean_length*8 < max_length) {
        split_by_space = false;
        split_plane = mean_length*5;
      }

      return split_by_space;
    }

    void make_leaf(
        Node & node,
        typename std::vector<RangeBoundingVolume>::iterator begin,
        typename std::vector<RangeBoundingVolume>::iterator end )
    {
      int num = end - begin;
      if (num > 0) {
        node.make_leaf_node(&(*begin),num);
      }
      else {
        node.make_leaf_node(0,0);
      }
    }

  private:

    //member variables
    std::vector<RangeBoundingVolume> m_range;
    Node     m_root;
    unsigned m_depth;
    unsigned m_num_nodes;


    //disable copy constructor and assignment operator
    explicit BihTree( const BihTree &);
    void operator = ( const BihTree &);

    friend class ::stk_classic::search::bih::UnitTestBihTree;
};


template <class RangeBoundingVolume, int COMPARE_MODE>
void BihTree<RangeBoundingVolume, COMPARE_MODE>::build_tree(
    Node & node,
    typename std::vector<RangeBoundingVolume>::iterator begin,
    typename std::vector<RangeBoundingVolume>::iterator end,
    unsigned current_depth)
{
  //update tree depth
  m_depth = current_depth > m_depth ? current_depth : m_depth;

  // make leaf node
  if (current_depth >= MAX_DEPTH || end-begin <= static_cast<int>(ITEMS_PER_LEAF)) {
    make_leaf(node,begin,end);
    return;
  }

  typename std::vector<RangeBoundingVolume>::iterator k = end;
  typename RangeBoundingVolume::Data split_plane = 0;
  int axis = 0;

  bool split_by_space = find_split_plane(begin,end, axis, split_plane);

  if (split_by_space) { // split by space
    k =  std::partition(begin,end,stk_classic::search::box::compare::Partition<RangeBoundingVolume,COMPARE_MODE>(axis,split_plane));

    float ratio = static_cast<float>(k-begin) / (end-begin);
    ratio -= 0.5f;
    ratio = ratio >= 0 ? ratio : -1.0f * ratio;

    if (ratio > 0.35f) {
      k = begin + (end - begin)/2;
      std::nth_element(begin,k,end,stk_classic::search::box::compare::Compare<RangeBoundingVolume,COMPARE_MODE>(axis));
    }
  }
  else { // split by length
    k =  std::partition(begin,end,stk_classic::search::box::compare::Partition<RangeBoundingVolume,stk_classic::search::box::compare::LENGTH>(axis,split_plane));
  }

  float clip_planes[2] = {-FLT_MAX,FLT_MAX};

  if (k - begin >= 1) {// left child non-empty
    clip_planes[0] = find_max(begin,k,axis);
  }

  if (end - k >= 1) { // right child non-empty
    clip_planes[1] = find_min(k,end,axis);
  }

  //make the node an interior node
  node.make_interior_node(axis,clip_planes);

  m_num_nodes += 2;

  build_tree(
      node.left(),
      begin,
      k,
      current_depth+1
      );

  build_tree(
      node.right(),
      k,
      end,
      current_depth+1
      );
}

template <class RangeBoundingVolume, int COMPARE_MODE>
template<class DomainBoundingVolume>
void BihTree<RangeBoundingVolume, COMPARE_MODE>::intersect_helper(
    const Node & node,
    const DomainBoundingVolume & domain,
    std::vector<std::pair<typename DomainBoundingVolume::Key, RangeKey> >  & relations) const
{
  //if leaf add items to vector
  if (node.is_leaf()) {
    if (node.size() > 0 ) {
      for( const RangeBoundingVolume * i = node.begin(); i != node.end(); ++i) {
	relations.push_back( std::pair<typename DomainBoundingVolume::Key, RangeKey>(domain.key,i->key) );
      }
    }
    return;
  }

  int axis = node.axis();

  if (domain.lower(axis) <= node.left_max()) { //go left
    intersect_helper( node.left(), domain, relations );
  }

  if (domain.upper(axis) >= node.right_min()) { //go right
    intersect_helper( node.right(), domain, relations );
  }

}

template <class RangeBoundingVolume, int COMPARE_MODE>
template<class DomainBoundingVolume>
void BihTree<RangeBoundingVolume, COMPARE_MODE>::intersect_helper(
    const Node & node,
    const DomainBoundingVolume & domain,
    std::set<std::pair<typename DomainBoundingVolume::Key, RangeKey> >  & relations) const
{
  //if leaf add items to vector
  if (node.is_leaf()) {
    if (node.size() > 0 ) {
      for( const RangeBoundingVolume * i = node.begin(); i != node.end(); ++i) {
	relations.insert( std::pair<typename DomainBoundingVolume::Key, RangeKey>(domain.key,i->key) );
      }
    }
    return;
  }

  int axis = node.axis();

  if (domain.lower(axis) <= node.left_max()) { //go left
    intersect_helper( node.left(), domain, relations );
  }

  if (domain.upper(axis) >= node.right_min()) { //go right
    intersect_helper( node.right(), domain, relations );
  }

}


template <class RangeBoundingVolume, int COMPARE_MODE>
template<class DomainBoundingVolume>
void BihTree<RangeBoundingVolume, COMPARE_MODE>::intersect_helper(
    const Node & node,
    const DomainBoundingVolume & domain,
    std::vector<RangeBoundingVolume> & range) const
{
  //if leaf add items to vector
  if (node.is_leaf()) {
    if (node.size() > 0 ) {
      for( const RangeBoundingVolume * i = node.begin(); i != node.end(); ++i) {
	range.push_back(*i);
      }
    }
    return;
  }

  int axis = node.axis();

  if (domain.lower(axis) <= node.left_max()) { //go left
    intersect_helper( node.left(), domain, range );
  }

  if (domain.upper(axis) >= node.right_min()) { //go right
    intersect_helper( node.right(), domain, range );
  }

}

} // namespace bih
} // namespace search
} // namespace stk_classic

#endif // stk_search_BihTree_hpp
