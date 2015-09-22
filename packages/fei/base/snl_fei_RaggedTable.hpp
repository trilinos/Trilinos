#ifndef _snl_fei_RaggedTable_hpp_
#define _snl_fei_RaggedTable_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <snl_fei_SetTraits_specialize.hpp>
#include <snl_fei_MapTraits_specialize.hpp>

#include <fei_IndexTable.hpp>
#include <fei_Pool_alloc.hpp>

namespace snl_fei {

/** Data-structure that accumulates row-column indices into a ragged table,
  useful for building a matrix-graph and other concepts where keys are mapped
  to lists of values. This class can use various maps as the underlying
  data-holders. Useful because specialized maps and sets are defined in fei to
  take advantage of data that contains significant chunks of contiguous indices.
*/
template<typename MAP_TYPE, typename SET_TYPE>
class RaggedTable : public fei::IndexTable {
 public:
  /** Constructor */
  RaggedTable(int firstKey,
	      int lastKey);

  /** Copy constructor */
  RaggedTable(const RaggedTable<MAP_TYPE,SET_TYPE>& src);

  virtual ~RaggedTable();

  /** alias for MAP_TYPE */
  typedef MAP_TYPE map_type;

  /** alias for SET_TYPE */
  typedef SET_TYPE row_type;

  /** add entries to the diagonal of the table */
  void addDiagonals(int numIndices,
		    const int* indices);

  /** add a list of indices to a specified row */
  void addIndices(int row,
                  int numIndices,
                  const int* indices);

  /** add a list of indices to several specified rows */
  void addIndices(int numRows,
                  const int* rows,
                  int numIndices,
                  const int* indices);

  /** obtain internal map attribute */
  MAP_TYPE& getMap();

  /** obtain internal map attribute */
  const MAP_TYPE& getMap() const;

  /** obtain specified row from internal map attribute */
  SET_TYPE* getRow(int row);

  /** let 'iterator' be an alias for MAP_TYPE's iterator
  */
  typedef typename MAP_TYPE::iterator iterator;

  /** 'first' row of table */
  iterator begin();

  /** just past the 'last' row of the table */
  iterator end();

  /** Test for equality of two RaggedTable objects. */
  bool equal(const RaggedTable<MAP_TYPE,SET_TYPE>& rhs, bool quiet=true) const;

 private:
  MAP_TYPE map_;
  fei_Pool_alloc<SET_TYPE> poolAllocatorSet_;
  SET_TYPE dummy;
}; //class RaggedTable

template<typename MAP_TYPE, typename SET_TYPE>
inline RaggedTable<MAP_TYPE,SET_TYPE>::RaggedTable(int firstKey,
                                          int lastKey)
  : map_(),
  poolAllocatorSet_(),
  dummy()
{
}

template<typename MAP_TYPE, typename SET_TYPE>
inline RaggedTable<MAP_TYPE,SET_TYPE>::RaggedTable(const RaggedTable<MAP_TYPE,SET_TYPE>& src)
  : map_(src.map_),
    poolAllocatorSet_()
{
}

template<typename MAP_TYPE, typename SET_TYPE>
RaggedTable<MAP_TYPE,SET_TYPE>::~RaggedTable()
{
  iterator it = begin();
  iterator it_end = end();
  for(; it!=it_end; ++it) {
    poolAllocatorSet_.destroy( it->second );
    poolAllocatorSet_.deallocate( it->second, 1 );
  }
}

template<typename MAP_TYPE, typename SET_TYPE>
inline void RaggedTable<MAP_TYPE,SET_TYPE>::addIndices(int row,
                                              int numIndices,
                                              const int* indices)
{
  iterator m_end = map_.end();
  iterator m_iter = MapTraits<MAP_TYPE>::lower_bound(map_, row);

  SET_TYPE* mapped_indices = NULL;

  bool found_row = false;
  if (m_iter != m_end) {
    if ((*m_iter).first == row) {
      mapped_indices = (*m_iter).second;
      found_row = true;
    }
  }

  if (!found_row) {
    mapped_indices = poolAllocatorSet_.allocate(1);
    poolAllocatorSet_.construct(mapped_indices, dummy);
    typename MAP_TYPE::value_type val(row, mapped_indices);
    MapTraits<MAP_TYPE>::insert(map_, m_iter, val);
  }

  for(int i=0; i<numIndices; ++i) {
    SetTraits<SET_TYPE>::insert(mapped_indices, indices[i]);
  }
}

template<typename MAP_TYPE, typename SET_TYPE>
inline void RaggedTable<MAP_TYPE,SET_TYPE>::addIndices(int numRows,
                             const int* rows,
                             int numIndices,
                             const int* indices)
{
  iterator m_end = map_.end();
  iterator m_iter;
  SET_TYPE* mapped_indices = NULL;

  for(int i=0; i<numRows; ++i) {
    int row = rows[i];
    m_iter = MapTraits<MAP_TYPE>::lower_bound(map_, row);

    bool found_row = false;
    if (m_iter != m_end) {
      const typename MAP_TYPE::value_type& m_pair = *m_iter;
      if (m_pair.first == row) {
        mapped_indices = m_pair.second;
        found_row = true;
      }
    }

    if (!found_row) {
      mapped_indices = poolAllocatorSet_.allocate(1);
      poolAllocatorSet_.construct(mapped_indices, dummy);
      typename MAP_TYPE::value_type val(row, mapped_indices);
      MapTraits<MAP_TYPE>::insert(map_, m_iter, val);
    }

    for(int j=0; j<numIndices; ++j) {
      SetTraits<SET_TYPE>::insert(mapped_indices, indices[j]);
    }
  }
}

template<typename MAP_TYPE, typename SET_TYPE>
inline MAP_TYPE& RaggedTable<MAP_TYPE,SET_TYPE>::getMap()
{
  return(map_);
}

template<typename MAP_TYPE, typename SET_TYPE>
inline const MAP_TYPE& RaggedTable<MAP_TYPE,SET_TYPE>::getMap() const
{
  return(map_);
}

template<typename MAP_TYPE, typename SET_TYPE>
inline typename RaggedTable<MAP_TYPE,SET_TYPE>::row_type*
RaggedTable<MAP_TYPE,SET_TYPE>::getRow(int row)
{
  iterator m_end = map_.end();
  iterator m_iter = map_.find(row);
  return( m_end == m_iter ? NULL : (*m_iter).second );
}

template<typename MAP_TYPE, typename SET_TYPE>
inline typename RaggedTable<MAP_TYPE,SET_TYPE>::iterator
RaggedTable<MAP_TYPE,SET_TYPE>::begin()
{
  return(map_.begin());
}

template<typename MAP_TYPE, typename SET_TYPE>
inline typename RaggedTable<MAP_TYPE,SET_TYPE>::iterator
RaggedTable<MAP_TYPE,SET_TYPE>::end()
{
  return(map_.end());
}

template<typename MAP_TYPE, typename SET_TYPE>
inline void RaggedTable<MAP_TYPE,SET_TYPE>::addDiagonals(int numIndices,
                                                const int* indices)
{
  for(int i=0; i<numIndices; ++i) {
    int ind = indices[i];
    addIndices(ind, 1, &ind);
  }
}

template<typename MAP_TYPE, typename SET_TYPE>
bool RaggedTable<MAP_TYPE,SET_TYPE>::equal(const RaggedTable<MAP_TYPE,SET_TYPE>& rhs, bool quiet) const
{
  if (map_.size() != rhs.getMap().size()) {
    if (!quiet) {
      FEI_COUT << "RaggedTable::equal sizes don't match." << FEI_ENDL;
    }
    return(false);
  }

  typename map_type::const_iterator
   m_iter = map_.begin(),
   m_end  = map_.end();

  typename map_type::const_iterator
   rhs_iter = rhs.getMap().begin(),
   rhs_end  = rhs.getMap().end();

  for(; m_iter != m_end; ++m_iter, ++rhs_iter) {
    if (rhs_iter->first != m_iter->first) {
      if (!quiet) {
        FEI_COUT << "RaggedTable::equal keys don't match." << FEI_ENDL;
      }
      return(false);
    }

    if (*(rhs_iter->second) != *(m_iter->second)) {
      if (!quiet) {
        FEI_COUT << "RaggedTable::equal row-values don't match." << FEI_ENDL;
      }
      return(false);
    }
  }

  return(true);
}

}//namespace snl_fei

#endif

