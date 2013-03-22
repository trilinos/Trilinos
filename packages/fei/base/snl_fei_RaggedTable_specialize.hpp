#ifndef _snl_fei_RaggedTable_specialize_hpp_
#define _snl_fei_RaggedTable_specialize_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <snl_fei_RaggedTable.hpp>
#include <fei_ctg_set.hpp>
#include <snl_fei_MapContig.hpp>

namespace snl_fei {

/** specialization for MapContig<fei::ctg_set<int> > */
template<>
class RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >
  : public fei::IndexTable {
public:
  RaggedTable(int firstKey, int lastKey);

  RaggedTable(const RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >& src);

  typedef MapContig<fei::ctg_set<int>*> map_type;
  typedef fei::ctg_set<int> row_type;
  typedef map_type::iterator iterator;

  virtual ~RaggedTable()
  {
    iterator it = begin();
    iterator it_end = end();
    for(; it!=it_end; ++it) {
      poolAllocatorSet_.destroy( (*it).second );
      poolAllocatorSet_.deallocate( (*it).second, 1 );
    }
  }

  void addDiagonals(int numIndices,
                    const int* indices);

  void addIndices(int row,
                  int numIndices,
                  const int* indices);

  void addIndices(int numRows,
                  const int* rows,
                  int numIndices,
                  const int* indices);

  map_type& getMap() { return( map_ ); }

  row_type* getRow(int row)
  {
    iterator m_end = map_.end();
    iterator m_iter = map_.find(row);
    return( m_end == m_iter ? NULL : (*m_iter).second );
  }

  iterator begin() { return( map_.begin() ); }

  iterator end() { return( map_.end() ); }

 private:
  map_type map_;
  fei_Pool_alloc<row_type> poolAllocatorSet_;
  row_type dummy;
};//RaggedTable<MapContig<fei::ctg_set<int>*> >

}//namespace snl_fei
#endif

