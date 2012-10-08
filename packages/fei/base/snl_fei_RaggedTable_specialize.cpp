/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <snl_fei_RaggedTable_specialize.hpp>

namespace snl_fei {

/** specialization for MapContig<fei::ctg_set<int>*> */
RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::RaggedTable(int firstKey, int lastKey)
  : map_(firstKey, lastKey),
    poolAllocatorSet_(),
    dummy()
{
  int len = lastKey-firstKey+1;
  if (len > 0) {
    map_type::value_type val;
    for(int i=0; i<len; ++i) {
      val.first = firstKey+i;
      row_type* row = poolAllocatorSet_.allocate(1);
      poolAllocatorSet_.construct(row,dummy);
      val.second = row;
      map_.insert(val);
    }
  }
}

RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::RaggedTable(const RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >& src)
 : map_(src.map_),
   poolAllocatorSet_()
{
}

void RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::addIndices(int row,
                             int numIndices,
                             const int* indices)
{
  iterator m_end = map_.end();
  iterator m_iter = map_.lower_bound(row);

  map_type::mapped_type mapped_indices = (*m_iter).second;

  if (mapped_indices == NULL) {
    throw std::runtime_error("RaggedTable<MapContig>, NULL row.");
  }

  for(int i=0; i<numIndices; ++i) {
    mapped_indices->insert2(indices[i]);
  }
}

void
RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::addIndices(int numRows,
                             const int* rows,
                             int numIndices,
                             const int* indices)
{
  iterator m_end = map_.end();
  map_type::mapped_type mapped_indices = NULL;

  for(int i=0; i<numRows; ++i) {
    int row = rows[i];
    iterator m_iter = map_.lower_bound(row);

    mapped_indices = (*m_iter).second;
    if (mapped_indices == NULL) {
      throw std::runtime_error("RaggedTable<MapContig>, NULL row.");
    }

    for(int j=0; j<numIndices; ++j) {
      mapped_indices->insert2(indices[j]);
    }
  }
}

void RaggedTable<MapContig<fei::ctg_set<int>*>,fei::ctg_set<int> >::addDiagonals(int numIndices,
                                                         const int* indices)
{
  for(int i=0; i<numIndices; ++i) {
    int ind = indices[i];
    addIndices(ind, 1, &ind);
  }
}

}//namespace snl_fei

