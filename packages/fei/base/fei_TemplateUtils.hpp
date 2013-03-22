#ifndef _fei_TemplateUtils_hpp_
#define _fei_TemplateUtils_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_iosfwd.hpp>
#include <fei_SparseRowGraph.hpp>
#include <snl_fei_RaggedTable.hpp>
#include <snl_fei_Utils.hpp>

#include <fei_ErrMacros.hpp>

#include <set>
#include <vector>
#include <map>

namespace fei {

  //------------------------------------------------------------------------
  /** dangerous function to copy a set to an
   array, assuming the set contents are of type int
  */
  template<typename SET_TYPE>
    void copySetToArray(const SET_TYPE& set_obj,
			int lenList,
			int* list)
    {
      int setsize = set_obj.size();
      int len = setsize > lenList ? lenList : setsize;

      typename SET_TYPE::const_iterator
        s_iter = set_obj.begin(),
               s_end = set_obj.end();
      for(int i=0; i<len; ++i, ++s_iter) {
        list[i] = *s_iter;
      }
    }

  /** copy a set to a vector
  */
  template<typename T>
    void copySetToVector(const std::set<T>& set_obj, std::vector<T>& vec)
    {
      vec.clear();
      vec.reserve(set_obj.size());

      typename std::set<T>::const_iterator
        s_iter = set_obj.begin(), s_end = set_obj.end();
      for(; s_iter != s_end; ++s_iter) {
        vec.push_back(*s_iter);
      }
    }

  /** dangerous function to copy map keys to an
   array, assuming the keys are of type int
  */
  template<typename MAP_TYPE>
    void copyKeysToArray(const MAP_TYPE& map_obj,
			 unsigned lenList,
			 int* list)
    {
      unsigned i = 0;
      typename MAP_TYPE::const_iterator
        iter = map_obj.begin(), iter_end = map_obj.end();

      for(; iter != iter_end; ++iter) {
        if (i == lenList) break;
        list[i++] = iter->first;
      }
    }

  /** function to copy map keys to a
   vector, assuming the keys are of type int
  */
  template<typename MAP_TYPE>
    void copyKeysToVector(const MAP_TYPE& map_obj,
                          std::vector<int>& keyvector)
    {
      keyvector.resize(map_obj.size());
      if (map_obj.empty()) return;
      copyKeysToArray<MAP_TYPE>(map_obj, map_obj.size(), &keyvector[0]);
    }

  template<typename T, typename U>
    void copyMapOfSetsToVectorOfVectors(const std::map<T,std::set<U> >& mapset,
                                     std::vector<T>& keys,
                                     std::vector<std::vector<U> >& values)
    {
      typedef std::map<T,std::set<U> > mapsettype;
      keys.resize(mapset.size());
      values.resize(mapset.size());
      typename mapsettype::const_iterator
       ms_iter = mapset.begin(), ms_end = mapset.end();
      for(size_t i=0; ms_iter!=ms_end; ++ms_iter, ++i) {
        keys[i] = ms_iter->first;
        typename std::set<U>::const_iterator
          s_iter = ms_iter->second.begin(), s_end = ms_iter->second.end();
        values[i].resize(ms_iter->second.size());
        for(size_t j=0; s_iter!=s_end; ++s_iter, ++j) {
          values[i][j] = *s_iter;
        }
      }
    }

  /** dangerous function to copy a map object to a
   pair of arrays, assuming the keys and values of
   the map are of type int.
  */
  template<typename MAP_TYPE>
    void copyToArrays(MAP_TYPE& map_obj,
		      int lenList,
		      int* keylist,
		      int* vallist)
    {
      int i = 0;
      typename MAP_TYPE::iterator
        iter = map_obj.begin(),
             iter_end = map_obj.end();

      for(; iter != iter_end; ++iter) {
        if (i == lenList) break;
        keylist[i] = (*iter).first;
        vallist[i++] = (*iter).second;
      }
    }

  /** iterate a map object, destroying its contents */
  template<typename MAP_TYPE>
    void destroyValues(MAP_TYPE& map_obj)
    {
      typename MAP_TYPE::iterator
        m_iter = map_obj.begin(),
      m_end  = map_obj.end();

      for(; m_iter != m_end; ++m_iter) {
        delete (*m_iter).second;
      }
    }

  /** write a ragged-table to a specified ostream */
  template<typename MAP_TYPE, typename SET_TYPE>
  void writeToStream(snl_fei::RaggedTable<MAP_TYPE,SET_TYPE>& table,
		     FEI_OSTREAM& os,
		     const char* lineprefix=NULL)
  {
    MAP_TYPE& map_obj = table.getMap();
    typename MAP_TYPE::iterator
      m_iter = map_obj.begin(),
      m_end = map_obj.end();

    for(; m_iter != m_end; ++m_iter) {
      if (lineprefix != NULL) {
        os << lineprefix;
      }

      os << " row "<<(*m_iter).first<<": ";

      typename SET_TYPE::const_iterator
        s_iter = (*m_iter).second->begin(),
               s_end = (*m_iter).second->end();

      for(; s_iter != s_end; ++s_iter) {
        os << *s_iter << " ";
      }

      os << FEI_ENDL;
    }
  }

  template<typename MAP_TYPE, typename SET_TYPE>
    void packRaggedTable(snl_fei::RaggedTable<MAP_TYPE,SET_TYPE>& table,
                       std::vector<int>& intdata)
    {
      MAP_TYPE& map_obj = table.getMap();
      int numRows = map_obj.size();

      typename MAP_TYPE::iterator
        m_iter = map_obj.begin(),
        m_end  = map_obj.end();

      int nnz = 0;

      for(; m_iter != m_end; ++m_iter) {
        typename MAP_TYPE::value_type m_pair = *m_iter;

        int rowLen = m_pair.second->size();
        nnz += rowLen;
      }

      intdata.resize(1+2*numRows+nnz);
      intdata[0] = numRows;
      int* rowNumbers = &intdata[1];
      int* rowLengths = rowNumbers+numRows;
      int* packedCols = rowLengths+numRows;

      m_iter = map_obj.begin();
      unsigned offset = 0;
      for(unsigned i=0; m_iter != m_end; ++m_iter, ++i) {
        typename MAP_TYPE::value_type m_pair = *m_iter;
        rowNumbers[i] = m_pair.first;
        rowLengths[i] = m_pair.second->size();

        int* colInds = &packedCols[offset];
        copySetToArray(*(m_pair.second), rowLengths[i], colInds);
        offset += rowLengths[i];
      }
    }

  /** create fei::SparseRowGraph object from a vector of ragged-tables. user is
    responsible for destroying the fei::SparseRowGraph */
  template<typename MAP_TYPE, typename SET_TYPE>
  fei::SharedPtr<fei::SparseRowGraph> createSparseRowGraph(const std::vector<snl_fei::RaggedTable<MAP_TYPE,SET_TYPE>* >& tables)
    {
      int numRows = 0;
      int nnz = 0;
      fei::SharedPtr<fei::SparseRowGraph> srg(new fei::SparseRowGraph);

      for(unsigned i=0; i<tables.size(); ++i) {
        MAP_TYPE& map_obj = tables[i]->getMap();
        numRows += map_obj.size();

        typename MAP_TYPE::iterator
          m_iter = map_obj.begin(),
          m_end  = map_obj.end();
        for(; m_iter != m_end; ++m_iter) {
          typename MAP_TYPE::value_type m_pair = *m_iter;
          nnz += m_pair.second->size();
        }
      }

      srg->rowNumbers.resize(numRows);
      srg->rowOffsets.resize(numRows+1);
      srg->packedColumnIndices.resize(nnz);

      unsigned offset1 = 0;
      unsigned rowOffset = 0;
      for(unsigned i=0; i<tables.size(); ++i) {
        MAP_TYPE& map_obj = tables[i]->getMap();

        typename MAP_TYPE::iterator
          m_iter = map_obj.begin(),
          m_end  = map_obj.end();
        for(; m_iter != m_end; ++m_iter) {
          typename MAP_TYPE::value_type m_pair = *m_iter;
          srg->rowNumbers[offset1] = m_pair.first;
          int rowLen = m_pair.second->size();
          srg->rowOffsets[offset1++] = rowOffset;
          int* cols = &srg->packedColumnIndices[rowOffset];
          copySetToArray(*(m_pair.second), rowLen, cols);
          rowOffset += rowLen;
        }
      }

      srg->rowOffsets[offset1] = rowOffset;

      return(srg);
    }

  /** copy a ragged-table to an existing fei::SparseRowGraph object.  */
  template<typename MAP_TYPE, typename SET_TYPE>
  void copyToSparseRowGraph(snl_fei::RaggedTable<MAP_TYPE,SET_TYPE>& table,
                            fei::SparseRowGraph& srg)
    {
      MAP_TYPE& map_obj = table.getMap();
      int numRows = map_obj.size();

      srg.rowNumbers.resize(numRows);
      srg.rowOffsets.resize(numRows+1);

      int* rowNumPtr = &(srg.rowNumbers[0]);
      int* rowOffsPtr = &(srg.rowOffsets[0]);

      typename MAP_TYPE::iterator
        m_iter = map_obj.begin(),
        m_end  = map_obj.end();

      int offset = 0;
      int nnz = 0;

      for(; m_iter != m_end; ++m_iter) {
        typename MAP_TYPE::value_type m_pair = *m_iter;

        rowNumPtr[offset] = m_pair.first;
        rowOffsPtr[offset++] = nnz;
        int rowLen = m_pair.second->size();
        nnz += rowLen;
      }
      rowOffsPtr[offset] = nnz;

      srg.packedColumnIndices.resize(nnz);
      int* colPtr = &(srg.packedColumnIndices[0]);
      offset = 0;
      m_iter = map_obj.begin();
      int i = 0;
      for(; m_iter != m_end; ++m_iter, ++i) {
        typename MAP_TYPE::value_type m_pair = *m_iter;

        int rowLen = rowOffsPtr[i+1]-rowOffsPtr[i];
        int* colInds = &(colPtr[offset]);
        copySetToArray(*(m_pair.second), rowLen, colInds);
        offset += rowLen;
      }
    }

  /** create fei::SparseRowGraph object from a ragged-table.
    */
  template<typename MAP_TYPE, typename SET_TYPE>
  fei::SharedPtr<fei::SparseRowGraph>
    createSparseRowGraph(snl_fei::RaggedTable<MAP_TYPE,SET_TYPE>& table)
    {
      fei::SharedPtr<fei::SparseRowGraph> srg(new fei::SparseRowGraph);

      copyToSparseRowGraph<MAP_TYPE, SET_TYPE>(table, *srg);

      return( srg );
    }

  /** function to count the "nonzeros" in a ragged-table */
  template<typename MAP_TYPE, typename SET_TYPE>
  int countNonzeros(snl_fei::RaggedTable<MAP_TYPE,SET_TYPE>& table)
    {
      int nnz = 0;
      MAP_TYPE& map_obj = table.getMap();
      typename MAP_TYPE::iterator
        m_iter = map_obj.begin(),
               m_end = map_obj.end();

      for(; m_iter != m_end; ++m_iter) {
        nnz += (*m_iter).second->size();
      }

      return(nnz);
    }

} //namespace fei

#endif // _fei_TemplateUtils_hpp_

