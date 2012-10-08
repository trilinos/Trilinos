/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_CommUtils.hpp>
#include <fei_iostream.hpp>
#include <fei_impl_utils.hpp>
#include <fei_FillableMat.hpp>
#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>
#include <fei_Graph.hpp>
#include <fei_Matrix.hpp>
#include <fei_Reducer.hpp>

namespace fei {
namespace impl_utils {

//----------------------------------------------------------------------------
void find_offsets(const std::vector<int>& sources,
                  const std::vector<int>& targets,
                  std::vector<int>& offsets)
{
  offsets.assign(sources.size(), -1);

  size_t ssize = sources.size(), tsize = targets.size();
  size_t si = 0, ti = 0;

  const int* sourcesptr = &sources[0];
  const int* targetsptr = &targets[0];

  while(si<ssize && ti<tsize) {
    if (sourcesptr[si] == targetsptr[ti]) {
      offsets[si++] = ti++;
      continue;
    }

    while(sourcesptr[si] < targetsptr[ti] && si<ssize) {
      ++si;
    }

    while(sourcesptr[si] > targetsptr[ti] && ti<tsize) {
      ++ti;
    }
  }
}

//----------------------------------------------------------------------------
size_t num_bytes_FillableMat(const fei::FillableMat& mat)
{
  int nrows = mat.getNumRows();
  int nnz = fei::count_nnz(mat);

  int num_chars_int = (2 + nrows*2 + nnz)*sizeof(int);
  int num_chars_double = nnz*sizeof(double);

  return num_chars_int + num_chars_double;
}

//----------------------------------------------------------------------------
void pack_FillableMat(const fei::FillableMat& mat, 
                      char* buffer)
{
  int nrows = mat.getNumRows();
  int nnz = fei::count_nnz(mat);

  int num_chars_int = (2 + nrows*2 + nnz)*sizeof(int);

  int* intdata = reinterpret_cast<int*>(buffer);
  double* doubledata = reinterpret_cast<double*>(buffer+num_chars_int);

  int ioffset = 0;
  int doffset = 0;

  intdata[ioffset++] = nrows;
  intdata[ioffset++] = nnz;

  int ioffsetcols = 2+nrows*2;

  fei::FillableMat::const_iterator
    r_iter = mat.begin(),
    r_end = mat.end();

  for(; r_iter!=r_end; ++r_iter) {
    int rowNumber = r_iter->first;
    const fei::CSVec* row = r_iter->second;

    intdata[ioffset++] = rowNumber;
    const int rowlen = row->size();
    intdata[ioffset++] = rowlen;

    const std::vector<int>& rowindices = row->indices();
    const std::vector<double>& rowcoefs = row->coefs();
    for(int i=0; i<rowlen; ++i) {
      intdata[ioffsetcols++] = rowindices[i];
      doubledata[doffset++] = rowcoefs[i];
    }
  }
}

//----------------------------------------------------------------------------
void unpack_FillableMat(const char* buffer_begin, const char* buffer_end,
                        fei::FillableMat& mat,
                        bool clear_mat_on_entry,
                        bool overwrite_entries)
{
  if (clear_mat_on_entry) {
    mat.clear();
  }

  if (buffer_end == buffer_begin) {
    return;
  }

  const int* intdata = reinterpret_cast<const int*>(buffer_begin);
  int ioffset = 0;
  int nrows = intdata[ioffset++];
  int nnz = intdata[ioffset++];

  int ioffsetcols = 2+nrows*2;

  int num_chars_int = (2+nrows*2 + nnz)*sizeof(int);
  const double* doubledata = reinterpret_cast<const double*>(buffer_begin+num_chars_int);

  int doffset = 0;

  for(int i=0; i<nrows; ++i) {
    int row = intdata[ioffset++];
    int rowlen = intdata[ioffset++];

    for(int j=0; j<rowlen; ++j) {
      int col = intdata[ioffsetcols++];
      double coef = doubledata[doffset++];

      if (overwrite_entries) {
        mat.putCoef(row, col, coef);
      }
      else {
        mat.sumInCoef(row, col, coef);
      }
    }
  }

  if (doffset != nnz) {
    throw std::runtime_error("fei::impl_utils::unpack_FillableMat: failed, sizes don't agree.");
  }
}

//----------------------------------------------------------------------------
bool unpack_CSRMat(const char* buffer_begin, const char* buffer_end, fei::CSRMat& mat)
{
  bool all_zeros = true;
  if (buffer_end == buffer_begin) {
    return all_zeros;
  }

  const int* intdata = reinterpret_cast<const int*>(buffer_begin);
  int ioffset = 0;
  int nrows = intdata[ioffset++];
  int nnz = intdata[ioffset++];

  fei::SparseRowGraph& srg = mat.getGraph();
  srg.rowNumbers.resize(nrows);
  srg.rowOffsets.resize(nrows+1);
  srg.packedColumnIndices.resize(nnz);
  std::vector<double>& packed_coefs = mat.getPackedCoefs();
  packed_coefs.resize(nnz);

  int ioffsetcols = 2+nrows*2;

  int num_chars_int = (2+nrows*2 + nnz)*sizeof(int);
  const double* doubledata = reinterpret_cast<const double*>(buffer_begin+num_chars_int);

  int doffset = 0;

  for(int i=0; i<nrows; ++i) {
    int row = intdata[ioffset++];
    int rowlen = intdata[ioffset++];

    srg.rowNumbers[i] = row;
    srg.rowOffsets[i] = doffset;

    for(int j=0; j<rowlen; ++j) {
      int col = intdata[ioffsetcols++];
      double coef = doubledata[doffset];
      if (coef != 0.0) all_zeros = false;
      srg.packedColumnIndices[doffset] = col;
      packed_coefs[doffset++] = coef;
    }
  }
  srg.rowOffsets[nrows] = nnz;
  return all_zeros;
}

size_t num_bytes_indices_coefs(const std::vector<int>& indices,
                        const std::vector<double>& coefs)
{
  int num = indices.size();
  int num_chars_int = (1+num)*sizeof(int);
  int num_chars = num_chars_int + num*sizeof(double);
  return num_chars;
}

void pack_indices_coefs(const std::vector<int>& indices,
                        const std::vector<double>& coefs,
                        std::vector<char>& buffer,
                        bool resize_buffer)
{
  if (indices.size() != coefs.size()) {
    throw std::runtime_error("fei::impl_utils::pack_indices_coefs failed, sizes don't match.");
  }

  int num = indices.size();
  int num_chars_int = (1+num)*sizeof(int);
  int num_chars = num_chars_int + num*sizeof(double);
  if (resize_buffer) {
    buffer.resize(num_chars);
  }

  int* intdata = reinterpret_cast<int*>(&buffer[0]);
  double* doubledata = reinterpret_cast<double*>(&buffer[0]+num_chars_int);

  int ioffset = 0;
  int doffset = 0;
  intdata[ioffset++] = num;
  for(int i=0; i<num; ++i) {
    intdata[ioffset++] = indices[i];
    doubledata[doffset++] = coefs[i];
  }
}

void unpack_indices_coefs(const std::vector<char>& buffer,
                          std::vector<int>& indices,
                          std::vector<double>& coefs)
{
  if (buffer.size() == 0) return;

  const int* intdata = reinterpret_cast<const int*>(&buffer[0]);
  int ioffset = 0;
  int num = intdata[ioffset++];
  int num_chars_int = (1+num)*sizeof(int);
  const double* doubledata = reinterpret_cast<const double*>(&buffer[0]+num_chars_int);

  indices.resize(num);
  coefs.resize(num);

  int doffset = 0;
  for(int i=0; i<num; ++i) {
    indices[i] = intdata[ioffset++];
    coefs[i] = doubledata[doffset++];
  }
}

//----------------------------------------------------------------------------
void separate_BC_eqns(const fei::FillableMat& mat,
                    std::vector<int>& bcEqns,
                    std::vector<double>& bcVals)
{
  fei::FillableMat::const_iterator
    m_iter = mat.begin(),
    m_end = mat.end();

  for(; m_iter != m_end; ++m_iter) {
    int eqn = m_iter->first;
    const fei::CSVec* row = m_iter->second;

    std::vector<int>::iterator
      iter = std::lower_bound(bcEqns.begin(), bcEqns.end(), eqn);
    if (iter == bcEqns.end() || *iter != eqn) {
      size_t offset = iter - bcEqns.begin();
      bcEqns.insert(iter, eqn);
      std::vector<double>::iterator viter = bcVals.begin();
      viter += offset;
      try {
        double val = get_entry(*row, eqn);
        bcVals.insert(viter, val);
      }
      catch(...) {
        fei::console_out() << "separate_BC_eqns ERROR, failed to find coef for eqn " << eqn;
        throw;
      }
    }
  }
}

//----------------------------------------------------------------------------
void create_col_to_row_map(const fei::FillableMat& mat,
                           std::multimap<int,int>& crmap)
{
  crmap.clear();

  if (mat.getNumRows() == 0) return;

  std::vector<int> rowNumbers;
  fei::get_row_numbers(mat, rowNumbers);

  fei::FillableMat::const_iterator
    m_iter = mat.begin(),
    m_end = mat.end();

  for(; m_iter != m_end; ++m_iter) {
    int rowNum = m_iter->first;
    const fei::CSVec* rowvec = m_iter->second;

    const std::vector<int>& rowindices = rowvec->indices();

    for(size_t i=0; i<rowindices.size(); ++i) {
      int colNum = rowindices[i];

      crmap.insert(std::make_pair(colNum, rowNum));
    }
  }
}

//----------------------------------------------------------------------------
int remove_couplings(fei::FillableMat& mat)
{
  int levelsOfCoupling = 0;

  std::vector<int> rowNumbers;
  fei::get_row_numbers(mat, rowNumbers);

  bool finished = false;
  while(!finished) {
    std::multimap<int,int> crmap;
    create_col_to_row_map(mat, crmap);

    typedef std::multimap<int,int>::iterator MM_Iter;

    fei::FillableMat::iterator
      m_iter = mat.begin(),
      m_end = mat.end();

    bool foundCoupling = false;
    for(; m_iter != m_end; ++m_iter) {
      int rownum = m_iter->first;
      fei::CSVec* mrow = m_iter->second;

      //now find which rows contain 'rownum' as a column-index:
      std::pair<MM_Iter,MM_Iter> mmi = crmap.equal_range(rownum);

      //loop over the rows which contain 'rownum' as a column-index:
      for(MM_Iter cri = mmi.first; cri != mmi.second; ++cri) {
        int cri_row = cri->second;

        fei::CSVec* frow = mat.create_or_getRow(cri_row);

        double coef = get_entry(*frow,rownum);

        remove_entry(*frow, rownum);

        std::vector<int>& indices = mrow->indices();
        std::vector<double>& coefs = mrow->coefs();

        size_t rowlen = mrow->size();
        for(size_t ii=0; ii<rowlen; ++ii) {
          coefs[ii] *= coef;
        }

        add_entries(*frow, rowlen, &indices[0], &coefs[0]);
        foundCoupling = true;
      }
    }

    if (foundCoupling == true) ++levelsOfCoupling;
    else finished = true;
  }

  if (levelsOfCoupling > 1) {
    fei::console_out() << "fei::removeCouplings WARNING, levelsOfCoupling="
      << levelsOfCoupling << " (Each level of coupling means that a slave DOF "
      << "depends on a master which is itself a slave.) Behavior is not well "
      << "understood for levelsOfCoupling greater than 1..."<<FEI_ENDL;
  }

  return levelsOfCoupling;
}

//----------------------------------------------------------------------------
void global_union(MPI_Comm comm,
                  const fei::FillableMat& localMatrix,
                  fei::FillableMat& globalUnionMatrix)
{
  globalUnionMatrix = localMatrix;

  int localProc = fei::localProc(comm);
  int numProcs = fei::numProcs(comm);

  if (numProcs < 2) {
    return;
  }

  //first pack the local matrix into a pair of std::vector objects

  size_t num_bytes = num_bytes_FillableMat(localMatrix);
  std::vector<char> localchardata(num_bytes);

  pack_FillableMat(localMatrix, &localchardata[0]);

  //next use Allgatherv to place every processor's packed arrays onto every
  //other processor.

  std::vector<int> recvdatalengths;
  std::vector<char> recvdata;
  int err = fei::Allgatherv(comm, localchardata, recvdatalengths, recvdata);
  if (err != 0) {
    throw std::runtime_error("fei::impl_utils::global_union: Allgatherv-int failed.");
  }

  //finally unpack the received arrays into matrix objects and combine them
  //into the globalUnionMatrix object.
  bool clearMatOnEntry = false;
  bool overwriteEntries = true;

  int ioffset = 0;
  for(size_t p=0; p<recvdatalengths.size(); ++p) {
    int len = recvdatalengths[p];

    if (len > 1 && (int)p != localProc) {
      unpack_FillableMat(&recvdata[ioffset], &recvdata[ioffset]+len,
                globalUnionMatrix, clearMatOnEntry, overwriteEntries);
    }

    ioffset += len;
  }
}

//----------------------------------------------------------------------------
void global_union(MPI_Comm comm,
                  const fei::CSVec& localVec, fei::CSVec& globalUnionVec)
{
  globalUnionVec.indices().clear();
  globalUnionVec.coefs().clear();

  const std::vector<int>& localindices = localVec.indices();
  const std::vector<double>& localcoefs = localVec.coefs();

  std::vector<int> localintdata;
  if (localindices.size() > 0) {
    localintdata.assign(&localindices[0], &localindices[0]+localindices.size());
  }

  std::vector<double> localdoubledata;
  if (localcoefs.size() > 0) {
    localdoubledata.assign(&localcoefs[0], &localcoefs[0]+localcoefs.size());
  }

  //use Allgatherv to place every processor's arrays onto every
  //other processor.

  std::vector<int> recvintdatalengths;
  std::vector<int> recvintdata;
  int err = fei::Allgatherv(comm, localintdata, recvintdatalengths, recvintdata);
  if (err != 0) {
    throw std::runtime_error("snl_fei::globalUnion csvec: Allgatherv-int failed.");
  }

  std::vector<int> recvdoubledatalengths;
  std::vector<double> recvdoubledata;
  err = fei::Allgatherv(comm, localdoubledata, recvdoubledatalengths, recvdoubledata);
  if (err != 0) {
    throw std::runtime_error("snl_fei::globalUnion csvec: Allgatherv-double failed.");
  }

  if (recvintdatalengths.size() != recvdoubledatalengths.size()) {
    throw std::runtime_error("snl_fei::globalUnion csvec: inconsistent lengths from Allgatherv");
  }

  //now unpack the received arrays into the globalUnionVec object.

  unsigned len = recvintdata.size();
  if (len > 0) {
    int* recvintPtr = &recvintdata[0];
    double* recvdoublePtr = &recvdoubledata[0];

    for(unsigned i=0; i<len; ++i) {
      fei::put_entry(globalUnionVec, recvintPtr[i], recvdoublePtr[i]);
    }
  }
}

//----------------------------------------------------------------------------
void translate_to_reduced_eqns(const fei::Reducer& reducer, fei::CSRMat& mat)
{
  fei::SparseRowGraph& srg = mat.getGraph();

  std::vector<int>& rowNumbers = srg.rowNumbers;
  for(size_t i=0; i<rowNumbers.size(); ++i) {
    rowNumbers[i] = reducer.translateToReducedEqn(rowNumbers[i]);
  }

  std::vector<int>& colIndices = srg.packedColumnIndices;
  for(size_t i=0; i<colIndices.size(); ++i) {
    colIndices[i] = reducer.translateToReducedEqn(colIndices[i]);
  }
}

//----------------------------------------------------------------------------
void translate_to_reduced_eqns(const fei::Reducer& reducer, fei::CSVec& vec)
{
  std::vector<int>& indices = vec.indices();
  for(size_t i=0; i<indices.size(); ++i) {
    indices[i] = reducer.translateToReducedEqn(indices[i]);
  }
}

//----------------------------------------------------------------------------
void add_to_graph(const fei::CSRMat& inmat, fei::Graph& graph)
{
  const std::vector<int>& rowNumbers = inmat.getGraph().rowNumbers;
  const std::vector<int>& rowOffsets = inmat.getGraph().rowOffsets;
  const std::vector<int>& pckColInds = inmat.getGraph().packedColumnIndices;

  for(size_t i=0; i<rowNumbers.size(); ++i) {
    int row = rowNumbers[i];
    int offset = rowOffsets[i];
    int rowlen = rowOffsets[i+1]-offset;
    const int* indices = &pckColInds[offset];

    if (graph.addIndices(row, rowlen, indices) != 0) {
      throw std::runtime_error("fei::impl_utils::add_to_graph ERROR in graph.addIndices.");
    }
  }
}

//----------------------------------------------------------------------------
void add_to_matrix(const fei::CSRMat& inmat, bool sum_into, fei::Matrix& matrix)
{
  const std::vector<int>& rowNumbers = inmat.getGraph().rowNumbers;
  const std::vector<int>& rowOffsets = inmat.getGraph().rowOffsets;
  const std::vector<int>& pckColInds = inmat.getGraph().packedColumnIndices;
  const std::vector<double>& pckCoefs = inmat.getPackedCoefs();

  for(size_t i=0; i<rowNumbers.size(); ++i) {
    int row = rowNumbers[i];
    int offset = rowOffsets[i];
    int rowlen = rowOffsets[i+1]-offset;
    const int* indices = &pckColInds[offset];
    const double* coefs = &pckCoefs[offset];

    if (sum_into) {
      if (matrix.sumIn(1, &row, rowlen, indices, &coefs, FEI_DENSE_ROW) != 0) {
        throw std::runtime_error("fei::impl_utils::add_to_matrix ERROR in matrix.sumIn.");
      }
    }
    else {
      if (matrix.copyIn(1, &row, rowlen, indices, &coefs, FEI_DENSE_ROW) != 0) {
        throw std::runtime_error("fei::impl_utils::add_to_matrix ERROR in matrix.copyIn.");
      }
    }
  }
}

}//namespace impl_utils
}//namespace fei

