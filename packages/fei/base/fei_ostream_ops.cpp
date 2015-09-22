/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_ostream_ops.hpp>

#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <fei_FillableMat.hpp>
#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>


FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Vector& vec)
{
  vec.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Matrix& mat)
{
  mat.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::FillableMat& mat)
{
  os << "num rows: " << mat.getNumRows() << FEI_ENDL;
  fei::FillableMat::iterator
    iter = mat.begin(), iter_end = mat.end();

  for(; iter!=iter_end; ++iter) {
    int row = iter->first;
    const fei::CSVec* v = iter->second;
    const std::vector<int>& v_ind = v->indices();
    const std::vector<double>& v_coef = v->coefs();
    os << row << ": ";
    for(size_t i=0; i<v_ind.size(); ++i) {
      os << "("<<v_ind[i]<<","<<v_coef[i]<<") ";
    }
    os << FEI_ENDL;
  }

  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::CSVec& vec)
{
  size_t len = vec.size();

  os << "   numEntries: " << len << FEI_ENDL;

  for(size_t i=0; i<len; ++i) {
    os << "     " << vec.indices()[i]<< ": "<<vec.coefs()[i] << FEI_ENDL;
  }

  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::CSRMat& mat)
{
  os << "num rows: " << mat.getNumRows() << FEI_ENDL;

  const std::vector<int>& rows = mat.getGraph().rowNumbers;
  const int* rowoffs = &(mat.getGraph().rowOffsets[0]);
  const std::vector<int>& cols = mat.getGraph().packedColumnIndices;
  const double* coefs = &(mat.getPackedCoefs()[0]);

  for(size_t i=0; i<rows.size(); ++i) {
    int row = rows[i];

    os << row << ": ";
    for(int j=rowoffs[i]; j<rowoffs[i+1]; ++j) {
      os << "("<<cols[j]<<","<<coefs[j]<<") ";
    }
    os << FEI_ENDL;
  }

  return(os);
}

