/*@HEADER
// ***********************************************************************
//
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_DenseRow_HPP
#define TIFPACK_DenseRow_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"

#define HAVE_TIFPACK_ARRAY_BOUNDSCHECK

#ifdef HAVE_TIFPACK_ARRAY_BOUNDSCHECK
#include <sstream>
#include <stdexcept>
#endif

namespace Tifpack {

/** Hold a sparse matrix-row in dense memory.
 * This allows for fast lookup (as opposed to holding the sparse row
 * in a std::map or a hash-table).
 */
template<class Scalar,class GlobalOrdinal>
class DenseRow {
public:
  /** Specify the range of indices that will be used to access the
   * contents of the row. The DenseRow object allocates memory of
   * length last_index-first_index+1 to hold the row's coefficients.
   */
  DenseRow(GlobalOrdinal first_index, GlobalOrdinal last_index);

  /** Destructor */
  ~DenseRow();

  /** Access the scalar corresponding to column 'idx'.
   * (Internally returns m_values[idx - first_index].)
   */
  Scalar& operator[](GlobalOrdinal idx);
  /** Access the scalar corresponding to column 'idx'.
   * (Internally returns m_values[idx - first_index].)
   */
  const Scalar& operator[](GlobalOrdinal idx) const;

  /** Return the number of distinct entries that have been referenced
   * by 'operator[]'.
   */
  size_t getNumEntries() const;

  /** Set all internal values to 0, and flag them as unused.
   */
  void reset();

  /** Copy indices and values to arrays.
   * (Only copies the ones that have been referenced using the
   *  'operator[]' accessor.)
   */
  void copyToArrays(Teuchos::Array<GlobalOrdinal>& indices,
                    Teuchos::Array<Scalar>& values) const;

private:
  GlobalOrdinal m_first_index;
  GlobalOrdinal m_last_index;
  Teuchos::Array<Scalar> m_values;
  typedef typename Teuchos::Array<Scalar>::size_type tsize_t;
  bool* m_used;
};//class DenseRow


template<class Scalar,class GlobalOrdinal>
DenseRow<Scalar,GlobalOrdinal>::DenseRow(GlobalOrdinal first_index, GlobalOrdinal last_index)
 : m_first_index(first_index),
   m_last_index(last_index),
   m_values(last_index-first_index+1, 0),
   m_used(NULL)
{
  m_used = new bool[m_values.size()];
  for(tsize_t i=0; i<m_values.size(); ++i) {
    m_used[i] = false;
  }
}

template<class Scalar,class GlobalOrdinal>
DenseRow<Scalar,GlobalOrdinal>::~DenseRow()
{
  delete [] m_used; m_used = NULL;
}

template<class Scalar,class GlobalOrdinal>
inline Scalar&
DenseRow<Scalar,GlobalOrdinal>::operator[](GlobalOrdinal idx)
{
#ifdef HAVE_TIFPACK_ARRAY_BOUNDSCHECK
  if (idx < m_first_index || idx > m_last_index) {
    std::ostringstream os;
    os << "Tifpack::DenseRow ERROR, input idx("<<idx<<") is outside "
       << "the bounds m_first_index("<<m_first_index<<") .. m_last_index("
       << m_last_index << ").";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
#endif

  GlobalOrdinal offset = idx-m_first_index;
  m_used[offset] = true;
  return m_values[offset];
}

template<class Scalar,class GlobalOrdinal>
inline const Scalar&
DenseRow<Scalar,GlobalOrdinal>::operator[](GlobalOrdinal idx) const
{
#ifdef HAVE_TIFPACK_ARRAY_BOUNDSCHECK
  if (idx < m_first_index || idx > m_last_index) {
    std::ostringstream os;
    os << "Tifpack::DenseRow ERROR, input idx("<<idx<<") is outside "
       << "the bounds m_first_index("<<m_first_index<<") .. m_last_index("
       << m_last_index << ").";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
#endif

  GlobalOrdinal offset = idx-m_first_index;
  m_used[offset] = true;
  return m_values[offset];
}

template<class Scalar,class GlobalOrdinal>
size_t DenseRow<Scalar,GlobalOrdinal>::getNumEntries() const
{
  size_t num_entries = 0;
  for(tsize_t i=0; i<m_values.size(); ++i) {
    if (m_used[i]) ++num_entries;
  }
  return num_entries;
}

template<class Scalar,class GlobalOrdinal>
void DenseRow<Scalar,GlobalOrdinal>::reset()
{
  for(tsize_t i=0; i<m_values.size(); ++i) {
    m_values[i] = 0;
    m_used[i] = false;
  }
}

template<class Scalar,class GlobalOrdinal>
void DenseRow<Scalar,GlobalOrdinal>::copyToArrays(Teuchos::Array<GlobalOrdinal>& indices,
                    Teuchos::Array<Scalar>& values) const
{
  size_t num_entries = getNumEntries();
  indices.resize(num_entries);
  values.resize(num_entries);
  tsize_t offset = 0;
  for(tsize_t i=0; i<m_values.size(); ++i) {
    if (m_used[i]) {
      values[offset] = m_values[i];
      indices[offset++] = m_first_index+i;
    }
  }
}

}//namespace Tifpack

#endif // TIFPACK_DenseRow_HPP

