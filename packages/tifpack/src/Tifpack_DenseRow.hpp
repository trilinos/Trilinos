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

#ifdef HAVE_TIFPACK_DEBUG
#define HAVE_TIFPACK_ARRAY_BOUNDSCHECK
#endif

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
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;

  /** Specify the max index that will be used to access the
   * contents of the row. The DenseRow object allocates memory of
   * length max_index+1 to hold the row's coefficients.
   */
  DenseRow(GlobalOrdinal last_index);

  /** Destructor */
  ~DenseRow();

  /** Access the scalar corresponding to column 'idx'.
   */
  Scalar& operator[](GlobalOrdinal idx);
  /** Access the scalar corresponding to column 'idx'.
   */
  const Scalar& operator[](GlobalOrdinal idx) const;

  /** Return the number of non-zero entries.
   */
  size_t getNumEntries() const;

  /** Set all internal values to 0.
   */
  void reset();

  /** reduce indices and values to arrays.
   * (Discards values to satisfy fill (maxlen) and drop-tolerance)
   */
  void reduceToArraysL(
           Teuchos::Array<magnitudeType>& magRow,
           int maxlen,
           const magnitudeType& dropTol,
           GlobalOrdinal diag,
           Teuchos::Array<GlobalOrdinal>& indices,
           Teuchos::Array<Scalar>& values) const;

  /** reduce indices and values to arrays.
   * (Discards values to satisfy fill (maxlen) and drop-tolerance)
   */
  Scalar reduceToArraysU(
           Teuchos::Array<magnitudeType>& magRow,
           int maxlen,
           const magnitudeType& dropTol,
           GlobalOrdinal diag,
           Teuchos::Array<GlobalOrdinal>& indices,
           Teuchos::Array<Scalar>& values) const;

private:
  GlobalOrdinal m_last_index;
  Teuchos::Array<Scalar> m_values;
  typedef typename Teuchos::Array<Scalar>::size_type tsize_t;
  Scalar m_zero;
};//class DenseRow


template<class Scalar,class GlobalOrdinal>
DenseRow<Scalar,GlobalOrdinal>::DenseRow(GlobalOrdinal last_index)
 : m_last_index(last_index),
   m_values(last_index+1, 0),
   m_zero(Teuchos::ScalarTraits<Scalar>::zero())
{
}

template<class Scalar,class GlobalOrdinal>
DenseRow<Scalar,GlobalOrdinal>::~DenseRow()
{
}

template<class Scalar,class GlobalOrdinal>
inline Scalar&
DenseRow<Scalar,GlobalOrdinal>::operator[](GlobalOrdinal idx)
{
#ifdef HAVE_TIFPACK_ARRAY_BOUNDSCHECK
  if (idx > m_last_index) {
    std::ostringstream os;
    os << "Tifpack::DenseRow ERROR, input idx("<<idx<<") is outside "
       << "the bounds m_last_index(" << m_last_index << ").";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
#endif

  return m_values[idx];
}

template<class Scalar,class GlobalOrdinal>
inline const Scalar&
DenseRow<Scalar,GlobalOrdinal>::operator[](GlobalOrdinal idx) const
{
#ifdef HAVE_TIFPACK_ARRAY_BOUNDSCHECK
  if (idx > m_last_index) {
    std::ostringstream os;
    os << "Tifpack::DenseRow ERROR, input idx("<<idx<<") is outside "
       << "the bounds m_last_index(" << m_last_index << ").";
    std::string str = os.str();
    throw std::runtime_error(str);
  }
#endif

  return m_values[idx];
}

template<class Scalar,class GlobalOrdinal>
size_t DenseRow<Scalar,GlobalOrdinal>::getNumEntries() const
{
  size_t num_entries = 0;
  tsize_t len = m_values.size();
  for(tsize_t i=0; i<len; ++i) {
    if (m_values[i] != m_zero) ++num_entries;
  }
  return num_entries;
}

template<class Scalar,class GlobalOrdinal>
void DenseRow<Scalar,GlobalOrdinal>::reset()
{
  tsize_t len = m_values.size();
  for(tsize_t i=0; i<len; ++i) {
    m_values[i] = m_zero;
  }
}

template<class Scalar,class GlobalOrdinal>
void DenseRow<Scalar,GlobalOrdinal>::reduceToArraysL(
           Teuchos::Array<magnitudeType>& magRow,
           int maxlen,
           const magnitudeType& dropTol,
           GlobalOrdinal diag,
           Teuchos::Array<GlobalOrdinal>& indices,
           Teuchos::Array<Scalar>& values) const
{
  magRow.resize(0);
  tsize_t len = diag;
  for(tsize_t i=0; i<len; ++i) {
    magnitudeType mag = Teuchos::ScalarTraits<Scalar>::magnitude(m_values[i]);
    if (mag > dropTol) magRow.push_back(mag);
  }

  magnitudeType cutoff = dropTol;

  if (magRow.size() > maxlen) {
    std::nth_element(magRow.begin(), magRow.begin()+maxlen, magRow.begin()+magRow.size(), std::greater<magnitudeType>());
    cutoff = magRow[maxlen];
  }


  indices.resize(0);
  values.resize(0);
  for(tsize_t i=0; i<len; ++i) {
    magnitudeType mag = Teuchos::ScalarTraits<Scalar>::magnitude(m_values[i]);
    if (mag >= cutoff) {
      values.push_back(m_values[i]);
      indices.push_back(i);
    }
    //else: should discarded values be added to the diagonal?
    //in Ifpack's ILUT, discarded values are added to the diagonal.
  }
}

template<class Scalar,class GlobalOrdinal>
Scalar DenseRow<Scalar,GlobalOrdinal>::reduceToArraysU(
           Teuchos::Array<magnitudeType>& magRow,
           int maxlen,
           const magnitudeType& dropTol,
           GlobalOrdinal diag,
           Teuchos::Array<GlobalOrdinal>& indices,
           Teuchos::Array<Scalar>& values) const
{
  magRow.resize(0);
  tsize_t len = m_values.size();
  for(tsize_t i=diag; i<len; ++i) {
    magnitudeType mag = Teuchos::ScalarTraits<Scalar>::magnitude(m_values[i]);
    if (mag > dropTol) magRow.push_back(mag);
  }

  magnitudeType cutoff = dropTol;

  if (magRow.size() > maxlen) {
    std::nth_element(magRow.begin(), magRow.begin()+maxlen, magRow.begin()+magRow.size(), std::greater<magnitudeType>());
    cutoff = magRow[maxlen];
  }


  indices.resize(0);
  values.resize(0);
  values.push_back(m_values[diag]);
  indices.push_back(diag);
  for(tsize_t i=diag+1; i<len; ++i) {
    magnitudeType mag = Teuchos::ScalarTraits<Scalar>::magnitude(m_values[i]);
    if (mag >= cutoff) {
      values.push_back(m_values[i]);
      indices.push_back(i);
    }
    //else: should discarded values be added to the diagonal?
    //in Ifpack's ILUT, discarded values are added to the diagonal.
  }

  return m_values[diag];
}

}//namespace Tifpack

#endif // TIFPACK_DenseRow_HPP

