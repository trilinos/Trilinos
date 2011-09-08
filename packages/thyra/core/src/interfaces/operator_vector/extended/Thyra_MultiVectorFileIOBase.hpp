// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_MULTI_VECTOR_FILE_IO_BASE_HPP
#define THYRA_MULTI_VECTOR_FILE_IO_BASE_HPP

#include "Thyra_OperatorVectorTypes.hpp"
#include "Teuchos_ParameterListAcceptor.hpp"
#include "Teuchos_VerboseObject.hpp"


namespace Thyra {


/** \brief Abstract strategy interface for reading and writing
 * (multi)vector objects to and from files.
 *
 * The concept of a file is every general and really can be implemented as any
 * type of object data base that is keyed on a string name
 * (i.e. <tt>fileNameBase</tt>).  In that sense, this interface is really an
 * interface to a general multi-vector serialization/deserialization
 * implementation, but file-based implementations are expected to be the most
 * common.
 *
 * This interface currently requires the client to know the correct vector
 * space and to pre-create the multi-vectors with the right number of columns
 * before they can be read in.  In all current use cases where this interface
 * is used, the client knows what it needs to read so this is fine.
 *
 * ToDo: Add a form of readMultiFromFile(...) that will accept just a vector
 * space and will create a multi-vector with as many columns as is specified
 * in the file.  Right now I don't know this functionality so I am not going
 * to implement this.  However, if an important use case is found where this
 * functionality is needed, then we can add this implementation without much
 * trouble.
 *
 * ToDo: Derive this interface from Teuchos::ParameterListAcceptor so that we
 * can set some options on how the reading and writing gets done (e.g. use
 * binary or ASCII formating).
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class MultiVectorFileIOBase
  : virtual public Teuchos::VerboseObject<MultiVectorFileIOBase<Scalar> >
{
public:

  /** \brief Return if the given multi-vector is compatible with this
   * implementation.*/
  virtual bool isCompatible( const MultiVectorBase<Scalar> &mv ) const = 0;

  /** \brief Read a (multi)vector from a file given the file base name.
   *
   * \param fileNameBase [in] The base name of the file(s) that will be used
   * to read the multi-vector from.
   *
   * \param mv [in/out] On output, this multi-vector will be filled with the
   * values from the given file(s).  This multi-vector must have already been
   * created and structured in such a way that is compatible with the format
   * of the multi-vector stored in the given file and the implementation of
   * this interface.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>mv!=NULL</tt>
   * <li><tt>this->isCompatible(*mv)==true</tt>.
   * </ul>
   */
  virtual void readMultiVectorFromFile(
    const std::string &fileNameBase,
    Thyra::MultiVectorBase<Scalar> *mv
    ) const = 0;
  
  /** \brief Write a (multi)vector to a file given the file base name.
   *
   * \param mv [in] The multi-vector that will be written to file(s).
   *
   * \param fileNameBase [in] The base name of the file(s) that will written
   * to with the values of the multi-vector.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->isCompatible(mv)==true</tt>.
   * </ul>
   */
  virtual void writeMultiVectorToFile(
    const Thyra::MultiVectorBase<Scalar> &mv,
    const std::string &fileNameBase
    ) const = 0;
  
  //@}
  
};


/** \brief Read a vector from file(s) given the file base name and a vector
 * space.
 *
 * \relates MultiVectorFileIOBase
 */
template<class Scalar>
Teuchos::RCP<VectorBase<Scalar> >
readVectorFromFile(
  const MultiVectorFileIOBase<Scalar> &fileIO,
  const std::string &fileNameBase,
  const VectorSpaceBase<Scalar> &vecSpc
  )
{
  Teuchos::RCP<VectorBase<Scalar> > v = createMember(vecSpc);
  fileIO.readMultiVectorFromFile(fileNameBase,&*v);
  return v;
}


} // namespace Thyra


#endif // THYRA_MULTI_VECTOR_FILE_IO_BASE_HPP
