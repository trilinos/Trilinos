// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// @HEADER

#ifndef THYRA_SPMD_MULTI_VECTOR_FILE_IO_HPP
#define THYRA_SPMD_MULTI_VECTOR_FILE_IO_HPP

#include "Thyra_SpmdMultiVectorSerializer.hpp"
#include "Teuchos_Utils.hpp"

namespace Thyra {

/** \brief Utility class for reading and writing parallel (or any serial)
 * Thyra vectors to and from parallel files.
 */
template<class Scalar>
class SpmdMultiVectorFileIO {
public:

  /** \brief . */
  SpmdMultiVectorFileIO(
    int   procRank  = -1
    ,int  numProcs  = -1
    );
  
  /** \brief Set the processor rank and the total number of processors.
   *
   * If <tt>numProcs < 0</tt> then procRank and numProcs will be determined
   * from <tt>Teuchos::GlobalMPISession</tt>.
   */
  void setProcRankAndSize(
    int   procRank  = -1
    ,int  numProcs  = -1
    );
  
  /** brief . */
  std::string getParallelFileName( const std::string &fileNameBase ) const;

  /** brief . */
  Teuchos::RefCountPtr<VectorBase<Scalar> >
  readVectorFromFile(
    const std::string                                            &fileNameBase
    ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &vs
    ,const bool                                                  binary = false
    ) const;

  /** brief . */
  void writeToFile(
    const MultiVectorBase<Scalar>    &mv
    ,const std::string               &fileNameBase
    ,const bool                      binary         = false
    ) const;

private:

  std::string parallelExtension_;

};

// ////////////////////////////
// Implementations!

template<class Scalar>
SpmdMultiVectorFileIO<Scalar>::SpmdMultiVectorFileIO(
  int   procRank
  ,int  numProcs
  )
{
  setProcRankAndSize(procRank,numProcs);
}

template<class Scalar>
void SpmdMultiVectorFileIO<Scalar>::setProcRankAndSize(
  int   procRank
  ,int  numProcs
  )
{
  parallelExtension_ = Teuchos::Utils::getParallelExtension(procRank,numProcs);
}

template<class Scalar>
std::string
SpmdMultiVectorFileIO<Scalar>::getParallelFileName( 
  const std::string &fileNameBase
  ) const
{
  std::ostringstream parallelFileName;
  parallelFileName << fileNameBase;
  if(parallelExtension_.length())
    parallelFileName << "." << parallelExtension_;
  return parallelFileName.str();
}

template<class Scalar>
Teuchos::RefCountPtr<VectorBase<Scalar> >
SpmdMultiVectorFileIO<Scalar>::readVectorFromFile(
  const std::string                                            &fileNameBase
  ,const Teuchos::RefCountPtr<const VectorSpaceBase<Scalar> >  &vs
  ,const bool                                                  binary
  ) const
{
  const std::string fileName = getParallelFileName(fileNameBase);
  std::ifstream in_file(fileName.c_str());
  TEST_FOR_EXCEPTION(
    in_file.eof(), std::logic_error
    ,"Error, the file \""<<fileName<<"\" could not be opened for input!"
    );
  Teuchos::RefCountPtr<VectorBase<Scalar> >
    vec = createMember(vs);
  SpmdMultiVectorSerializer<Scalar> mvSerializer(binary);
  mvSerializer.deserialize(in_file,&*vec);
  return vec;
}

template<class Scalar>
void SpmdMultiVectorFileIO<Scalar>::writeToFile(
  const MultiVectorBase<Scalar>    &mv
  ,const std::string               &fileNameBase
  ,const bool                      binary
  ) const
{
  const std::string fileName = getParallelFileName(fileNameBase);
  std::ofstream out_file(fileName.c_str());
  SpmdMultiVectorSerializer<Scalar> mvSerializer(binary);
  mvSerializer.serialize(mv,out_file);
}

} // namespace Thyra

#endif // THYRA_SPMD_MULTI_VECTOR_FILE_IO_HPP
