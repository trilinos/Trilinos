// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_SPMD_MULTI_VECTOR_FILE_IO_HPP
#define THYRA_DEFAULT_SPMD_MULTI_VECTOR_FILE_IO_HPP

#include "Thyra_MultiVectorFileIOBase.hpp"
#include "Thyra_SpmdMultiVectorSerializer.hpp"
#include "Teuchos_Utils.hpp"

namespace Thyra {

/** \brief Concrete implementation of <tt>MultiVectorFileIO</tt> that reads
 * and writes SPMD-based (multi)vectors to and from files.
 *
 * The file accessed by each process process is
 * <tt>fileNameBase.extentionTagName.numProcs.procRank</tt>.  If
 * <tt>extentionTagName==""</tt> then the file with the name
 * <tt>fileNameBase.numProcs.procRank</tt> is accessed in each process.  By
 * setting up different file name extension information (see
 * <tt>setFileNameExtension()</tt>), the client can carefully control how file
 * base names are mapped into actual sets of files.
 *
 * ToDo: This implementation will have to be refactored once I refactor how
 * SPMD-based vectors and multi-vectors can be accessed in a general way.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class Scalar>
class DefaultSpmdMultiVectorFileIO
  : public MultiVectorFileIOBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors */
  //@{

  /** \brief Construct with file extension information (calls
   * <tt>setFileNameExtension()</tt>).
   */
  DefaultSpmdMultiVectorFileIO(
    const std::string      &extensionTagName  = ""
    ,const  int            numProcs           = -1
    ,const  int            procRank           = -1
    );
  
  /** \brief Set file name extension information to disambiguate files on
   * different processes and from other files.
   *
   * \param  extensionTagName
   *           [in] An extension name string that will be appended the beginning
   *           of full file extension.  Default is "".
   * \param  numProcs
   *           [in] The total number of processes in the communicator.
   *           Default value is <tt>-1</tt>.
   * \param  procRank
   *           [in] The rank of this process.
   *           Default value is <tt>-1</tt>.
   *
   * If <tt>numProcs < 0</tt> then <tt>procRank</tt> and <tt>numProcs</tt>
   * will be determined from <tt>Teuchos::GlobalMPISession</tt>.
   */
  void setFileNameExtension(
    const std::string      &extensionTagName  = ""
    ,const  int            numProcs           = -1
    ,const  int            procRank           = -1
    );
  
  /** brief Return the file name that is used in this process. */
  std::string getLocalFileName( const std::string &fileNameBase ) const;

  //@}

  /** \name Overridden from MultiVectorFileIOBase */
  //@{
  /** \brief . */
  bool isCompatible( const MultiVectorBase<Scalar> &mv ) const;
  /** \brief . */
  void readMultiVectorFromFile(
    const std::string                 &fileNameBase
    ,Thyra::MultiVectorBase<Scalar>   *mv
    ) const;
  /** \brief . */
  void writeMultiVectorToFile(
    const Thyra::MultiVectorBase<Scalar>   &mv
    ,const std::string                     &fileNameBase
    ) const;
  //@}

private:

  std::string    localFileNameExtension_;
  bool           useBinaryMode_;

  mutable SpmdMultiVectorSerializer<Scalar>   mvSerializer_;
  
};

// ///////////////////////////
// Implementations

template<class Scalar>
DefaultSpmdMultiVectorFileIO<Scalar>::DefaultSpmdMultiVectorFileIO(
  const std::string      &extensionTagName
  ,const  int            numProcs
  ,const  int            procRank
  )
  :useBinaryMode_(false) // ToDo: Make this adjustable!
  ,mvSerializer_(useBinaryMode_)
{
  setFileNameExtension(extensionTagName,numProcs,procRank);
}

template<class Scalar>
void DefaultSpmdMultiVectorFileIO<Scalar>::setFileNameExtension(
  const std::string      &extensionTagName
  ,const  int            numProcs
  ,const  int            procRank
  )
{
  const std::string
    endExtension = Teuchos::Utils::getParallelExtension(procRank,numProcs);
  if(extensionTagName.length())
    localFileNameExtension_ = extensionTagName+"."+endExtension;
  else
    localFileNameExtension_ = endExtension;
}

template<class Scalar>
std::string
DefaultSpmdMultiVectorFileIO<Scalar>::getLocalFileName( 
  const std::string &fileNameBase
  ) const
{
  std::ostringstream parallelFileName;
  parallelFileName << fileNameBase;
  if(localFileNameExtension_.length())
    parallelFileName << "." << localFileNameExtension_;
  return parallelFileName.str();
}

template<class Scalar>
bool DefaultSpmdMultiVectorFileIO<Scalar>::isCompatible(
  const MultiVectorBase<Scalar> &mv
  ) const
{
  return mvSerializer_.isCompatible(mv);
}

template<class Scalar>
void DefaultSpmdMultiVectorFileIO<Scalar>::readMultiVectorFromFile(
  const std::string                 &fileNameBase
  ,Thyra::MultiVectorBase<Scalar>   *mv
  ) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!mv);
  const std::string fileName = getLocalFileName(fileNameBase);
  std::ifstream in_file(fileName.c_str());
  TEUCHOS_TEST_FOR_EXCEPTION(
    in_file.eof(), std::logic_error
    ,"Error, the file \""<<fileName<<"\" could not be opened for input!"
    );
  mvSerializer_.binaryMode(useBinaryMode_);
  mvSerializer_.deserialize(in_file,mv);
}

template<class Scalar>
void DefaultSpmdMultiVectorFileIO<Scalar>::writeMultiVectorToFile(
  const Thyra::MultiVectorBase<Scalar>   &mv
  ,const std::string                     &fileNameBase
  ) const
{
  const std::string fileName = getLocalFileName(fileNameBase);
  std::ofstream out_file(fileName.c_str());
  mvSerializer_.binaryMode(useBinaryMode_);
  mvSerializer_.serialize(mv,out_file);
}

} // namespace Thyra

#endif // THYRA_DEFAULT_SPMD_MULTI_VECTOR_FILE_IO_HPP
