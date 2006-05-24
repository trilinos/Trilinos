#include "Thyra_MultiVectorSerialization.hpp"
#include "Teuchos_Utils.hpp"

namespace Thyra {

/** \brief Utility class for reading and writing parallel (or any serial)
 * Thyra vectors to and from parallel files.
 */
template<class Scalar>
class ParallelMultiVectorFileIO {
public:

  /** \brief . */
  ParallelMultiVectorFileIO(
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
ParallelMultiVectorFileIO<Scalar>::ParallelMultiVectorFileIO(
  int   procRank
  ,int  numProcs
  )
{
  setProcRankAndSize(procRank,numProcs);
}

template<class Scalar>
void ParallelMultiVectorFileIO<Scalar>::setProcRankAndSize(
  int   procRank
  ,int  numProcs
  )
{
  parallelExtension_ = Teuchos::Utils::getParallelExtension(procRank,numProcs);
}

template<class Scalar>
std::string
ParallelMultiVectorFileIO<Scalar>::getParallelFileName( 
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
ParallelMultiVectorFileIO<Scalar>::readVectorFromFile(
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
  MultiVectorSerialization<Scalar> mvSerializer(binary);
  mvSerializer.unserialize(in_file,&*vec);
  return vec;
}

template<class Scalar>
void ParallelMultiVectorFileIO<Scalar>::writeToFile(
  const MultiVectorBase<Scalar>    &mv
  ,const std::string               &fileNameBase
  ,const bool                      binary
  ) const
{
  const std::string fileName = getParallelFileName(fileNameBase);
  std::ofstream out_file(fileName.c_str());
  MultiVectorSerialization<Scalar> mvSerializer(binary);
  mvSerializer.serialize(mv,out_file);
}

} // namespace Thyra
