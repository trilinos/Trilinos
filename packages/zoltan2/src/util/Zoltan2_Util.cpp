// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_Util.cpp
 *  \brief Useful namespace methods.
 */

#include <Zoltan2_Util.hpp>
#include <Teuchos_OpaqueWrapper.hpp>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#include <unistd.h>

namespace Zoltan2{

/* On a linux node, find the total memory currently allocated
 * to this process.
 * Return the number of kilobytes allocated to this process.
 * Return 0 if it is not possible to determine this.
 */
long getProcessKilobytes()
{
long pageSize;

#ifdef _SC_PAGESIZE
  pageSize = sysconf(_SC_PAGESIZE);
#else
#warning "Page size query is not possible.  No per-process memory stats."
  return 0;
#endif

  pid_t pid = getpid();
  std::ostringstream fname;
  fname << "/proc/" << pid << "/statm";
  std::ifstream memFile;

  try{
    memFile.open(fname.str().c_str());
  }
  catch (...){
    return 0;
  }

  char buf[128];
  memset(buf, 0, 128);
  while (memFile.good()){
    memFile.getline(buf, 128);
    break;
  }

  memFile.close();

  std::istringstream sbuf(buf);
  long totalPages;
  sbuf >> totalPages;

  long pageKBytes = pageSize / 1024;
  totalPages = atol(buf);

  return totalPages * pageKBytes;
}



#ifdef HAVE_ZOLTAN2_MPI

/*! \brief Convert an MPI communicator to a MpiComm object.
 */

RCP<Teuchos::MpiComm<int> >
  MPI2Teuchos(const MPI_Comm &comm)
{
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
  RCP<mpiWrapper_t>handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);
  RCP<Teuchos::MpiComm<int> > tcommPtr(
    new Teuchos::MpiComm<int>(handle));

  return tcommPtr;
}

RCP<const Teuchos::MpiComm<int> >
  MPI2TeuchosConst(const MPI_Comm &comm)
{
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
  RCP<mpiWrapper_t>handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);
  RCP<const Teuchos::MpiComm<int> > tcommPtr(
    new Teuchos::MpiComm<int>(handle));

  return tcommPtr;
}

/*! \brief Convert a Teuchos::MpiComm object to the underlying MPI_Comm.
 */

MPI_Comm
  Teuchos2MPI(const RCP<Comm<int> > &comm)
{
  MPI_Comm mpiComm;
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;

  Comm<int> *c = comm.get();
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm = MPI_COMM_SELF;   // or would this be an error?
  }

  return mpiComm;
}

MPI_Comm
  TeuchosConst2MPI(const RCP<const Comm<int> > &comm)
{
  MPI_Comm mpiComm;
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;

  const Comm<int> *cConst = comm.get();
  Comm<int> *c = const_cast<Comm<int> *>(cConst);
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm = MPI_COMM_SELF;   // or would this be an error?
  }

  return mpiComm;
}

#endif

} // namespace Zoltan2
