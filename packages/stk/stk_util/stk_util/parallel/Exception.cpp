/**   ------------------------------------------------------------
 *    Copyright 2003 - 2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <stdexcept>
#include <exception>
#include <new>
#include <typeinfo>
#include <ios>
#include <string>
#include <sstream>
#include <iostream>

#include <assert.h>

#include <stk_util/diag/Env.hpp>
#include <stk_util/diag/Platform.hpp>
#include <stk_util/parallel/Exception.hpp>
#include <stk_util/parallel/ExceptionReport.hpp>
#include <stk_util/parallel/ExceptionIos.hpp>
#include <stk_util/diag/String.hpp>
#include <stk_util/diag/Trace.hpp>

#include <stk_util/parallel/mpih.hpp>

namespace sierra {

void
sierra_exception_throw()
{}


ParallelThrowRegistry &
ParallelThrowRegistry::instance()
{
  static ParallelThrowRegistry s_parallelThrowRegistry;

  return s_parallelThrowRegistry;
}


ParallelThrowRegistry::Registry::Registry()
{}


ParallelThrowRegistry::Registry::~Registry()
{
  // Truely sick.  Each is registered twice, once for the parallel version of the
  // exception and once for the <stdexcept> base class version.  The double increment
  // keeps from deleting it twice.  See ParallelThrowRegistry::registerException.
  for (iterator it = begin(); it != end(); ++it, ++it)
    delete (*it).second;
}


ExParallel &
ParallelThrowRegistry::register_exception_a(
  const std::type_info &	exception_type,
  ExParallel *			exception)
{
  if (!findException(exception_type)) {
    m_registry.push_back(Registry::value_type(&exception_type, exception));
    mpih::Add_Handle(*exception);
  }
  return *exception;
}

ExParallel *
ParallelThrowRegistry::findException(
  const std::type_info &		exception_type)
{
  for (Registry::iterator it = m_registry.begin(); it != m_registry.end(); ++it)
    if (*(*it).first == exception_type)
      return (*it).second;

  return NULL;
}


void
ExParallel::parallel_handler()
{}


void
throw_copy(
  const std::exception &	x,
  const std::string &		append_message)
{
  ExParallel *exception = ParallelThrowRegistry::instance().findException(typeid(x));
  if (!exception)
    exception = ParallelThrowRegistry::instance().findException(typeid(Exception));

  exception->clear();
  *exception << x.what() << append_message;

  exception->throw_copy();
}


void
set_exception()
{
  BadException x;
  x << "Unknown exception";
  set_exception(static_cast<ExParallel &>(x));
}


void
set_exception(
  std::exception &		x)
{
  ExParallel *registered_exception = ParallelThrowRegistry::instance().findException(typeid(x));

  if (!registered_exception)
    registered_exception = ParallelThrowRegistry::instance().findException(typeid(Exception));

  registered_exception->setDescription(x.what());
  registered_exception->setTraceback(Diag::Traceback::printTraceback(Diag::Traceback::snapshot()));

//  std::cerr << "Exception " << demangle(typeid(*registered_exception).name()) << " will be thrown from processor " << Env::parallel_rank() << " on the next MPIH function:" << std::endl
//	    << registered_exception->getDescription() << std::endl
//	    << registered_exception->getTraceback() << std::endl;

  mpih::Set_Local_Handle(const_cast<ExParallel &>(*registered_exception));
}


void
set_exception(
  ExParallel &			x)
{
  ExParallel *registered_exception = ParallelThrowRegistry::instance().findException(typeid(x));

  if (!registered_exception)
    registered_exception = ParallelThrowRegistry::instance().findException(typeid(Exception));

  registered_exception->setDescription(x.getDescription());
  registered_exception->setTraceback(Diag::Traceback::printTraceback(Diag::Traceback::snapshot()));

//  std::cerr << "Exception " << demangle(typeid(*registered_exception).name()) << " will be thrown from processor " << Env::parallel_rank() << " on the next MPIH function:" << std::endl
//	    << registered_exception->getDescription() << std::endl
//	    << registered_exception->getTraceback() << std::endl;

  mpih::Set_Local_Handle(const_cast<ExParallel &>(*registered_exception));
}


void
register_stl_parallel_exceptions()
{
  mpih::Enable();

  Exception::registerException();
  BadAlloc::registerException();
  BadCast::registerException();
  BadTypeid::registerException();
  LogicError::registerException();
  DomainError::registerException();
  InvalidArgument::registerException();
  LengthError::registerException();
  OutOfRange::registerException();
  RuntimeError::registerException();
  RangeError::registerException();
  OverflowError::registerException();
  UnderflowError::registerException();
  BadException::registerException();

  mpih::Activate_Handles();
}


void
parallel_throw(
  MPI_Comm		mpi_comm)
{
  int nprocs;
  MPI_Comm_size(mpi_comm, &nprocs);

  ExParallel **handles = new ExParallel* [nprocs];

  mpih::Get_Global_Handles(handles);

  MPIH_Handler_compete handler_compete_fn;
  MPIH_Handler_execute handler_execute_fn;
  mpih::Get_Functions(&handler_compete_fn ,
		      &handler_execute_fn);

  /* Now that we have the handles,
   * reset the handles so we don't throw again.  This way
   * whatever function catches the exception we are about to
   * throw can call mpih and mpih will not just throw again.
   */
  mpih::Reset_Local_Handle();

  /* First iterate through all of the exceptions thrown on all of
   * the processors, and if any of them were thrown on this
   * processor, print an error message and a traceback.
   * only the owning processor will have the traceback information.
   */
  /* Iterate through all of the exceptions thrown on all of the processors
   * and call the parallel_handler() function defined by any derived from
   * ExParallel.  This is done across all processors so that
   * it is valid to do collective communication inside of parallel_handler()
   */
  for (int i = 0; i < nprocs; ++i) {
    if (handles[i]) {
      ExParallel *x = dynamic_cast<ExParallel *>(handles[i]);
      if (x)
	x->parallel_handler();
    }
  }

  /* Iterate through all of the exceptions thrown on all of the processors
   * and select the one to throw in parallel on all processors.  We would
   * like to find one derived from ExParallel.
   */

  ExParallel *the_exception = NULL;
  int originating_processor = -1;

  for (int i = 0; i < nprocs; ++i) {
    if (handles[i]) {
      ExParallel *x = dynamic_cast<ExParallel *>(handles[i]);
      if (x) {
	if (handler_compete_fn)
	  (handler_compete_fn) (reinterpret_cast<void **>(&handles[i]), the_exception);
	if ( handles[i] != the_exception ) {
	  the_exception = x;
	  originating_processor = i;
	}
      }
    }
  }

  delete [] handles;

  /* Since this function is called in parallel, it is possible
   * to perform collective communication.  Here the traceback
   * and error messages are broadcast and set on all processors.
   * These are the only two fields that are guarenteeded to be
   * in each exception class.  Other data stored in specialized
   * derived classes will have to be communicated seperately.
   * If needed this communication could be added to a virtual
   * base class.  That is a future enhancements depending on
   * the demand.
   */
  if (the_exception) {
    // Copy the description from the originating process to everywhere.
    std::string description(the_exception->getDescriptionStream().str());
    int description_len = description.length();
    MPI_Bcast(&description_len,
	      1,
	      MPI_INT,
	      originating_processor,
	      mpi_comm);

    char *description_buf = new char[description_len];
    description.copy(description_buf, description_len);

    MPI_Bcast(description_buf,
	      description_len,
	      MPI_CHAR,
	      originating_processor,
	      mpi_comm);

    // Copy the traceback stack from the originating process to everywhere.
    const std::string &traceback(the_exception->getTraceback());
    int traceback_len = traceback.length();
    MPI_Bcast(&traceback_len,
	      1,
	      MPI_INT,
	      originating_processor,
	      mpi_comm);

    char *traceback_buf = new char[traceback_len];
    traceback.copy(traceback_buf, traceback_len);

    MPI_Bcast(traceback_buf,
	      traceback_len,
	      MPI_CHAR,
	      originating_processor,
	      mpi_comm);

    // Rebuild the exception from the broadcasted data
    the_exception->setDescription(std::string(description_buf, description_len));
    the_exception->setTraceback(std::string(traceback_buf, traceback_len));
    the_exception->setParallel(originating_processor);

//     std::cerr << "Throwing exception " << demangle(typeid(*the_exception).name()) << " in parallel" << std::endl
//               << the_exception->getDescription() << std::endl
//               << the_exception->getTraceback() << std::endl;

#ifdef SIERRA_MPIH_VERBOSE
    Env::outputP0()
      <<"*************** Exception handling ***************"<<endl
      <<" A parallel exception of type "<< typeid(*the_exception).name()<<endl
      <<" will be thrown on all processors."<<endl;
#endif

    delete [] traceback_buf;
    delete [] description_buf;
    the_exception->throw_copy();
  }
  else {
#ifdef SIERRA_MPIH_VERBOSE
    Env::outputP0()
      <<"*************** Exception handling ***************"<<endl
      <<" A parallel exception of type Unknown_Exception"<<endl
      <<" will be thrown on all processors."<<endl;
#endif
    throw Exception();
  }
}

} // namespace sierra
