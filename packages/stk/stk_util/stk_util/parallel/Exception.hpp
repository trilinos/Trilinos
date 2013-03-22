/*--------------------------------------------------------------------*/
/*    Copyright 2003 - 2008 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_UTIL_PARALLEL_Exception_hpp
#define STK_UTIL_PARALLEL_Exception_hpp

/*--------------------------------------------------------------------*/
/** @file
 *  Exception Classes for SIERRA
 *
 *  This file defines all of the exception classes thrown in and
 *  managed by the SIERRA the framework.  A SIERRA exception class
 *  is defined for each of the C++ standard exceptions, except for
 *  ios_base::failure.  The ios_base::failure is omitted because
 *  it would require inclusion of the heavyweight <ios> header file.
 *
 *  The SIERRA exception types, their corresponding standard
 *  exception base type, and source of that standard exception
 *  base type are as follows.
 *
 *  It is important to note that the sierra exceptions do not follow the
 *  same hierachy as the std exceptions.  All the sierra execeptions are
 *  derived from ExParallel and not from each other.  Therefore, it is
 *  suggested that you catch only the std::exception & to catch all std
 *  and sierra or catch ExParallel & to catch all sierra exceptions.  The
 *  throw_copy() function promotes std exceptions to the corresponding
 *  sierra exception.
 *
 *  @li  sierra::Exception         std::exception          <exception>
 *  @li  sierra::BadException      std::bad_exception      <exception>
 *  @li  sierra::BadAlloc          std::bad_alloc          <new>
 *  @li  sierra::BadCast           std::bad_cast           <typeinfo>
 *  @li  sierra::BadTypeid         std::bad_typeid         <typeinfo>
 *  @li  sierra::LogicError        std::logic_error        <stdexcept>
 *  @li  sierra::DomainError       std::domain_error       <stdexcept>
 *  @li  sierra::InvalidArgument   std::invalid_argument   <stdexcept>
 *  @li  sierra::LengthError       std::length_error       <stdexcept>
 *  @li  sierra::OutOfRange        std::out_of_range       <stdexcept>
 *  @li  sierra::RuntimeError      std::runtime_error      <stdexcept>
 *  @li  sierra::RangeError        std::range_error        <stdexcept>
 *  @li  sierra::OverflowError     std::overflow_error     <stdexcept>
 *  @li  sierra::UnderflowError    std::underflow_error    <stdexcept>
 *
 *  @li  typedef ExTemp1< std::ios_base::failure >  IosBaseFailure;
 *  IosBaseFailure is declared in ExceptionIos.h in order to
 *  avoid including the ios header file in this header.
 *
 *  Any additional exception types should be derived from one of the
 *  above SIERRA exception classes.
 *
 * To create a new exception, DO NOT subclass from a parallel exception.
 * Create your new exception from derived from the STL exception class,
 * then typedef your new exception to a parallel exception using the
 * ExTemp or ExTemp1 template.  Then, be sure to call the
 * MyException::registerException() during application initialization.
 */

#include <stdexcept>
#include <exception>
#include <new>
#include <typeinfo>
#include <string>
#include <vector>
// #include <iostream> // for std::cerr
#include <sstream>

#include <stk_util/stk_config.h>
#if defined( STK_HAS_MPI )
#include <mpi.h>
#endif

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_util/util/FeatureTest.hpp>
#include <stk_util/diag/String.hpp>
#include <stk_util/diag/StringUtil.hpp>
#include <stk_util/diag/Trace.hpp>

namespace sierra {

///
/// @addtogroup ExceptionDetail
/// @{
///

class ExParallel;

/**
 * @brief Member function <b>register_stl_parallel_exceptions</b> registers the stl
 * exceptions with the parallel exception handler.
 *
 */
void register_stl_parallel_exceptions();

/**
 * @brief Function <b>set_exception</b> is called on a single processor when an
 * exception is caught.  The next collective communication will propogate the exception to
 * all processors.  This flavor is called when an unknown exception (...) is caught.
 *
 */
void set_exception();

/**
 * @brief Function <b>set_exception</b> is called on a single processor when an
 * exception is caught.  The next collective communication will propogate the exception to
 * all processors.
 *
 * @param x		a <b>std::exception</b> reference to the base exception.
 */
void set_exception(std::exception &x);

/**
 * @brief Function <b>set_exception</b> is called on a single processor when an
 * exception is caught.  The next collective communication will propogate the exception to
 * all processors.
 *
 * @param x		a <b>std::exception</b> reference to the parallel
 *			exception.
 */
void set_exception(ExParallel &x);

/**
 * @brief Member function <b>sierra_exception_throw</b> is called whenever a
 * parallel exception is constructed.  This acn be used as a breakpoint in a debugger if
 * finding the magic throw symbol is too tricky.
 *
 */
void sierra_exception_throw();

/**
 * @brief Function <b>throw_copy</b> throws a copy of the exception.  The exception
 * is located in the parallel exception registry.  This exception is cloned and thrown
 * using the parallel exception's <b>throw_copy</b> member function.
 *
 * If the exception cannot be located, an <b>Exception</b> parallel exception is
 * thrown instead contained the message from the original exception from
 * <b>what()</b>.
 *
 * @param x		a <b>std::exception</b> reference to the exception to be
 *			copied and thrown.
 *
 * @param message	a <b>std::string</b> const reference to a message to be
 *			appended to the exception's description.
 *
 */
void throw_copy(const std::exception &x, const std::string &message);

/**
 * @brief Function <b>parallel_throw</b> throws a consistant exception in
 * parallel.  parallel_throw is called after the exception has been propagated to all
 * processors.  parallel_throw will decide which exception to rethrow in parallel on all
 * processors.
 *
 * @param mpi_comm	a <b>MPI_Comm</b> value of the communicator to talk over.
 *
 */
void parallel_throw(MPI_Comm mpi_comm);


/**
 * @brief Class <b>ParallelThrowRegistry</b> is a registry of known parallel
 * exceptions.  For the negotiation of parallel exceptions, each parallel exception must
 * be registered.  The class serves as the registry and registrar for the exceptions.
 *
 * When an exception is registered, a new exception of that type is provided which is
 * cloned whenever and exception of that type is required.
 *
 * Two different exception types are registered.  The parallel exception type and the base
 * exception type which is discovered using the BaseExceptionType type of the parallel
 * exception class.  When the parallel exception registry is searched, the parallel
 * exception to be cloned is returned when either type is found.
 *
 * During parallel exception processing, the information from exception to be thrown in
 * parallel is stored in the registered exception object.
 *
 */
class ParallelThrowRegistry
{
public:
  /**
   * @brief Member function <b>instance</b> returns the singleton instance for the
   * parallel exception registry.
   *
   * @return a <b>ParallelThrowRegistry</b> ...
   */
  static ParallelThrowRegistry &instance();

  /**
   * @brief Static member function <b>setExceptionsRegistered</b> sets the
   * exceptions have been registered flag.  This flag is used to control the calling of
   * the function <b>sierra_exception_thrown</b> when a parallel exception is
   * created.
   *
   */
  static bool areExceptionsRegistered() {
    return !ParallelThrowRegistry::instance().m_registry.empty();
  }

  /**
   * @brief Member template function <b>registerException</b> registers an exception
   * of the specified type with the parallel registry.  An exception is provided which is
   * to be cloned whenever a parallel exception of that type is the be created for
   * parallel throwing.
   *
   * @return			an <b>ExParallel</b> reference to the
   *				<b>x_to_clone</b> exception.
   */
  template<class T>
  ExParallel &registerException() {
    T *x = new T();
    register_exception_a(typeid(T), x);
    register_exception_a(typeid(typename T::BaseExceptionType), x);
    return *x;
  }

  /**
   * @brief Member function <b>findException</b> returns a pointer to the matching exception
   * in parallel exception registry.
   *
   * @param exception_type	a <b>std::type_info</b> reference to the exception
   *				type being searched.
   *
   * @return			an <b>ExParallel</b> pointer to the exception to be cloned.
   */
  ExParallel *findException(const std::type_info &exception_type);

private:
  /**
   * @brief Member function <b>register_exception_a</b> performs the actual registration.
   *
   * @param exception_type		a <b>std::type_info</b> reference of the
   *					base exception type obtained via typeid of the
   *					BaseExceptionType type of the parallel exception type.
   *
   * @param exception			an <b>ExParallel</b> pointer to the parallel
   *					exception to be cloned.
   *
   */
  ExParallel &register_exception_a(const std::type_info &parallel_exception_type, ExParallel *exception);

private:
  /**
   * @brief Typedef <b>Registry</b> is a registry of pairing of each exception type
   * to its registered parallel throw exception to be cloned.
   *
   */
  class Registry : public std::vector<std::pair<const std::type_info *, ExParallel *> > 
  {
  public:
    Registry();

    ~Registry();
  };

  Registry		m_registry;			///< Registry of exceptions
};


/**
 * @brief Class <b>ExParallel</b> implements the features of a parallel exception.
 * It is a <b>std::string</b> which stores the exception description.  It also
 * provides "put to" (operator;lt&;lt&) functions for appending the description with
 * information.
 *
 */
class ExParallel
{
protected:
  /**
   * Creates a new <b>ExParallel</b> instance which is not about to be thrown in
   * parallel and has no description.
   *
   */
  ExParallel()
    : m_descriptionStream(),
      m_whatBuffer(),
      m_traceback(Diag::Trace::printTraceback(Diag::Trace::Traceback::snapshot())),
      m_parallel(-1)
  {
    if (ParallelThrowRegistry::areExceptionsRegistered())
      sierra_exception_throw();
  }

  /**
   * Creates a new <b>ExParallel</b> instance with an initial description and may be
   * about to be thrown in parallel.
   *
   * @param message		an <b>std::string</b> const reference to the initial
   *				description.
   *
   * @param parallel		an <b>bool</b> value of true if this exception is
   *				about to be thrown in parallel.
   */
  explicit ExParallel(const std::string & message, int parallel = -1)
    : m_descriptionStream(),
      m_whatBuffer(),
      m_traceback(Diag::Trace::printTraceback(Diag::Trace::Traceback::snapshot())),
      m_parallel(parallel)
  {
    m_descriptionStream << message;
    if (ParallelThrowRegistry::areExceptionsRegistered())
      sierra_exception_throw();
  }

  /**
   * Creates a new <b>ExParallel</b> copy.
   *
   * @param x			an <b>ExParallel</b> const reference to the
   *				exception to copy.
   */
  ExParallel(const ExParallel &x)
    : m_descriptionStream(),
      m_whatBuffer(),
      m_traceback(x.m_traceback),
      m_parallel(x.m_parallel)
  {
    m_descriptionStream << x.m_descriptionStream.str();
  }

public:
  /**
   * Destroys a <b>ExBase</b> instance.
   *
   */
  virtual ~ExParallel()
  {}

  /**
   * @brief Member function <b>what</b> returns the exception's description.
   *
   * @return			a <b>char</b> const pointer to the exception's
   *				description.
   */
  virtual const char *what() const throw() {
    try {
      m_whatBuffer = m_descriptionStream.str();
      return m_whatBuffer.c_str();
    }
    catch(...) {
      return NULL;
    }
  }

  /**
   * @brief Member function <b>clear</b> clears the contents of the exception.
   *
   * @return			an <b>ExParallel</b> reference to the exception.
   */
  ExParallel &clear() {
    m_descriptionStream.str("");
    m_traceback.clear();
    m_parallel = -1;
    return *this;
  }

  /**
   * @brief Member function <b>setDescription</b> sets the value of the exception's
   * description.
   *
   * @param description		a <b>std::string</b> const reference to the
   *				description.
   *
   * @return			an <b>ExParallel</b> reference to the exception.
   */
  ExParallel &setDescription(const std::string &description) {
    clear();
    m_descriptionStream << description;
    return *this;
  }

  /**
   * @brief Member function <b>getDescription</b> returns the exception's
   * description.
   *
   * @return			a <b>std::string</b> value of the exception's
   *				description.
   */
  std::string getDescription() const{
    return m_descriptionStream.str();
  }

  /**
   * @brief Member function <b>getDescriptionStream</b> returns the stream used to
   * assemble the description.
   *
   * @return			a <b>std::ostringstream</b> reference to the
   *				exception's description stream.
   */
  std::ostringstream &getDescriptionStream() {
    return m_descriptionStream;
  }

  /**
   * @brief Member function <b>getDescriptionStream</b> returns the stream used to
   * assemble the description.
   *
   * @return			a <b>std::ostringstream</b> const reference to the
   *				exception's description stream.
   */
  const std::ostringstream &getDescriptionStream() const {
    return m_descriptionStream;
  }

  /**
   * @brief Member function <b>setTraceback</b> sets the exception's traceback
   * to the caller generating the exception.
   *
   * @param traceback		a <b>std::string</b> const reference to the
   *				traceback.
   *
   * @return			an <b>ExParallel</b> reference to the exception.
   */
  ExParallel &setTraceback(const std::string &traceback) {
    m_traceback = traceback;
    return *this;
  }

  /**
   * @brief Member function <b>getTraceback</b> returns the exception's traceback string.
   *
   * @return			a <b>std::string</b> const reference to the
   *				exception's traceback string.
   */
  const std::string &getTraceback() const {
    return m_traceback;
  }

  /**
   * @brief Member function <b>setParallel</b> sets the originating processor for an
   * exception that is being thrown in parallel.
   *
   * @param parallel			a <b>int</b> value to set the originating
   *					processor for this parallel exception.
   *
   * @return				an <b>ExParallel</b> reference to this
   *					parallel exception.
   */
  ExParallel &setParallel(int parallel) {
    m_parallel = parallel;
    return *this;
  }

  /**
   * @brief Member function <b>getParallel</b> returns the originating processor for
   * a parallel exception or -1 if the exception was not thrown in parallel.
   *
   * @return				a <b>int</b> value of the originating processor
   *					for the parallel exception or -1.
   */
  int getParallel() const {
    return m_parallel;
  }


  /**
   * @brief Member function <b>isParallel</b> returns true if the exception is
   * being thrown in parallel.
   *
   * @return				a <b>bool</b> value of true if the exception
   *					is being thrown in parallel.
   */
  bool isParallel() const {
    return m_parallel != -1;
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the std manipilator
   * functions to the <b>ExParallel</b> object.  This allows the manipulators to
   * modify the description.  (Currently in a limited fashion).
   *
   * @return				an <b>ExParallel</b> reference to this
   *					parallel exception.
   */
  ExParallel &operator<<(std::ostream& (*f)(std::ostream&)) {
    f(m_descriptionStream);
    return *this;
  }

  /**
   * Member function <b>operator&lt;&lt;</b> passes any data type to the exception
   * string class for conversion to a string.
   *
   * @param t			a <b>T</b> const reference that is to be converted
   *				to a string.
   *
   * @return			a <b>ExParallel</b> reference to this object;
   */
  template<class U>
  ExParallel &operator<<(const U &t) {
    m_descriptionStream << t;
    return *this;
  }

  /**
   * @brief Member function <b>throw_copy</b> is a pure virtual function which is allows
   * the copying and throwing of the parallel exception using a pointer to the base
   * exception object.  This function should never return, it should throw the copied
   * exception.
   *
   * The exception is being throw in parallel so be sure to set the parallel thrown
   * processor via <b>setParallel(int)</b> before throwing the copy.
   *
   */
  virtual void throw_copy() const = 0;

  /**
   * @brief Member function <b>parallel_handler</b> is called just before a parallel
   * exception is thrown.  It is guaranteed to be called in parallel on all processors, so
   * collective communication is allowed inside Parallel_Handler.  This function might be
   * used to copy information to all processors.  The default is to do nothing.
   *
   */
  virtual void parallel_handler();

private:
  std::ostringstream	m_descriptionStream;		///< Description from original throw
  mutable std::string	m_whatBuffer;			///< what() buffer
  std::string		m_traceback;			///< Traceback from original throw
  int			m_parallel;			///< True if being thrown in parallel
};


/**
 * @brief Template <b>ExTemp</b> takes a zero argument exception and makes it into a
 * parallel throwable and put-to-able (<<) exception.  This exception may be caught with
 * either the base class <b>T</b> type or the template <b>ExTemp&lt;T&gt;</b>
 * type.
 *
 */
template<class T>
class ExTemp : public ExParallel, public T
{
public:
  typedef ExTemp<T> ParallelExceptionType;		///< Parallel exception type
  typedef T BaseExceptionType;				///< Base exception type

  /**
   * Creates a new <b>ExTemp</b> instance.
   *
   */
  ExTemp()
    : ExParallel(),
      T()
  {}

  /**
   * Creates a new <b>ExTemp</b> instance with an initial description.
   *
   * @param message		a <b>std::string</b> const reference to the initial
   *				exception description.
   */
  explicit ExTemp(const std::string &message)
    : ExParallel(message),
      T()
  {}

  /**
   * Creates a new <b>ExTemp</b> copy.
   *
   * @param x			an <b>ExTemp</b> variable ...
   */
  ExTemp(const ExTemp & x)
    : ExParallel(static_cast<const ExParallel &>(x)),
      T(static_cast<const T &>(x))
  {}

  /**
   * Destroys a <b>ExTemp</b> instance.
   *
   */
  virtual ~ExTemp() throw()
  {}

  /**
   * @brief Member function <b>what</b> returns the exception's description.
   *
   * @return			a <b>char</b> const pointer to the exception's
   *				description.
   */
  virtual const char *what() const throw() {
    return ExParallel::what();
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the std manipilator
   * functions to the <b>ExTemp</b> object.  This allows the manipulators to modify the
   * description.  (Currently in a limited fashion).
   *
   * @return				an <b>ExParallel</b> reference to this
   *					parallel exception.
   */
  ExTemp &operator<<(std::ostream& (*f)(std::ostream&)) {
    f(getDescriptionStream());
    return *this;
  }

  /**
   * Member function <b>operator&lt;&lt;</b> passes any data type to the exception
   * string class for conversion to a string.
   *
   * @param t			a <b>T</b> const reference that is to be converted
   *				to a string.
   *
   * @return			a <b>ExTemp</b> reference to this object;
   */
  template<class U>
  ExTemp &operator<<(const U &t) {
    getDescriptionStream() << t;
    return *this;
  }

  /**
   * @brief Member function <b>copy</b> throws a copy of the original exception.  It
   * copies the original message, sets the parallel thrown flag and throws the new
   * exception.
   *
   */
  virtual void throw_copy() const {
    ParallelExceptionType t(*this);

//    std::cerr << "throwing " << this->what() << std::endl
//              << "      as " << t.what();
    throw t;
  }

  /**
   * @brief Member function <b>registerException</b> registers the exception with
   * the parallel exception registry.
   *
   * @return			an <b>ExParallel</b> reference to the exception to
   *				be cloned when thrown.
   */
  static void registerException() {
#ifdef SIERRA_TEMPLATE_CALL_BUG
    ParallelThrowRegistry::instance().template registerException<ExTemp>();
#else
    ParallelThrowRegistry::instance().registerException<ExTemp>();
#endif
  }
};


template<class T>
class ExTemp1 : public ExParallel, public T
{
public:
  typedef ExTemp1<T> ParallelExceptionType;		///< Parallel exception type
  typedef T BaseExceptionType;				///< Base exception type

  /**
   * Creates a new <b>ExTemp1</b> instance.
   *
   */
  ExTemp1()
    : ExParallel(),
      T(std::string())
  {}

  /**
   * Creates a new <b>ExTemp1</b> instance with an initial description.
   *
   * @param message		a <b>std::string</b> const reference to the initial
   *				exception description.
   */
  explicit ExTemp1(const std::string & message)
    : ExParallel(message),
      T(message)
  {}

  /**
   * Creates a new <b>ExTemp</b> copy.
   *
   * @param x			an <b>ExTemp</b> variable ...
   */
  ExTemp1(const ExTemp1 & x)
    : ExParallel(static_cast<const ExParallel &>(x)),
      T(static_cast<const T &>(x))
  {}

  /**
   * Destroys a <b>ExTemp1</b> instance.
   *
   */
  virtual ~ExTemp1() throw()
  {}

  /**
   * @brief Member function <b>what</b> returns the exception's description.
   *
   * @return			a <b>char</b> const pointer to the exception's
   *				description.
   */
  virtual const char *what() const throw() {
    return ExParallel::what();
  }

  /**
   * @brief Member function <b>operator&lt;&lt;</b> passes the std manipilator
   * functions to the <b>ExTemp1</b> object.  This allows the manipulators to modify the
   * description.  (Currently in a limited fashion).
   *
   * @return				an <b>ExParallel</b> reference to this
   *					parallel exception.
   */
  ExTemp1 &operator<<(std::ostream& (*f)(std::ostream&)) {
    f(getDescriptionStream());
    return *this;
  }


  /**
   * Member function <b>operator&lt;&lt;</b> passes any data type to the exception
   * string class for conversion to a string.
   *
   * @param t			a <b>T</b> const reference that is to be converted
   *				to a string.
   *
   * @return			a <b>ExTemp1</b> reference to this object;
   */
  template<class U>
  ExTemp1 &operator<<(const U &t) {
    getDescriptionStream() << t;
    return *this;
  }

  /**
   * @brief Member function <b>throw_copy</b> throws a copy of the original
   * exception.  It copies the original message, appends the message to the exception
   * description, sets the parallel thrown flag and throws the new exception.
   *
   * @param message		a <b>std::string</b> const reference to a message to
   *				be appended to the exception's description prior to being
   *				thrown.
   */
  virtual void throw_copy() const {
    ParallelExceptionType t(*this);

//    std::cerr << "throwing " << this->what() << std::endl
//              << "      as " << t.what();
    throw t;
  }

  /**
   * @brief Member function <b>registerException</b> registers the exception with
   * the parallel exception registry.
   *
   * @return			an <b>ExParallel</b> reference to the exception to
   *				be cloned when thrown.
   */
  static ExParallel &registerException() {
#ifdef SIERRA_TEMPLATE_CALL_BUG
    return ParallelThrowRegistry::instance().template registerException<ExTemp1>();
#else
    return ParallelThrowRegistry::instance().registerException<ExTemp1>();
#endif
  }
};

typedef ExTemp<std::exception>		Exception;		///< Defined in <exception>
typedef ExTemp<std::bad_exception>	BadException;		///< Defined in <exception>
typedef ExTemp<std::bad_alloc>		BadAlloc;		///< Defined in <new>
typedef ExTemp<std::bad_typeid>		BadTypeid;		///< Defined in <typeinfo>
typedef ExTemp<std::bad_cast>		BadCast;		///< Defined in <typeinfo>
typedef ExTemp1<std::ios_base::failure>	IosBaseFailure;		///< Defined in <ios>
typedef ExTemp1<std::logic_error>	LogicError;		///< Defined in <stdexcept>
typedef ExTemp1<std::domain_error>	DomainError;		///< Defined in <stdexcept>
typedef ExTemp1<std::invalid_argument>	InvalidArgument;	///< Defined in <stdexcept>
typedef ExTemp1<std::length_error>	LengthError;		///< Defined in <stdexcept>
typedef ExTemp1<std::out_of_range>	OutOfRange;		///< Defined in <stdexcept>
typedef ExTemp1<std::runtime_error>	RuntimeError;		///< Defined in <stdexcept>
typedef ExTemp1<std::range_error>	RangeError;		///< Defined in <stdexcept>
typedef ExTemp1<std::overflow_error>	OverflowError;		///< Defined in <stdexcept>
typedef ExTemp1<std::underflow_error>	UnderflowError;		///< Defined in <stdexcept>

//----------------------------------------------------------------------

/**
 * @brief Class <b>runtime_user_error</b> ...
 *
 */
class runtime_user_error : public std::runtime_error
{
public:
  explicit runtime_user_error(const std::string &message) throw()
   : std::runtime_error(message)
  {}

  runtime_user_error(const runtime_user_error &x) throw()
    : std::runtime_error(x)
  {}

  virtual ~runtime_user_error() throw ()
  {}
};

typedef ExTemp1<runtime_user_error> RuntimeUserError;

///
/// @}
///

} // namepace sierra

// DO NOT USE ParallelStackTrace, unless you know that the exception will be thrown in
// parallel.  This is a magic string used by the parallel_throw routine that makes
// sure the printout will only be done on processor 0.  Since the exception is
// parallel there is no reason to print the same message many times.
#define StackTraceMessage "  exception thrown from "
#define ParallelStackTraceMessage "  parallel exception thrown from "
#define ParallelStackTrace std::string(std::string(ParallelStackTraceMessage) + stk::source_relative_path(STR_TRACE))

#endif // STK_UTIL_PARALLEL_Exception_hpp
