#ifndef KOKKOS_BUFFER_MACROS
#define KOKKOS_BUFFER_MACROS

#include <Teuchos_TestForException.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <string>
#include <stdexcept>

#ifdef HAVE_KOKKOS_DEBUG
#define MARK_COMPUTE_BUFFER(buff) \
  { \
    std::string ptrtype("compute"); \
    Teuchos::set_extra_data(ptrtype, "BufferType", Teuchos::inOutArg(buff)); \
  }
#define CHECK_COMPUTE_BUFFER(buff) \
  { \
    TEST_FOR_EXCEPTION( Teuchos::get_extra_data<std::string>(buffSrc, "BufferType") != "compute", std::logic_error, \
        Teuchos::typeName(*this) << ": argument buffer is not a compute buffer as expected." ); \
  }
#else
#define MARK_COMPUTE_BUFFER(buff)
#define CHECK_COMPUTE_BUFFER(buff)
#endif

#endif
