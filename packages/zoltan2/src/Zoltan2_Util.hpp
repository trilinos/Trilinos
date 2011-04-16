// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_UTIL_HPP_
#define _ZOLTAN2_UTIL_HPP_

/*! \file Zoltan2_Util.hpp
*/

#include <exception>
#include <iostream>
#include <Zoltan2.hpp>

// Exception handling for client codes that call Zoltan2

#define HANDLE_ZOLTAN2_LIBRARY_EXCEPTIONS(msg) \
catch (const Zoltan2::memory_error &e){ \
  std::cerr << "Zoltan2 memory: "<< e.what() << std::endl; \
  return 1; \
} \
catch (const Zoltan2::usage_error &e){ \
  std::cerr << "Zoltan2 usage: "<< e.what() << std::endl; \
  return 1; \
} \
catch (const Zoltan2::data_type_error &e){ \
  std::cerr << "Zoltan2 data type: " << e.what() << std::endl; \
  return 1; \
} \
catch (const Zoltan2::possible_bug &e){ \
  std::cerr << "Zoltan2 bug: " << e.what() << std::endl; \
  return 1; \
} \
catch (const std::exception &e){ \
  std::cerr << "Caught by Zoltan2: " << e.what() << std::endl; \
  return 1; \
} \
catch (const Zoltan2::unknown_error_type &e){ \
  std::cerr << "Unknown error thrown to Zoltan2" << e.what() << std::endl; \
  return 1; \
}

// Exception handling for Zoltan2 when it calls standard library, Teuchos, etc.

#define Z2_HANDLE_UNKNOWN_EXCEPTIONS(msg) \
catch (const std::exception &e){ \
  std::cerr << __FILE__ << " (" << __LINE__ << ") " << msg << std::endl; \
  throw(e); \
} \
catch (...) { \
  throw(Zoltan2::unknown_error_type(msg, __FILE__, __LINE__)); \
}

// Exception handling for Zoltan2 when it calls within Zoltan2.  If
//   we're catching all exceptions we wouldn't need catch(...) at
//   the end.

#define Z2_HANDLE_MY_EXCEPTIONS \
catch (const Zoltan2::Zoltan2_error &e){ \
  std::cerr << __FILE__ << " (" << __LINE__ << ") " << std::endl; \
  throw(e); \
} \
catch (const std::exception &e){ \
  std::cerr << __FILE__ << " (" << __LINE__ << ") " << std::endl; \
  throw(e); \
} \
catch (...) { \
  throw(Zoltan2::unknown_error_type("", __FILE__, __LINE__)); \
}

#endif
