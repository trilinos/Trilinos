#ifndef MUELU_BOOSTGRAPHVIZ_HPP
#define MUELU_BOOSTGRAPHVIZ_HPP

// This header file can be used in place of <boost/graph/graphviz.hpp>. It disable the warnings present in boost.

// Note: pragma warnings available since gcc 4.2
//       pragma push/pop available since gcc 4.6

#include "MueLu_ConfigDefs.hpp"
#ifdef HAVE_MUELU_BOOST

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#endif // GCC_VERSION
#endif // __GNUC__

#include <boost/graph/graphviz.hpp>

#ifdef __GNUC__
# if ((__GNUC__ * 100) + __GNUC_MINOR__) >= 406
#pragma warning(pop)
#endif // GCC_VERSION
#endif // __GNUC__

#endif // HAVE_MUELU_BOOST

#endif // MUELU_BOOSTGRAPHVIZ_HPP
