#ifndef STK_SEARCH_BOOST_RTREE_INTERFACE_HPP
#define STK_SEARCH_BOOST_RTREE_INTERFACE_HPP

#ifdef __INTEL_COMPILER
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/utility/enable_if.hpp>
#else
// this disables the checking of shadowed variables on GCC only
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/utility/enable_if.hpp>
#pragma GCC diagnostic pop
#endif

#endif // STK_SEARCH_BOOST_RTREE_INTERFACE_HPP
