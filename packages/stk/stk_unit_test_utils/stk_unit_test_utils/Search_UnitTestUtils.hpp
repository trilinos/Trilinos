// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_SEARCH_UNITTESTUTILS_HPP
#define STK_SEARCH_UNITTESTUTILS_HPP

#include <stk_util/environment/memory_util.hpp>
#include <stk_util/environment/WallTime.hpp>

#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoxIdent.hpp>

#include <stk_search/CoarseSearch.hpp>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommSparse.hpp>

#include <stk_unit_test_utils/getOption.h>

#include <iomanip>

typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<float> FloatBox;
typedef stk::search::Box<double> StkBox;

typedef int Ident;
typedef stk::search::IdentProc<int,int> IdentProc;

typedef std::pair<FloatBox, Ident> FloatBoxIdentPair;
typedef std::pair<StkBox,Ident> StkBoxIdentPair;
typedef std::pair<Sphere,Ident> SphereIdentPair;
typedef std::pair<FloatBox, IdentProc> FloatBoxIdentProcPair;
typedef std::pair<StkBox,IdentProc> StkBoxIdentProcPair;

typedef std::vector<FloatBoxIdentProcPair> FloatBoxIdentProcVector;
typedef std::vector<StkBoxIdentProcPair> StkBoxIdentProcVector;
typedef std::vector<FloatBoxIdentPair> FloatBoxIdentVector;
typedef std::vector<StkBoxIdentPair> StkBoxIdentVector;

typedef stk::search::BoxIdent<FloatBox, Ident> FloatBoxIdent;
typedef stk::search::BoxIdent<StkBox, Ident> StkBoxIdent;
typedef stk::search::BoxIdent<Sphere, Ident> SphereIdent;
typedef stk::search::BoxIdentProc<FloatBox, IdentProc> FloatBoxIdentProc;
typedef stk::search::BoxIdentProc<StkBox, IdentProc> StkBoxIdentProc;
typedef stk::search::BoxIdentProc<Sphere, IdentProc> SphereIdentProc;

typedef stk::search::IdentProcIntersection<IdentProc, IdentProc> IdentProcIntersection;
typedef stk::search::IdentIntersection<Ident, Ident> IdentIntersection;
typedef std::vector<std::pair<IdentProc,IdentProc>> SearchResults;
typedef std::vector<std::pair<Ident,Ident>> LocalSearchResults;

namespace stk {
namespace unit_test_util {

template<class VolumeType>
VolumeType generateBoundingVolume(double x, double y, double z, double radius);

template<>
inline
Point generateBoundingVolume<Point>(double x, double y, double z, double /*radius*/)
{
  return Point(x,y,z);
}

template<>
inline
Sphere generateBoundingVolume<Sphere>(double x, double y, double z, double radius)
{
  return Sphere(Point(x,y,z),radius);
}

//       ------------
//      |            |
//      |      radius|
//      |      ------|
//      |            |
//      |            |
//       ------------
// width = 2*radius
template<>
inline
StkBox generateBoundingVolume< StkBox >(double x, double y, double z, double radius)
{
  Point min_corner(x-radius,y-radius,z-radius);
  Point max_corner(x+radius,y+radius,z+radius);
  return StkBox(min_corner,max_corner);
}

template <typename VolumeType>
std::pair<VolumeType, IdentProc> generateBoundingVolume(double x, double y, double z, double radius, int id, int proc)
{
  return std::make_pair(generateBoundingVolume<VolumeType>(x,y,z,radius), IdentProc(id,proc));
}

template<class BoxType>
KOKKOS_FUNCTION
BoxType device_generateBox(double x, double y, double z, double radius);

template<>
KOKKOS_INLINE_FUNCTION
Point device_generateBox<Point>(double x, double y, double z, double /*radius*/)
{
  return Point(x,y,z);
}

template<>
KOKKOS_INLINE_FUNCTION
Sphere device_generateBox<Sphere>(double x, double y, double z, double radius)
{
  return Sphere(Point(x, y, z), radius);
}

template<>
KOKKOS_INLINE_FUNCTION
StkBox device_generateBox<StkBox>(double x, double y, double z, double radius)
{
  Point min_corner(x-radius, y-radius, z-radius);
  Point max_corner(x+radius, y+radius, z+radius);
  return StkBox(min_corner, max_corner);
}

template <typename BoxType, typename IdentProcType>
KOKKOS_FUNCTION
stk::search::BoxIdentProc<BoxType, IdentProcType> device_generateBoxIdentProc(double x, double y, double z,
                                                                              double radius, int id, int proc)
{
  return stk::search::BoxIdentProc<BoxType, IdentProcType>{device_generateBox<BoxType>(x, y, z, radius),
                                                           IdentProcType{id, proc}};
}

template <typename BoxType, typename IdentType>
KOKKOS_FUNCTION
stk::search::BoxIdent<BoxType, IdentType> device_generateBoxIdent(double x, double y, double z,
                                                                  double radius, IdentType id)
{
  return stk::search::BoxIdent<BoxType, IdentType>{device_generateBox<BoxType>(x, y, z, radius), id};
}

template <typename BoxIdent>
auto box_ident_to_pair(BoxIdent const& ident) {
  return std::make_pair(ident.box, ident.ident);
}

//======================

inline size_t getGoldValueForTest()
{
    std::string goldValue = stk::unit_test_util::get_option("-gold");
    if ( goldValue == "skip" ) return 0u;
    std::istringstream ss(goldValue);
    size_t goldValueNumber=0;
    ss >> goldValueNumber;
    return goldValueNumber;
}

inline void gatherResultstoProcZero(MPI_Comm comm, SearchResults& boxIdPairResults)
{
    int procId=-1;
    MPI_Comm_rank(comm, &procId);

    int numProc=-1;
    MPI_Comm_size(comm, &numProc);

    int procIdDestination = 0;
    stk::CommSparse commSparse(comm);
    for (int phase=0; phase<2; ++phase)
    {
        if ( procId != procIdDestination )
        {
            for (size_t j=0;j<boxIdPairResults.size();++j)
            {
                commSparse.send_buffer(procIdDestination).pack< std::pair<IdentProc, IdentProc> >(boxIdPairResults[j]);
            }
        }

        if (phase == 0) { //allocation phase
          commSparse.allocate_buffers();
        }
        else { // communication phase
          commSparse.communicate();
        }
    }

    if ( procId == procIdDestination )
    {
        for ( int p = 0 ; p < numProc ; ++p )
        {
            stk::CommBuffer &buf = commSparse.recv_buffer( p );
            while ( buf.remaining() )
            {
                std::pair<IdentProc, IdentProc> temp;
                buf.unpack< std::pair<IdentProc, IdentProc> >( temp );
                boxIdPairResults.push_back(temp);
            }
        }
    }
}

inline void printPeformanceStats(double elapsedTime, MPI_Comm comm)
{
    size_t maxHwm = 0, minHwm = 0, avgHwm = 0;
    stk::get_memory_high_water_mark_across_processors(comm, maxHwm, minHwm, avgHwm);

    int proc=-1;
    MPI_Comm_rank(comm, &proc);

    int numProcs=0;
    MPI_Comm_size(comm, &numProcs);

    double minTime = 0, maxTime = 0, avgTime = 0;
    MPI_Allreduce(&elapsedTime, &maxTime, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&elapsedTime, &minTime, 1, MPI_DOUBLE, MPI_MIN, comm);
    double elapsedTimeDivided = elapsedTime/numProcs;
    MPI_Allreduce(&elapsedTimeDivided, &avgTime, 1, MPI_DOUBLE, MPI_SUM, comm);

    if (proc == 0)
    {
      const double bytesInMegabyte = 1024.0*1024.0;
      std::cout << "Max time: "  << maxTime << ", Min time: " << minTime << ", Avg time: " << avgTime << std::endl;
      std::cout << std::setw(6) << std::fixed << std::setprecision(1) << "Max HWM: "<<maxHwm/bytesInMegabyte
        <<", Min HWM: "<<minHwm/bytesInMegabyte<<", Avg HWM: "<<avgHwm/bytesInMegabyte<<std::endl;
      std::cout<<"### Total Number of Steps Taken ###: 1"<<std::endl;
      std::cout<<"### Total Wall Clock Run Time Used ###: "<< maxTime <<std::endl;

      std::cout << "\nSTKPERF peak memory sum: " << maxHwm << std::endl;

    }
}

namespace simple_fields {

template<class VolumeType>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
VolumeType generateBoundingVolume(double x, double y, double z, double radius);

template<>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
inline
Point generateBoundingVolume<Point>(double x, double y, double z, double /*radius*/)
{
  return Point(x,y,z);
}

template<>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
inline
Sphere generateBoundingVolume<Sphere>(double x, double y, double z, double radius)
{
  return Sphere(Point(x,y,z),radius);
}

//       ------------
//      |            |
//      |      radius|
//      |      ------|
//      |            |
//      |            |
//       ------------
// width = 2*radius
template<>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
inline
StkBox generateBoundingVolume< StkBox >(double x, double y, double z, double radius)
{
  Point min_corner(x-radius,y-radius,z-radius);
  Point max_corner(x+radius,y+radius,z+radius);
  return StkBox(min_corner,max_corner);
}

template <typename VolumeType>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
std::pair<VolumeType, IdentProc> generateBoundingVolume(double x, double y, double z, double radius, int id, int proc)
{
  return std::make_pair(generateBoundingVolume<VolumeType>(x,y,z,radius), IdentProc(id,proc));
}


template<class BoxType>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
KOKKOS_FUNCTION
BoxType device_generateBox(double x, double y, double z, double radius);

template<>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
KOKKOS_INLINE_FUNCTION
Point device_generateBox<Point>(double x, double y, double z, double /*radius*/)
{
  return Point(x,y,z);
}

template<>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
KOKKOS_INLINE_FUNCTION
Sphere device_generateBox<Sphere>(double x, double y, double z, double radius)
{
  return Sphere(Point(x, y, z), radius);
}

template<>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
KOKKOS_INLINE_FUNCTION
StkBox device_generateBox<StkBox>(double x, double y, double z, double radius)
{
  Point min_corner(x-radius, y-radius, z-radius);
  Point max_corner(x+radius, y+radius, z+radius);
  return StkBox(min_corner, max_corner);
}

template <typename BoxType, typename IdentProcType>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
KOKKOS_FUNCTION
stk::search::BoxIdentProc<BoxType, IdentProcType> device_generateBoxIdentProc(double x, double y, double z,
                                                                              double radius, int id, int proc)
{
  return stk::search::BoxIdentProc<BoxType, IdentProcType>{device_generateBox<BoxType>(x, y, z, radius),
                                                           IdentProcType{id, proc}};
}

template <typename BoxType, typename IdentType>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
KOKKOS_FUNCTION
stk::search::BoxIdent<BoxType, IdentType> device_generateBoxIdent(double x, double y, double z,
                                                                  double radius, IdentType id)
{
  return stk::search::BoxIdent<BoxType, IdentType>{device_generateBox<BoxType>(x, y, z, radius), id};
}

template <typename BoxIdent>
STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
auto box_ident_to_pair(BoxIdent const& ident) {
  return std::make_pair(ident.box, ident.ident);
}

//======================

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
size_t getGoldValueForTest();

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void gatherResultstoProcZero(MPI_Comm comm, SearchResults& boxIdPairResults);

STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this function instead")
void printPeformanceStats(double elapsedTime, MPI_Comm comm);

} // namespace simple_fields

} // namespace unit_test_util
} // namespace stk

#endif
