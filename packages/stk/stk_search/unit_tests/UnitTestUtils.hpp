// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <Geom_AxisAlignedBB.h>
#include <Geom_Search.h>
#include <search/ContactRangeSearch.h>
#include <search/ContactCommunication.h>

#include "gtkTraitsForSearch.h"
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/OctTreeOps.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

#include <optionParsing/getOption.h>

typedef stk::search::IdentProc<int,int> Ident;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> StkBox;
typedef std::pair<StkBox,Ident> StkBoxWithId;
typedef std::vector< StkBoxWithId > StkBoxVector;

typedef gtk::AxisAlignedBB GtkBox;
typedef std::vector<std::pair<Ident,Ident> > SearchResults;
typedef std::pair<GtkBox,Ident> BoxWithId;
typedef std::vector< BoxWithId > GtkBoxVector;


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
std::pair<VolumeType, Ident> generateBoundingVolume(double x, double y, double z, double radius, int id, int proc)
{
  return std::make_pair(generateBoundingVolume<VolumeType>(x,y,z,radius), Ident(id,proc));
}

//======================

inline size_t getGoldValueForTest()
{
    std::string goldValue = unitTestUtils::getOption("-gold");
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
    stk::CommAll gather(comm);
    for (int phase=0; phase<2; ++phase)
    {
        if ( procId != procIdDestination )
        {
            for (size_t j=0;j<boxIdPairResults.size();++j)
            {
                gather.send_buffer(procIdDestination).pack< std::pair<Ident, Ident> >(boxIdPairResults[j]);
            }
        }

        if (phase == 0) { //allocation phase
          gather.allocate_buffers( numProc / 4 );
        }
        else { // communication phase
          gather.communicate();
        }
    }

    if ( procId == procIdDestination )
    {
        for ( int p = 0 ; p < numProc ; ++p )
        {
            stk::CommBuffer &buf = gather.recv_buffer( p );
            while ( buf.remaining() )
            {
                std::pair<Ident, Ident> temp;
                buf.unpack< std::pair<Ident, Ident> >( temp );
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

#endif
