/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions Contact: William McLendon (wcmclen@sandia.gov) or
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using std::cout;
using std::vector;

//#include "../../src/sparse/impl/KokkosSparse_partitioning_impl.hpp"

int main(int argc, char *argv[])
{
  /*
    const int device_id   = 0;

    Kokkos::initialize(Kokkos::InitArguments(1, -1, device_id));

    if(params.verbose)
    {
      Kokkos::print_configuration(std::cout);
    }
    */

    //Generate a square 2D mesh in serial, with edges in both diagonals
    const int xy = 31;
    const int nnodes = (xy + 1) * (xy + 1);
    vector<bool> adjMat(nnodes * nnodes, false);
    //Number of nodes is (n+1)^2
    for(int cellX = 0; cellX < xy; cellX++)
    {
      for(int cellY = 0; cellY < xy; cellY++)
      {
        int upLeft =    cellX     + (xy + 1) * cellY;
        int upRight =   cellX + 1 + (xy + 1) * cellY;
        int downLeft =  cellX     + (xy + 1) * (cellY + 1);
        int downRight = cellX + 1 + (xy + 1) * (cellY + 1);
        #define CONNECT(n1, n2) \
          adjMat[n1 + n2 * nnodes] = true; \
          adjMat[n2 + n1 * nnodes] = true;
        //Form this pattern in each cell:
        //
        //    +------+
        //    |\    /|
        //    | \  / |
        //    |  \/  |
        //    |  /\  |
        //    | /  \ |
        //    |/    \|
        //    +------+
        //
        CONNECT(upLeft, upRight);
        CONNECT(upLeft, downLeft);
        CONNECT(upLeft, downRight);
        CONNECT(upRight, downRight);
        CONNECT(downLeft, downRight);
        CONNECT(downLeft, upRight);
      }
    }
    
    //Build a sparse (CRS) graph from the dense adjacency matrix
    int numEdges = 0;
    for(size_t i = 0; i < adjMat.size(); i++)
      numEdges += (adjMat[i] ? 1 : 0);

    /*
    Kokkos::View<int*, Kokkos::HostSpace> rowmap("Rowmap", nnodes + 1);
    Kokkos::View<int*, Kokkos::HostSpace> entries("Entries", numEdges);
    int accum = 0;
    for(int r = 0; r <= nnodes; r++)
    {
      rowmap(r) = accum;
      if(r == nnodes)
        break;
      for(int c = 0; c < nnodes; c++)
      {
        if(adjMat[c + r * nnodes])
          entries(accum++) = c;
      }
    }
    */

    //Dump the graph to a graphviz file
    FILE* g = fopen("graph.dot", "w");
    fprintf(g, "graph {\n");
    for(int r = 0; r < nnodes; r++)
    {
      for(int c = r; c < nnodes; c++)
      {
        if(adjMat[c + r * nnodes])
          fprintf(g, "n%d -- n%d\n", r, c);
      }
    }
    fprintf(g, "}\n");
    fclose(g);
    //Kokkos::finalize();
    return 0;
}

