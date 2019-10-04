// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_GregoryPatchData_hpp
#define percept_GregoryPatchData_hpp

#include <percept/PerceptMesh.hpp>

#include <cmath>
#include <iostream>

#define MyPow(x,n) std::pow(x,n)
#define MyPow2(x) ((x)*(x))
#define MyInverse(x) (1.0/(x))

namespace percept {

  namespace gregory_patch {

    static int tri_edges[3][5] = {
      {0,1,2,3,4},
      {4,10,14,16,17},
      {17,15,11,5,0}
    };

    static int tri_ribbons[3][2][4][2] = {
      //ribbons: (convention is interior is on the left, 1st 4 are "q" boundary points, 2nd 4 are interior "p" points)
      //  Note: DL(q[0],q[1]) = ql[1] = degree lower of quartic given by q[i],i=0,4 = 1/3 (4 q[1] - q[0]) (Note: ql[0] = q[0])
      //  Note: DH(ql[0],ql[1]) = qh[1] = degree higher of cubic ql[i],i=0,3 = (ql[0] + 3 ql[1]) / 4

      //v = 0 (edge 0)
      //{ {0, DL(0,1), DL(4,3), 4}, {DL(0,5), 6, 8, DL(4,10)} }
      { { {0,0}, {0,1}, {4,3}, {4,4} }, {{0,5}, {6,6}, {8,8}, {4,10} } },

      //w = 0 (edge 1)
      //{ {4, DL(4,10), DL(17,16), 17}, {DL(4,3), 9, 12, DL(17,15) } }
      { {{4,4}, {4,10}, {17,16}, {17,17} }, {{4,3}, {9,9}, {12,12}, {17,15} } },

      //u = 0 (edge 2)
      //{ {17, DL(17,15), DL(0,5), 0}, {DL(17,16), 13, 7, DL(0,1)}}
      { {{17,17}, {17,15}, {0,5}, {0,0} }, {{17,16}, {13,13}, {7,7}, {0,1} } }
    };

    static int quad_edges[4][4] = {
      //v = 0
      {0,1,2,3},
      //u = 1
      {3,9,15,19},
      //v = 1
      {19,18,17,16},
      //u = 0
      {16, 10, 4, 0}
    };

    //ribbons: (convention is interior is on the left, 1st 4 are "q" boundary points, 2nd 4 are interior "p" points)
    static int quad_ribbons[4][2][4] = {
      //v = 0
      {{0,1,2,3}, {4,5,7,9}},
      //u = 1
      {{3,9,15,19}, {2, 8, 14, 18}},
      //v = 1
      {{19,18,17,16}, {15,13,11,10}},
      //u = 0
      {{16, 10, 4, 0}, {17, 12, 6, 1}}
    };

  }

}
#endif
