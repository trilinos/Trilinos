
//@HEADER
// ************************************************************************
// 
//               MeshingGenie: Fracture Meshing Services Package 
//                 Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

// R5.0

#ifndef HYBRID_MESHER_2D_H
#define HYBRID_MESHER_2D_H

#include "MeshingGenie_2d.h"

class HybridMesher_2d
{
 public:
  //! constructor
  HybridMesher_2d(double dm, std::vector<double> &VoronoiBoundaries,
		  std::vector< std::vector<double> > &Holes,
		  std::vector< std::vector<double> > &Cracks,
		  std::vector<double> &StructuredBoundary, int NumStr,
		  int use_fixed_seed)
    :_dm(dm), _VorBound(VoronoiBoundaries), _Holes(Holes),
    _Cracks(Cracks), _StrBound(StructuredBoundary), _NumStr(NumStr),
    _fixed_seed(use_fixed_seed)
  {};


  //! Destructor
  ~HybridMesher_2d(){ };

  int execute();

  void get_Tessellation(std::vector<double> &x, std::vector<double> &y,
			std::vector< std::vector<size_t> > &elements);

 private:
  void pave_region(size_t h, size_t i, size_t num_voronoi);

  std::vector< std::vector<double> > _Cracks;
  double _dm;
  std::vector< std::vector<size_t> > _elements;
  int _fixed_seed;
  std::vector< std::vector<double> > _Holes;
  int _NumStr;
  std::vector<double> _StrBound;
  std::vector<double> _VorBound;
  std::vector<double> _x;
  std::vector<double> _y;
};

#endif	


