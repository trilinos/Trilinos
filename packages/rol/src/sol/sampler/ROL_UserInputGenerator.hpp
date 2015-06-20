// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_USERINPUTGENERATOR_HPP
#define ROL_USERINPUTGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"
#include "ROL_BatchManager.hpp"
#include <fstream>
#include <iostream>
#include <string>

namespace ROL {

template<class Real> 
class UserInputGenerator : public SampleGenerator<Real> {
public:
  UserInputGenerator(std::string file_pt, std::string file_wt, int n, int dim, 
    Teuchos::RCP<BatchManager<Real> > &bman) : SampleGenerator<Real>(bman) {
    // Read in full point data and weight data
    std::fstream input_pt;
    input_pt.open(file_pt.c_str(),std::ios::in);
    std::fstream input_wt;
    input_wt.open(file_wt.c_str(),std::ios::in);
    if ( !input_pt.is_open() || !input_wt.is_open() ) {
      if ( !input_pt.is_open() ) {
        if ( bman->batchID() == 0 ) {
          std::cout << "CANNOT OPEN " << file_pt.c_str() << "\n";
        }
      }
      if ( !input_wt.is_open() ) {
        if ( bman->batchID() == 0 ) {
          std::cout << "CANNOT OPEN " << file_wt.c_str() << "\n";
        }
      }
    }
    else {
      std::vector<std::vector<Real> > pt(n);
      std::vector<Real> wt(n,0.0);
      std::vector<Real> point(dim,0.0);;
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < dim; j++) {
          input_pt >> point[j];
        } 
        pt[i] = point;
        input_wt >> wt[i];
      }
      // Get process rankd and number of processes
      int rank  = bman->batchID();
      int nProc = bman->numBatches();
      // Separate samples across processes
      int frac = n/nProc;
      int rem  = n%nProc;
      int N    = frac;
      if ( rank < rem ) {
        N++;
      }
      std::vector<std::vector<Real> > my_pt(N);
      std::vector<Real> my_wt(N,0.0);
      int index = 0;
      for (int i = 0; i < N; i++) {
        index = i*nProc + rank;
        my_pt[i] = pt[index];  
        my_wt[i] = wt[index];
      }
      SampleGenerator<Real>::setPoints(my_pt);
      SampleGenerator<Real>::setWeights(my_wt);
    }
    input_pt.close();
    input_wt.close();
  }

  void refine(void) {}
};

}

#endif
