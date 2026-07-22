// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_USERINPUTGENERATOR_HPP
#define ROL_USERINPUTGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_BatchManager.hpp"
#include <fstream>
#include <iostream>
#include <string>

namespace ROL {

template<typename Real> 
class UserInputGenerator : public SampleGenerator<Real> {
private:
  int nSamp_;

  void sample(std::string file_pt, std::string file_wt,
              int n, int dim,
              const ROL::Ptr<BatchManager<Real>> &bman) {
    nSamp_ = n;
    // Read in full point data and weight data
    std::fstream input_pt, input_wt;
    input_pt.open(file_pt.c_str(),std::ios::in);
    input_wt.open(file_wt.c_str(),std::ios::in);
    if ( !input_pt.is_open() || !input_wt.is_open() ) {
      if ( !input_pt.is_open() ) {
        if ( bman->batchID() == 0 ) {
          std::cout << "CANNOT OPEN " << file_pt.c_str() << std::endl;
        }
      }
      if ( !input_wt.is_open() ) {
        if ( bman->batchID() == 0 ) {
          std::cout << "CANNOT OPEN " << file_wt.c_str() << std::endl;
        }
      }
    }
    else {
      std::vector<std::vector<Real>> pt(n);
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
      if ( rank < rem ) N++;
      std::vector<std::vector<Real>> my_pt(N);
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

public:
  UserInputGenerator(ROL::ParameterList &parlist,
               const ROL::Ptr<BatchManager<Real>> &bman)
    : SampleGenerator<Real>(bman) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Sample Generator").sublist("User Input");
    if ( list.isParameter("Points File")  &&
         list.isParameter("Weights File") &&
         list.isParameter("Number of Samples") &&
         list.isParameter("Dimension") ) {
      std::string file_pt = list.get("Points File Name","points.txt");
      std::string file_wt = list.get("Weights File Name","weights.txt");
      int n = list.get("Number of Samples",100);
      int dim = list.get("Dimension",4);
      sample(file_pt,file_wt,n,dim,bman);
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::UserInputGenerator): ParameterList does not contain sufficient information.");
    }
  }

  UserInputGenerator(std::string file_pt,
                     std::string file_wt,
                     int n,
                     int dim,
                     const ROL::Ptr<BatchManager<Real>> &bman)
    : SampleGenerator<Real>(bman) {
    sample(file_pt,file_wt,n,dim,bman);
  }

  void refine(void) {}

  int numGlobalSamples(void) const {
    return nSamp_;
  }
};

}

#endif
