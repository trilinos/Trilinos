// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_STK_PeriodicBC_MatchConditions_hpp__
#define __Panzer_STK_PeriodicBC_MatchConditions_hpp__

#include "Teuchos_Tuple.hpp"

#include <vector>
#include <string>

namespace panzer_stk {

/** Match coordinates that share the same point on a line
  */ 
class CoordMatcher {
   double error_;
   int index_;
   char labels_[3];

   void buildLabels()
   { labels_[0] = 'x'; labels_[1] = 'y'; labels_[2] = 'z'; }

   void parseParams(const std::vector<std::string> & params) 
   { 
      std::string errStr = "CoordMatcher \"" + std::string(1,labels_[index_]) + "-coord\" takes at most one parameter <tol>";
      TEUCHOS_TEST_FOR_EXCEPTION(params.size()>1,std::logic_error,errStr);
 
      // read in string, get double
      if(params.size()==1) {
         std::stringstream ss;
         ss << params[0];
         ss >> error_; 
      }
      // else use default value for error
   }

public:
   CoordMatcher(int index) : error_(1e-8),index_(index) { buildLabels(); }
   CoordMatcher(int index,double error) : error_(error),index_(index) { buildLabels(); }
   CoordMatcher(int index,const std::vector<std::string> & params) : error_(1e-8),index_(index) 
   { buildLabels(); parseParams(params); }

   CoordMatcher(const CoordMatcher & cm) : error_(cm.error_),index_(cm.index_) { buildLabels(); }

   bool operator()(const Teuchos::Tuple<double,3> & a,
                   const Teuchos::Tuple<double,3> & b) const
   { return std::fabs(a[index_]-b[index_])<error_; /* I'm being lazy here! */ }

   std::string getString() const 
   { 
      std::stringstream ss;
      ss << labels_[index_] << "-coord <tol=" << error_ << ">";
      return ss.str();
   }
};

/** Match coordinates at the same point on a plane
  */ 
class PlaneMatcher {
   double error_;
   int index0_, index1_;
   char labels_[3];
  
   void buildLabels()
   { labels_[0] = 'x'; labels_[1] = 'y'; labels_[2] = 'z'; }

   void parseParams(const std::vector<std::string> & params) 
   { 
      std::string errStr = "PlaneMatcher \"" + std::string(1,labels_[index0_])+std::string(1,labels_[index1_]) 
                         + "-coord\" takes only one parameter <tol>";
      TEUCHOS_TEST_FOR_EXCEPTION(params.size()>1,std::logic_error,errStr);
 
      // read in string, get double
      if(params.size()==1) {
         std::stringstream ss;
         ss << params[0];
         ss >> error_; 
      }
      // else use default value for error
   }

public:
   PlaneMatcher(int index0,int index1) : error_(1e-8),index0_(index0), index1_(index1) 
   { TEUCHOS_ASSERT(index0!=index1); buildLabels(); }

   PlaneMatcher(int index0,int index1,double error) : error_(error),index0_(index0), index1_(index1) 
   { TEUCHOS_ASSERT(index0!=index1); buildLabels(); }

   PlaneMatcher(int index0,int index1,const std::vector<std::string> & params) 
      : error_(1e-8), index0_(index0), index1_(index1) 
   { TEUCHOS_ASSERT(index0!=index1); buildLabels(); parseParams(params); }

   PlaneMatcher(const PlaneMatcher & cm) : error_(cm.error_),index0_(cm.index0_), index1_(cm.index1_) 
   { buildLabels(); }

   bool operator()(const Teuchos::Tuple<double,3> & a,
                   const Teuchos::Tuple<double,3> & b) const
   { return (std::fabs(a[index0_]-b[index0_])<error_) 
         && (std::fabs(a[index1_]-b[index1_])<error_) ; /* I'm being lazy here! */ }

   std::string getString() const 
   { 
      std::stringstream ss;
      ss << labels_[index0_] << labels_[index1_] << "-coord <tol=" << error_ << ">";
      return ss.str();
   }
};

} // end panzer_stk

#endif
