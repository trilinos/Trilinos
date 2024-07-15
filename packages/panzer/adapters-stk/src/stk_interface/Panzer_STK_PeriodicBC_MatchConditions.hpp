// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_PeriodicBC_MatchConditions_hpp__
#define __Panzer_STK_PeriodicBC_MatchConditions_hpp__

#include "Teuchos_Tuple.hpp"

#include <vector>
#include <string>

namespace panzer_stk {

/** Match coordinates that share the same point on a line that is axis aligned
  */
class CoordMatcher {
   double error_;
   int index_;
   bool relative_; // compute error relative to length of domain
   char labels_[3];

   void buildLabels()
   { labels_[0] = 'x'; labels_[1] = 'y'; labels_[2] = 'z'; }

   void parseParams(const std::vector<std::string> & params)
   {
      std::string errStr = "CoordMatcher \"" + std::string(1,labels_[index_]) + "-coord\" takes at most two parameters <tol, relative>";
      TEUCHOS_TEST_FOR_EXCEPTION(params.size()>2,std::logic_error,errStr);

      // read in string, get double
      if(params.size()>0) {
         std::stringstream ss;
         ss << params[0];
         ss >> error_;
         if(params.size()==2){
           std::string errStr2 = params[1] + " is not a valid periodic option (try \"relative\")";
           TEUCHOS_TEST_FOR_EXCEPTION(params[1]!="relative",std::logic_error,errStr2);
           relative_ = true;
         }
      }
      // else use default value for error
   }

public:
   /// index is the coordinate direction that will be compared to find matching nodes.
   CoordMatcher(int index) : error_(1e-8),index_(index),relative_(false) { buildLabels(); }
   CoordMatcher(int index,double error) : error_(error),index_(index),relative_(false) { buildLabels(); }
   CoordMatcher(int index,const std::vector<std::string> & params) : error_(1e-8),index_(index),relative_(false)
   { buildLabels(); parseParams(params); }

   CoordMatcher(const CoordMatcher & cm) : error_(cm.error_),index_(cm.index_),relative_(cm.relative_) { buildLabels(); }

   bool operator()(const Teuchos::Tuple<double,3> & a,
                   const Teuchos::Tuple<double,3> & b) const
   {
     double error = error_;
     if(relative_) // scale error by length of domain
       error*=std::fabs(a[1-index_]-b[1-index_]);
     return std::fabs(a[index_]-b[index_])<error; /* I'm being lazy here! */
   }

   std::string getString() const
   {
      std::stringstream ss;
      ss << labels_[index_] << "-coord <tol=" << error_ << ">";
      return ss.str();
   }

  int getIndex() const {return index_;}

  int getPeriodicDirection() const
  {
    // This assumes a 2D x-y mesh even though the object supports the
    // z direction. Once in 3D you need the PlaneMatcher.
    TEUCHOS_ASSERT(index_ != 2);
    if (index_ == 0)
      return 1;
    return 0;
  }

  double getAbsoluteTolerance() const {return error_;}

  void transform(double * ptB, const std::vector<double> & centroidA) const
  {
    // Instead of matching directly, shift pt B given the centroid
    // of side A
    // For now, we assume at 2D x-y mesh as above so just need
    // to overwrite the coordinate in the periodic direction

    const int periodicIndex = this->getPeriodicDirection();
    ptB[periodicIndex] = centroidA[periodicIndex];

    return;
  }
};

/** Match coordinates at the same point on a plane
  */
class PlaneMatcher {
   double error_;
   int index0_, index1_;
   bool relative_; // compute error relative to length of domain
   char labels_[3];

   void buildLabels()
   { labels_[0] = 'x'; labels_[1] = 'y'; labels_[2] = 'z'; }

   void parseParams(const std::vector<std::string> & params)
   {
      std::string errStr = "PlaneMatcher \"" + std::string(1,labels_[index0_])+std::string(1,labels_[index1_])
                         + "-coord\" takes at most two parameter <tol, relative>";
      TEUCHOS_TEST_FOR_EXCEPTION(params.size()>2,std::logic_error,errStr);

      // read in string, get double
      if(params.size()>0) {
         std::stringstream ss;
         ss << params[0];
         ss >> error_;
         if(params.size()==2){
           if (params[1] == "3D") { 
            // Warn user but continue
            std::cout << "WARNING : Keyword " << params[1] << " not needed for PlaneMatcher" << std::endl; 
            return;
           }
           std::string errStr2 = params[1] + " is not a valid periodic option (try \"relative\")";
           TEUCHOS_TEST_FOR_EXCEPTION(params[1]!="relative",std::logic_error,errStr2);
           relative_ = true;
         }
      }
      // else use default value for error
      return;
   }

public:
   PlaneMatcher(int index0,int index1) : error_(1e-8),index0_(index0), index1_(index1), relative_(false)
   { TEUCHOS_ASSERT(index0!=index1); buildLabels(); }

   PlaneMatcher(int index0,int index1,double error) : error_(error),index0_(index0), index1_(index1), relative_(false)
   { TEUCHOS_ASSERT(index0!=index1); buildLabels(); }

   PlaneMatcher(int index0,int index1,const std::vector<std::string> & params)
      : error_(1e-8), index0_(index0), index1_(index1), relative_(false)
   { TEUCHOS_ASSERT(index0!=index1); buildLabels(); parseParams(params); }

   PlaneMatcher(const PlaneMatcher & cm) : error_(cm.error_),index0_(cm.index0_), index1_(cm.index1_), relative_(cm.relative_)
   { buildLabels(); }

   bool operator()(const Teuchos::Tuple<double,3> & a,
                   const Teuchos::Tuple<double,3> & b) const
   {
     double error = error_;
     if(relative_) // scale error by length of domain in normal direction
       error*=std::fabs(a[3-index0_-index1_]-b[3-index0_-index1_]);
     return (std::fabs(a[index0_]-b[index0_])<error_)
         && (std::fabs(a[index1_]-b[index1_])<error_) ; /* I'm being lazy here! */
   }

   std::string getString() const
   {
      std::stringstream ss;
      ss << labels_[index0_] << labels_[index1_] << "-coord <tol=" << error_ << ">";
      return ss.str();
   }

  int getIndex0() const {return index0_;}
  int getIndex1() const {return index1_;}
  int getPeriodicDirection() const
  {
    if (index0_ ==0) {
      if (index1_ == 1)
        return 2; // match x,y=periodic in z
      else
        return 1; // match x,z=periodic in y
    }
    else if (index0_ == 1) {
      if (index1_ == 0)
        return 2; // match y,x=periodic in z
      else
        return 0; // match y,z=periodic in x
    }
    else {
      if (index1_ == 0)
        return 1; // match z,x=periodic in y
      else
        return 0; // match z,y=periodic in x
    }
  }
  
  double getAbsoluteTolerance() const {return error_;}

  void transform(double * ptB, const std::vector<double> & centroidA) const
  {
    // Instead of matching directly, shift pt B given the centroid
    // of side A
    // For now, we assume the planes are aligned with one of the
    // coordinate axes so we just need to overwrite the coordinate
    // in the periodic direction

    const int periodicIndex = this->getPeriodicDirection();
    ptB[periodicIndex] = centroidA[periodicIndex];

    return;
  }
};

/** Match coordinates at the same point in two planes. This handles quarter symmetry.
  */
class QuarterPlaneMatcher {
   double error_;
   int index0a_, index0b_, index1_;
   char labels_[3];

   void buildLabels()
   { labels_[0] = 'x'; labels_[1] = 'y'; labels_[2] = 'z'; }

   void parseParams(const std::vector<std::string> & params)
   {
      std::string errStr = "QuarterPlaneMatcher \"(" + std::string(1,labels_[index0a_])+std::string(1,labels_[index0b_])+")"+std::string(1,labels_[index1_])
                         + "-quarter-coord\" takes only one parameter <tol>";
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
   QuarterPlaneMatcher(int index0a,int index0b,int index1)
      : error_(1e-8), index0a_(index0a), index0b_(index0b), index1_(index1)
   { TEUCHOS_ASSERT(index0a!=index1); TEUCHOS_ASSERT(index0b!=index1); buildLabels(); }

   QuarterPlaneMatcher(int index0a,int index0b,int index1,double error)
      : error_(error), index0a_(index0a), index0b_(index0b), index1_(index1)
   { TEUCHOS_ASSERT(index0a!=index1); TEUCHOS_ASSERT(index0b!=index1); buildLabels(); }

   QuarterPlaneMatcher(int index0a,int index0b,int index1,const std::vector<std::string> & params)
      : error_(1e-8), index0a_(index0a), index0b_(index0b), index1_(index1)
   { TEUCHOS_ASSERT(index0a!=index1); TEUCHOS_ASSERT(index0b!=index1); buildLabels(); parseParams(params); }

   QuarterPlaneMatcher(const QuarterPlaneMatcher & cm)
      : error_(cm.error_), index0a_(cm.index0a_), index0b_(cm.index0b_), index1_(cm.index1_)
   { buildLabels(); }

   bool operator()(const Teuchos::Tuple<double,3> & a,
                   const Teuchos::Tuple<double,3> & b) const
   { return (std::fabs(a[index0a_]-b[index0b_])<error_)
         && (std::fabs(a[index1_]-b[index1_])<error_) ; /* I'm being lazy here! */ }

   std::string getString() const
   {
      std::stringstream ss;
      ss << "(" << labels_[index0a_] << labels_[index0b_] << ")" << labels_[index1_] << "-quarter-coord <tol=" << error_ << ">";
      return ss.str();
   }

  double getAbsoluteTolerance() const {return error_;}

  void transform(double * ptB, const std::vector<double> & centroidA) const
  {
    // Instead of matching directly, shift pt B given the centroid
    // of side A
    // For now, we assume the planes are aligned with one of the
    // coordinate axes 
    // We leave ptB[index1_] alone,
    // put ptB[index0b_] in the index0a_ slot and replace it with 
    // centroidA[index0b_] which is the fixed value for plane B

    ptB[index0a_] = ptB[index0b_];
    ptB[index0b_] = centroidA[index0b_];

    return;
  }
};

/** Match coordinates for a 3D wedge. The wedge must be meshed such
    that it is mirrored/split over the xz or yz plane. The index is the
    coordinate index that is compared . If the mirror plane is xz, then
    index0=1. If the mirror plane is yz, then index0=0. A mirror plane
    of xz is specified as "wx" in the string parser. A mirror plane of
    yz is specified as "wy" in the string parser.
  */
class WedgeMatcher {
  double error_;
  /// index to compare - 0 for wy (mirrored over yz), 1 for wx (mirrored over xz)
  int index0_;
  /// Set to true if a 3D problem, set to false if 2D
  bool is_three_d_;
public:
  enum class MirrorPlane : int {
     XZ_PLANE=0,
     YZ_PLANE=1
  };
  WedgeMatcher(MirrorPlane mp,const std::vector<std::string> & params )
    : error_(1e-8),index0_(0),is_three_d_(true)
  {
    if (mp == MirrorPlane::XZ_PLANE)
      index0_ = 1;
    else // YZ_PLANE
      index0_ = 0;

    TEUCHOS_TEST_FOR_EXCEPTION(params.size() > 2,std::logic_error,"WedgeMatcher can only have one or two option parameters (tolerance and dimension)!");

    // read in string, get double
    if (params.size() > 0)
      error_ = std::stod(params[0]);

    if (params.size() > 1) {
      if (params[1] == "2D")
        is_three_d_ = false;
      else if (params[1] == "3D")
        is_three_d_ = true;
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,"ERROR: WedgeMatcher::parsParams() - the second params must be iether \"2D\" or \"3D\", param=" << params[1] << "\n");
      }
    }
  }
  WedgeMatcher(const WedgeMatcher & cm) = default;

  bool operator()(const Teuchos::Tuple<double,3> & a,
                  const Teuchos::Tuple<double,3> & b) const
  {
    if (is_three_d_) {
      return ( (std::fabs(a[index0_]+b[index0_])<error_) &&
               (std::fabs(a[1-index0_]-b[1-index0_])<error_) &&
               (std::fabs(a[2]-b[2])<error_) );
    }

    // else 2D
    return ( (std::fabs(a[index0_]+b[index0_])<error_) &&
             (std::fabs(a[1-index0_]-b[1-index0_])<error_) );
  }

  std::string getString() const
   {
     std::stringstream ss;
     if (index0_ == 0)
       ss << "wy-coord <tol=" << error_ << ">";
     else
       ss << "wx-coord <tol=" << error_ << ">";
     return ss.str();
   }

  int getIndex() const {return index0_;}

  WedgeMatcher::MirrorPlane getMirrorPlane() const
  {
    if (index0_ == 0)
      return MirrorPlane::YZ_PLANE;
    return MirrorPlane::XZ_PLANE;
  }

  bool isThreeD() const {return is_three_d_;}
  
  double getAbsoluteTolerance() const {return error_;}

  void transform(double * ptB, const std::vector<double> & centroidA) const
  {
    // Instead of matching directly, shift pt B given the centroid
    // of side A
    // For now, we assume the wedge is mirrored over the yz or xz plane
    // Then we just need to mirror over the plane

    ptB[index0_] = -ptB[index0_];

    return;
  }
};

} // end panzer_stk

#endif
