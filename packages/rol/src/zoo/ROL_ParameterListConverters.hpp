// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PARAMETER_LIST_CONVERTERS_H
#define ROL_PARAMETER_LIST_CONVERTERS_H

#include <map>
#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"


namespace ROL {

namespace StringList {

/// Helper function for making vectors of strings
inline std::vector<std::string> join( const std::string &s1,
                                      const std::string &s2 ) {
  std::vector<std::string> v; 
  v.push_back(s1);
  v.push_back(s2);

  return v;
}

/// Helper function for making vectors of strings
inline std::vector<std::string> join( const std::string &s1,
                                      const std::string &s2,
                                      const std::string &s3 ) {
  std::vector<std::string> v; 
  v.push_back(s1);
  v.push_back(s2);
  v.push_back(s3);

  return v;
}

/// Helper function for making vectors of strings
inline std::vector<std::string> join( const std::string &s1,
                                      const std::string &s2,
                                      const std::string &s3,
                                      const std::string &s4 ) {
  std::vector<std::string> v; 
  v.push_back(s1);
  v.push_back(s2);
  v.push_back(s3);
  v.push_back(s4);

  return v;
}

/// Helper function for making vectors of strings
inline std::vector<std::string> join( const std::string &s1,
                                      const std::string &s2,
                                      const std::string &s3,
                                      const std::string &s4,
                                      const std::string &s5 ) {
  std::vector<std::string> v; 
  v.push_back(s1);
  v.push_back(s2);
  v.push_back(s3);
  v.push_back(s4);
  v.push_back(s5);

  return v;
}

} // namespace StringList

template<class ParameterType>
void setParameter( ROL::ParameterList &parlist, 
                   const std::vector<std::string> &location,
                   const std::vector<std::string>::iterator iter,
                   ParameterType value ) {

  if( iter == location.end()-1 ) { // Key name 
    parlist.set(*iter,value); 
  }
  else { // sublist 
    ROL::ParameterList &sublist = parlist.sublist(*iter);
    setParameter(sublist,location,iter+1,value);
  }
 

}



/// Produce a heirarchical parameter list using the new names from a flat list of the old names
inline void tierParameterList( ROL::ParameterList &outList, 
                               const ROL::ParameterList &inList ) {

  using StringList::join;

  typedef std::string                  Str;
  typedef std::vector<Str>             Vec;
  typedef std::map<Str,Vec>            Map;
  typedef ParameterList::ConstIterator IterPL;
  typedef typename Vec::iterator       IterVec;
  typedef typename Map::iterator       IterMap;
  
  Map dict;

  // Original flat list name                                  heirarchical list name 
  dict["Use Inexact Gradient"]                              = join("General","Inexact Gradient");
  dict["Use Inexact Objective Function"]                    = join("General","Inexact Objective Function");
  dict["Use Inexact Hessian-Times-A-Vector"]                = join("General","Inexact Hessian-Times-A-Vector");
  dict["Use Projected Gradient Criticality Measure"]        = join("General","Projected Gradient Criticality Measure");
  dict["Scale for Epsilon Active Sets"]                     = join("General","Scale for Epsilon Active Sets");

  dict["Absolute Krylov Tolerance"]                         = join("General","Krylov","Absolute Tolerance");
  dict["Relative Krylov Tolerance"]                         = join("General","Krylov","Relative Tolerance");
  dict["Maximum Number of Krylov Iterations"]               = join("General","Krylov","Iteration Limit");
  dict["Krylov Type"]                                       = join("General","Krylov","Type");

  dict["Barzilai-Borwein"]                                  = join("General","Secant","Barzilai-Borwein");
  dict["Maximum Secant Storage"]                            = join("General","Secant","Maximum Storage");
  dict["Secant Type"]                                       = join("General","Secant","Type");
  dict["Use Secant Hessian-Times-A-Vector"]                 = join("General","Secant","Use as Hessian");
  dict["Use Secant Preconditioning"]                        = join("General","Secant","Use as Preconditioner");

  dict["Gradient Tolerance"]                                = join("Status Test","Gradient Tolerance");
  dict["Maximum Number of Iterations"]                      = join("Status Test","Iteration Limit");
  dict["Step Tolerance"]                                    = join("Status Test","Step Tolerance");

  dict["Accept Last Alpha"]                                 = join("Step","Line Search","Accept Last Alpha");
  dict["Accept Linesearch Minimizer"]                       = join("Step","Line Search","Accept Linesearch Minimizer");
  dict["Maximum Number of Function Evaluations"]            = join("Step","Line Search","Function Evaluation Limit");
  dict["Initial Linesearch Parameter"]                      = join("Step","Line Search","Initial Step Size");
  dict["Sufficient Decrease Parameter"]                     = join("Step","Line Search","Sufficient Decrease Tolerance");
  dict["User Defined Linesearch Parameter"]                 = join("Step","Line Search","User Defined Initial Step Size");

  dict["Curvature Conditions Parameter"]                    = join("Step","Line Search","Curvature Condition","General Parameter");
  dict["Curvature Conditions Parameter: Generalized Wolfe"] = join("Step","Line Search","Curvature Condition","Generalized Wolfe Parameter");
  dict["Linesearch Curvature Condition"]                    = join("Step","Line Search","Curvature Condition","Type");

  dict["Nonlinear CG Type"]                                 = join("Step","Line Search","Descent Method","Nonlinear CG Type");
  dict["Descent Type"]                                      = join("Step","Line Search","Descent Method","Type");

  dict["Backtracking Rate"]                                 = join("Step","Line Search","Line-Search Method","Backtracking Rate");
  dict["Bracketing Tolerance"]                              = join("Step","Line Search","Line-Search Method","Bracketing Tolerance");
  dict["Linesearch Type"]                                   = join("Step","Line Search","Line-Search Method","Type");

  dict["Initial Trust-Region Radius"]                       = join("Step","Trust Region","Initial Radius");
  dict["Maximum Trust-Region Radius"]                       = join("Step","Trust Region","Maximum Radius");
  dict["Radius Growing Threshold"]                          = join("Step","Trust Region","Radius Growing Threshold");
  dict["Radius Growing Rate"]                               = join("Step","Trust Region","Radius Growing Rate");
  dict["Radius Shrinking Threshold"]                        = join("Step","Trust Region","Radius Shrinking Threshold");
  dict["Trust-Region Safeguard"]                            = join("Step","Trust Region","Safeguard Size");
  dict["Trust-Region Subproblem Solver Type"]               = join("Step","Trust Region","Subproblem Solver");
  dict["Step Acceptance Parameter"]                         = join("Step","Trust Region","Step Acceptance Threshold");

  dict["Gradient Update Relative Tolerance"]                = join("Step","Trust Region","Gradient","Relative Tolerance");
  dict["Gradient Update Tolerance Scaling"]                 = join("Step","Trust Region","Gradient","Tolerance Scaling");
  dict["Value Update Exponent"]                             = join("Step","Trust Region","Inexact","Value","Exponent");
  dict["Value Update Forcing Sequence Initial Value"]       = join("Step","Trust Region","Inexact","Value","Forcing Sequence Initial Value");
  dict["Value Update Forcing Sequence Reduction Factor"]    = join("Step","Trust Region","Inexact","Value","Forcing Sequence Reduction Factor");
  dict["Value Update Forcing Sequence Update Frequency"]    = join("Step","Trust Region","Inexact","Value","Forcing Sequence Update Frequency");
  dict["Value Update Tolerance Scaling"]                    = join("Step","Trust Region","Inexact","Value","Tolerance Scaling");


  // Add duplicate entries with unformatted keys
  for(IterMap itmap = dict.begin(); itmap != dict.end(); ++itmap) {
    Str key = itmap->first;
    Vec value = itmap->second;
    dict[removeStringFormat(key)] = value;
  }
    

  // Iterate over parameter list
  for(IterPL itpl = inList.begin(); itpl != inList.end(); ++itpl) {

    // Get name of old key
    Str key( inList.name(itpl) ); 

    // Look up location/name of new key
    Vec location = dict[removeStringFormat(key)];
    
    // Skip if not found in map
    if(location.size() != 0) {

      IterVec itvec = location.begin();

      if( inList.isType<bool>(key) ) {
        bool value = inList.get<bool>( key );
        setParameter( outList, location, itvec, value );
      }
      else if( inList.isType<int>(key) ) {
        int value = inList.get<int>( key );
        setParameter( outList, location, itvec, value );
      }
      else if( inList.isType<double>(key) ) {
        double value = inList.get<double>( key );
        setParameter( outList, location, itvec, value );
      }
      else if( inList.isType<std::string>(key) ) {
        std::string value = inList.get<std::string>( key );
        setParameter( outList, location, itvec, value );
      }
      else {
        ROL_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>>Error: Unhandled parameter type." );  
      }
    } 
   
  }
}

} // namespace ROL

#endif 
