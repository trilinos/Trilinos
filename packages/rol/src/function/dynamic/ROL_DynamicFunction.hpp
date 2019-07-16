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


#pragma once
#ifndef ROL_DYNAMICFUNCTION_HPP
#define ROL_DYNAMICFUNCTION_HPP

#include <initializer_list>

#include "ROL_TimeStamp.hpp"

#include "ROL_PartitionedVector.hpp"
#include "ROL_VectorWorkspace.hpp"

/** @ingroup dynamic_group
    \class ROL::DynamicFunction
    \brief Provides update interface, casting and vector management to
           DynamicConstraint and DynamicObjective.

    Allowed strings to indicate terms with zero derivatives:
    {"uo","un","z","uo_uo","uo_un","uo_z","un_uo","un_un","un_z","z_uo","z_un", "z_z"}
    Note that is "X_Y" is specified, "Y_X" will be added automatically.
 
*/
namespace ROL {
template<typename Real> 
class DynamicFunction {
public:

  using V  = Vector<Real>;
  using PV = PartitionedVector<Real>;
  using TS = TimeStamp<Real>;

  DynamicFunction( std::initializer_list<std::string> zero_deriv_terms={} ) :
    zero_deriv_terms_(zero_deriv_terms) {
    
    if( is_zero_derivative("uo") )    add_terms("uo_uo","uo_un","uo_z");  
    if( is_zero_derivative("un") )    add_terms("un_uo","un_un","un_z");
    if( is_zero_derivative("z")  )    add_terms("z_uo", "z_un", "z_z" );
    if( is_zero_derivative("uo_un") ) add_terms("un_uo");  
    if( is_zero_derivative("uo_z") )  add_terms("z_uo");  
    if( is_zero_derivative("un_z") )  add_terms("z_un");  
    if( is_zero_derivative("un_uo") ) add_terms("uo_un");  
    if( is_zero_derivative("z_uo") )  add_terms("uo_z");  
    if( is_zero_derivative("z_un") )  add_terms("un_z");  

  }

  virtual ~DynamicFunction() {}

  // Update old state
  virtual void update_uo( const V& x, const TS& ts ) {}

  // Update new state
  virtual void update_un( const V& x, const TS& ts ) {}

  // Update control
  virtual void update_z( const V& x, const TS& ts ) {}

  bool is_zero_derivative( const std::string& key ) {
    return std::find( std::begin(zero_deriv_terms_),
                      std::end(zero_deriv_terms_),
                      key ) != std::end(zero_deriv_terms_); 
  }


protected:

  VectorWorkspace<Real>& getVectorWorkspace() const;

  PV& partition( V& x ) const;
  const PV& partition( const V& x ) const;

  V& getNew( V& x ) const;
  const V& getNew( const V& x ) const;

  V& getOld( V& x ) const;
  const V& getOld( const V& x ) const;  

private:

  mutable VectorWorkspace<Real> workspace_;
  std::vector<std::string>      zero_deriv_terms_;

  template<typename First>
  void add_terms( const First& first ) {
    if( !is_zero_derivative(first) ) 
      zero_deriv_terms_.push_back(first);
  }

  template<typename First, typename...Rest>
  void add_terms( const First& first, const Rest&... rest ) {
    if( !is_zero_derivative(first) ) 
      zero_deriv_terms_.push_back(first);
    add_terms(rest...);
  }

};


} // namespace ROL



#include "ROL_DynamicFunctionDef.hpp"

#endif // ROL_DYNAMICFUNCTION_HPP

