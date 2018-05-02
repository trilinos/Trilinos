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
#ifndef ROL_FINITEDIFFERENCE_HPP
#define ROL_FINITEDIFFERENCE_HPP


/** \class ROL::FiniteDifference
    \brief General interface for numerical approximation of Constraint 
           or Objective methods using the finite difference method. Also
           provides a unified interface for validation.

           ----------------------------------------------------------------

           The FiniteDifference::check() methods are defined for two
           situations:

           1) The tested function has the signature
           
              (void) ( Vector<Real>& y, const Vector<Real>& x ); 

              which is specified by the type FiniteDifference::vector_type
              and it will be tested using a reference function that has the 
              signature
 
              (Real) ( const Vector<Real>& x );
 
              which is specified by the type FiniteDifference::scalar_type.
               

           2) The tested function has the signature

              (void) ( Vector& dy_dv, const Vector& v, const Vector& x );

              which is specified by the type FiniteDifference::dderiv_type
              and it will be tested using a reference function of vector_type.
 
           The check methods require pointers to the two functions above,
           an update function, and the necessary input vectors and result vector, 
           although this vector is not used directly, but is cloned, if 
           necessary, for workspace.

           ----------------------------------------------------------------

           Appropriate pointers to the member functions of the classes 
           to be tested can easily be created using std::bind 
           (defined in the <functional> header).


           While currently only check() is implemented, methods could be
           added to provide finite difference approximations to methods
           as well.

*/

#include <functional>

#include "ROL_Vector.hpp"
#include "ROL_VectorWorkspace.hpp"

namespace ROL {



namespace details {

using namespace std;

template<typename Real>
using f_update_t = function<void( const Vector<Real>& )>;

template<typename Real>
using f_scalar_t = function<Real( const Vector<Real>& )>;

template<typename Real>
using f_vector_t = function<void( Vector<Real>&, const Vector<Real>& )>;

template<typename Real>
using f_dderiv_t = function<void( Vector<Real>&, const Vector<Real>&, const Vector<Real>& )>;


template<typename Real>
class FiniteDifference {
public:
 
  using V = ROL::Vector<Real>;

  FiniteDifference( Teuchos::ParameterList& pl,
                    ostream& os = cout );

  FiniteDifference( const int order = 1, 
                    const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                    const int width = 20,
                    const int precision = 11,
                    const bool printToStream = true,
                    ostream& os = cout );

  virtual ~FiniteDifference(){}

  virtual vector<vector<Real>> scalar_check( f_scalar_t<Real> f_ref, 
                                             f_vector_t<Real> f_test, 
                                             f_update_t<Real> f_update,
                                             const V& result, 
                                             const V& input,
                                             const V& direction,
                                             const string& label ) const;
                                        
  virtual vector<vector<Real>> vector_check( f_vector_t<Real> f_ref, 
                                             f_dderiv_t<Real> f_test, 
                                             f_update_t<Real> f_update,
                                             const V& result, 
                                             const V& input,
                                             const V& direction,
                                             const string& label ) const;
private:

  int          order_;         // Finite difference order (1,2,3, or 4)
  int          numSteps_;      // Number of evalutions of different step sizes
  int          width_;         // For print formatting
  int          precision_;     // Number of digits to display
  bool         printToStream_; // False will suppress output
  vector<Real> steps_;         // Set of step sizes of FD approximation
  
  ostream& os_;                // pointer to Output stream
  
  mutable VectorWorkspace<Real> workspace_;

}; // class FiniteDifference

} // namespace details

using details::FiniteDifference;
using details::f_scalar_t;
using details::f_vector_t;
using details::f_dderiv_t;
using details::f_update_t;


} // namespace ROL

#include "ROL_FiniteDifferenceDef.hpp"

#endif // ROL_FINITEDIFFERENCE_HPP

