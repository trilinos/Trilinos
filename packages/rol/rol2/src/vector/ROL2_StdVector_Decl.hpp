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

#ifndef ROL2_STDVECTOR_DECL_H
#define ROL2_STDVECTOR_DECL_H

/** \class ROL2::StdVector
    \brief Provides the std::vector implementation of the ROL2::Vector interface.
*/

namespace ROL2 {

template<class Real>
class StdVector : public Vector<Real> {
public:

  using size_type = typename std::vector<Real>::size_type;

  StdVector( const Ptr<std::vector<Real>> & std_vec );

  StdVector( int dim, Real val=0.0 );

  StdVector( std::initializer_list<Real> ilist );

  // Data access methods

  Real& operator[] ( int i );

  const Real& operator[] ( int i ) const;

  Ptr<const std::vector<Real>> getVector() const;

  Ptr<std::vector<Real>> getVector();


  // Overridden Methods
 
  void set( const Vector<Real>& x ) override;

  void plus( const Vector<Real>& x ) override;

  void axpy( Real alpha, const Vector<Real>& x ) override;

  void scale( Real alpha ) override;

  virtual Real dot( const Vector<Real> &x ) const override;

  Real norm() const override;

  virtual Ptr<Vector<Real>> clone() const override;

  Ptr<Vector<Real>> basis( int i ) const override;

  int dimension() const override;

  void applyUnary( const Elementwise::UnaryFunction<Real>& f ) override;

  void applyBinary( const Elementwise::BinaryFunction<Real>& f,
                    const Vector<Real>&                   x ) override;

  Real reduce( const Elementwise::ReductionOp<Real> &r ) const override;

  void setScalar( Real C ) override;

  void randomize( Real l = 0.0, Real u = 1.0 ) override;

  virtual void print( std::ostream &outStream ) const override;

  inline static std::vector<Real>& getData( Vector<Real>& x ) { 
    return *(static_cast<StdVector<Real>&>(x).std_vec_);
  }

  inline static const std::vector<Real>& getData( const Vector<Real>& x ) { 
    return *(static_cast<const StdVector<Real>&>(x).std_vec_);
  }

private:

  Ptr<std::vector<Real>> std_vec_;

}; // class StdVector

} // namespace ROL2

#endif // ROL2_STDVECTOR_DECL_HPP
