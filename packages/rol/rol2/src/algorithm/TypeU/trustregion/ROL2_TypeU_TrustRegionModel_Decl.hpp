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

#ifndef ROL2_TYPEU_TRUSTREGIONMODEL_DECL_H
#define ROL2_TYPEU_TRUSTREGIONMODEL_DECL_H

/** @ingroup func_group
    \class ROL2::TypeU::TrustRegionModel
    \brief Provides the interface to evaluate trust-region model functions.

    ROL2::TypeU::TrustRegionModel provides the interface to implement a number of
    trust-region models for unconstrained and constrained optimization.
    The default implementation is the standard quadratic trust region model
    for unconstrained optimization.

    -----
*/

namespace ROL2 {
namespace TypeU {

template<typename Real>
class TrustRegionModel : public Objective<Real> {
public:

  enum class Type : std::int16_t {
    ColemanLi = 0,
    KelleySachs,
    LineMore,
    Last
  };

  virtual ~TrustRegionModel() = default;

  TrustRegionModel( ParameterList&           list,
                    const Ptr<Secant<Real>>& secant = nullPtr,
                    Secant<Real>::Mode       mode   = Secant<Real>::Mode::Both );

  void initialize( const Vector<Real>& x, 
                   const Vector<Real>& g );

  // Some versions of Clang will issue a warning that update hides and 
  // overloaded virtual function without this using declaration
  using Objective<Real>::update;

  void validate(       Objective<Real>&        obj,
                 const Vector<Real>&           x,
                 const Vector<Real>&           g,
                       TrustRegion<Real>::Type etr );
        
  virtual void setData(      Objective<Real>& obj,
                       const Vector<Real>&    x,
                       const Vector<Real>&    g );

  void update( const Vector<Real>& x, 
               const Vector<Real>& s,
               const Vector<Real>& gold, 
               const Vector<Real>& gnew,
                     Real          snorm, 
                     int           iter );

  //-----------------------------------------------------------
  // Override Objective Interface

  virtual Real value( const Vector<Real>& s, 
                            Real&         tol ) override;

  virtual void gradient(       Vector<Real>& g, 
                         const Vector<Real>& s, 
                               Real &tol ) override;

  virtual void hessVec(       Vector<Real>& hv, 
                        const Vector<Real>& v, 
                        const Vector<Real>& s, 
                              Real&         tol ) override {
    applyHessian(hv,v,tol);
  }

  //-----------------------------------------------------------

  virtual void invHessVec(       Vector<Real>& hv, 
                           const Vector<Real>& v, 
                           const Vector<Real>& s, 
                                 Real&         tol ) override {
    applyInvHessian(hv,v,tol);
  }

  virtual void precond(       Vector<Real>& Pv, 
                        const Vector<Real>& v, 
                        const Vector<Real>& s, 
                              Real&         tol ) override {
    applyPrecond(Pv,v,tol);
  }

  void setGradient(  const Ptr<Vector<Real>>&    g )      { g_ = g;           }
  void setIterate(   const Ptr<Vector<Real>>&    x )      { x_ = x;           }
  void setObjective( const Ptr<Objective<Real>>& obj )    { obj_ = obj;       }
  void setSecant(    const Ptr<Secant<Real>>&    secant ) { secant_ = secant; }

  const Vector<Real>&    getGradient()  const { return g_;      }
  const Vector<Real>&    getIterate()   const { return x_;      }
  const Objective<Real>& getObjective() const { return obj_;    }
  const Secant<Real>&    getSecant()    const { return secant_; }

protected:

  Vector<Real>&    getGradient()  { return g_;      }
  Vector<Real>&    getIterate()   { return x_;      }
  Objective<Real>& getObjective() { return obj_;    }
  Secant<Real>&    getSecant()    { return secant_; }
 
  void applyHessian(       Vector<Real>& hv, 
                     const Vector<Real>& v, 
                           Real&         tol );

  void applyInvHessian(       Vector<Real>& hv, 
                        const Vector<Real>& v, 
                              Real&         tol );

  void applyPrecond(       Vector<Real>& Pv, 
                     const Vector<Real>& v, 
                           Real&         tol );

  static EnumMap<Type> type_dict;

private:

  Ptr<Objective<Real>> obj_;
  Ptr<Vector<Real>>    x_, g_, dual_;
  Ptr<Secant<Real>>    secant_;

  bool useSecantPrecond_;
  bool useSecantHessVec_;

}; // class TrustRegionModel

template<class Real>
EnumMap<TrustRegionModel<Real>::Type> 
TrustRegionModel<Real>::type_dict = { "Coleman-Li", 
                                      "Kelley-Sachs",
                                      "Lin-More" };


} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_TRUSTREGIONMODEL_DECL_H
