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

#ifndef ROL2_TYPEU_NONLINEARCG_DECL_H
#define ROL2_TYPEU_NONLINEARCG_DECL_H

/** @ingroup step_group
    \class ROL2::TypeU::NonlinearCG
    \brief Provides the interface to compute optimization steps
           with nonlinear CG.
*/

namespace ROL2 {
namespace TypeU {

template<typename Real>
class NonlinearCG : public DescentDirection<Real> {
public:

  enum class Type : std::int16_t {
    HestenesStiefel = 0,
    FletcherReeves,
    Daniel,
    PolakRibiere,
    FletcherConjDesc,
    LiuStorey,
    DaiYuan,
    HagerZhang,
    OrenLuenberger,
    UserDefined,
    Last
  };

  static EnumMap<Type> type_dict;

  /** \brief Constructor.

      Constructor to build a NonlinearCG object with a user-defined 
      nonlinear CG object.  Algorithmic specifications are passed in through 
      a ROL::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     nlcg       is a user-defined NonlinearCG object
  */
  NonlinearCG(       ParameterList&          parlist,
               const Ptr<NonlinearCG<Real>>& nlcg = nullPtr);

  void compute(       Vector<Real>&    s, 
                      Real&            snorm, 
                      Real&            sdotg, 
                      int&             iter, 
                      int&             flag,
                const Vector<Real>&    x, 
                const Vector<Real>&    g, 
                      Objective<Real>& obj) override;

  void writeName( std::ostream& os ) const override;

private:

  Ptr<NonlinearCG<Real>> nlcg_; ///< NonlinearCG object (used for quasi-Newton)
  ENonlinearCG enlcg_;
  std::string ncgName_;

}; // class NonlinearCG

template<typename Real>
EnumMap<NonlinearCG<Real>::Type> 
NonlinearCG<Real>::type_dict = { "HestenesStiefel",
                                 "FletcherReeves",
                                 "Daniel (uses Hessian)",
                                 "Polak-Ribiere",
                                 "Fletcher Conjugate Descent",
                                 "Liu-Storey",
                                 "Dai-Yuan",
                                 "Hager-Zhang",
                                 "Oren-Luenberger",
                                 "User Defined" };

} // namespace TypeU
} // namespace ROL2

#endif // ROL2_TYPEU_NONLINEARCG_DECL_H






