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

#ifndef ROL_TYPEU_LINESEARCHALGORITHM_H
#define ROL_TYPEU_LINESEARCHALGORITHM_H

#include "ROL_TypeU_Algorithm.hpp"
#include "ROL_LineSearch_U_Types.hpp"
#include "ROL_DescentDirection_U.hpp"
#include "ROL_LineSearch_U.hpp"

/** \class ROL::TypeU::LineSearchAlgorithm
    \brief Provides an interface to run unconstrained line search algorithms.

    Suppose \f$\mathcal{X}\f$ is a Hilbert space of 
    functions mapping \f$\Xi\f$ to \f$\mathbb{R}\f$.  For example, 
    \f$\Xi\subset\mathbb{R}^n\f$ and \f$\mathcal{X}=L^2(\Xi)\f$ or 
    \f$\Xi = \{1,\ldots,n\}\f$ and \f$\mathcal{X}=\mathbb{R}^n\f$. We 
    assume \f$f:\mathcal{X}\to\mathbb{R}\f$ is twice-continuously Fr&eacute;chet 
    differentiable and \f$a,\,b\in\mathcal{X}\f$ with \f$a\le b\f$ almost 
    everywhere in \f$\Xi\f$.  Note that these line-search algorithms will also work 
    with secant approximations of the Hessian. 
    This step applies to unconstrained optimization problems,
    \f[
        \min_x\quad f(x).
    \f]

    Given the \f$k\f$-th iterate \f$x_k\f$ and a descent direction
    \f$s_k\f$, the line search approximately minimizes the 1D objective function 
    \f$\phi_k(t) = f(x_k + t s_k)\f$.  The approximate minimizer \f$t\f$ must satisfy 
    sufficient decrease and curvature conditions into guarantee global convergence.  The 
    sufficient decrease condition (often called the Armijo condition) is 
    \f[
       \phi_k(t) \le \phi_k(0) + c_1 t \phi_k'(0) \qquad\iff\qquad
       f(x_k+ts_k) \le f(x_k) + c_1 t \langle \nabla f(x_k),s_k\rangle_{\mathcal{X}}
    \f]
    where \f$0 < c_1 < 1\f$.  The curvature conditions implemented in ROL include:

    <CENTER>
    | Name              | Condition                                                     | Parameters |
    | :---------------- | :-----------------------------------------------------------: | :---------------------: |
    | Wolfe             | \f$\phi_k'(t) \ge c_2\phi_k'(0)\f$                            | \f$c_1<c_2<1\f$         |
    | Strong Wolfe      | \f$\left|\phi_k'(t)\right| \le c_2 \left|\phi_k'(0)\right|\f$ | \f$c_1<c_2<1\f$         |
    | Generalized Wolfe | \f$c_2\phi_k'(0)\le \phi_k'(t)\le-c_3\phi_k'(0)\f$            | \f$0<c_3<1\f$           |
    | Approximate Wolfe | \f$c_2\phi_k'(0)\le \phi_k'(t)\le(2c_1-1)\phi_k'(0)\f$        | \f$c_1<c_2<1\f$         |
    | Goldstein         | \f$\phi_k(0)+(1-c_1)t\phi_k'(0)\le \phi_k(t)\f$               | \f$0<c_1<\frac{1}{2}\f$ |
    </CENTER>
    
    Note that \f$\phi_k'(t) = \langle \nabla f(x_k+ts_k),s_k\rangle_{\mathcal{X}}\f$.

    LineSearchStep implements a number of algorithms for unconstrained 
    optimization.  These algorithms are: Steepest descent; Nonlinear CG; Quasi-Newton methods;
    Inexact Newton methods; Newton's method. These methods are chosen through the EDescent enum.
*/

namespace ROL {
namespace TypeU {

template <class Real>
class LineSearchAlgorithm : public Algorithm<Real> {
private:

  Ptr<DescentDirection_U<Real>> desc_;       ///< Unglobalized step object
  Ptr<LineSearch_U<Real>>       lineSearch_; ///< Line-search object

  EDescentU            edesc_; ///< enum determines type of descent direction
  ELineSearchU         els_;   ///< enum determines type of line search
  ECurvatureConditionU econd_; ///< enum determines type of curvature condition

  bool acceptLastAlpha_;  ///< For backwards compatibility. When max function evaluations are reached take last step

  bool usePreviousAlpha_; ///< If true, use the previously accepted step length (if any) as the new initial step length

  int verbosity_;
  bool printHeader_;
  int ls_nfval_, ls_ngrad_;
  int SPflag_, SPiter_;

  std::string lineSearchName_, descentName_;

  using Algorithm<Real>::state_;
  using Algorithm<Real>::status_;
  using Algorithm<Real>::initialize;

public:
  /** \brief Constructor.

      Standard constructor to build a LineSearchStep object.  Algorithmic 
      specifications are passed in through a ROL::ParameterList.  The
      line-search type, secant type, Krylov type, or nonlinear CG type can
      be set using user-defined objects.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     lineSearch is a user-defined descent direction object
      @param[in]     lineSearch is a user-defined line search object
  */
  LineSearchAlgorithm( ParameterList &parlist,
                       const Ptr<DescentDirection_U<Real>> &descent = nullPtr,
                       const Ptr<LineSearch_U<Real>> &lineSearch = nullPtr );

  void initialize(const Vector<Real> &x, const Vector<Real> &g,
                  Objective<Real> &obj, std::ostream &outStream = std::cout);

  void run( Vector<Real>       &x,
            const Vector<Real> &g, 
            Objective<Real>    &obj,
            std::ostream       &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;
  
  void writeOutput( std::ostream& os, bool print_header = false ) const override;

}; // class ROL::TypeU::LineSearchAlgorithm

} // namespace TypeU
} // namespace ROL

#include "ROL_TypeU_LineSearchAlgorithm_Def.hpp"

#endif
