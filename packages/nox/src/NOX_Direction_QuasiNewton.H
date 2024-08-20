// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_DIRECTION_QUASINEWTON_H
#define NOX_DIRECTION_QUASINEWTON_H

#include "NOX_Common.H"

#ifdef WITH_PRERELEASE

#include "NOX_Direction_Generic.H" // base class
#include "Teuchos_ParameterList.hpp"       // class data element
#include "Teuchos_RCP.hpp" // class data element

// Forward declaration
namespace NOX {
  class Utils;
  class GlobalData;
}

namespace NOX {
namespace Direction {

//! %Quasi-Newton direction
/*!

  We will calculate a quasi-Newton direction for the related minimization problem

  &nbsp;&nbsp;
  \f$
  f(x) = \frac{1}{2} \|F(x)\|^2.
  \f$

  We compute the direction

  &nbsp;&nbsp;
  \f$d_k = -H_k g_k.\f$

  Here \f$H_k\f$ is a limited-memory BFGS (L-BFGS) approximation to
  the <em>inverse</em> Hessian of \f$f\f$ at \f$x_k\f$, and \f$g_k =
  J_k^T F_k\f$ represents the gradient of \f$f\f$ at \f$x_k\f$.

  The \f$H_k\f$ for L-BFGS is only stored implicitly using information
  from the \f$m\f$ most recent changes in \f$x\f$ and changes in
  \f$g\f$. Specifically, we store pairs \f$(s_i, y_i)\f$ for
  \f$i=0,\dots,m-1\f$.  Here \f$(s_0,y_0)\f$ denotes the <em>oldest</em>
  information, i.e.,

  &nbsp;&nbsp;
  \f$
  s_0 = x_{k+1-m} - x_{k-m}
  \quad\mbox{and}\quad
  y_0 = g_{k+1-m} - g_{k-m}.
  \f$

  Likewise, \f$(s_m,y_m)\f$ denotes the <em>newest</em> information, i.e.,

  &nbsp;&nbsp;
  \f$
  s_{m-1} = x_{k+1} - x_{k}
  \quad\mbox{and}\quad
  y_{m-1} = g_{k+1} - g_{k}.
  \f$

  The two-loop recursion formula for computing \f$ d = -H_k g_k\f$ is
  then given as follows.

  <ol>

  <li> \f$d \;\leftarrow\; -g_k\f$

  <li> For \f$ i=m-1,\dots,0 \f$

  &nbsp;&nbsp;
  \f$\alpha_i = \rho_i s_i^T d \quad \mbox{(store $\alpha_i$)}  \f$

  &nbsp;&nbsp;
  \f$d \;\leftarrow\; d - \alpha_i y_i\f$


  <li> \f$ d \;\leftarrow\; \displaystyle\frac{s_{m-1}^T y_{m-1}}{y_{m-1}^T y_{m-1}} d\f$

  <li> For \f$ i=0,\dots,m-1 \f$

  &nbsp;&nbsp;
  \f$\beta = \rho_i y_i^T d\f$

  &nbsp;&nbsp;
  \f$d \;\leftarrow\; d + (\alpha_i - \beta) s_i\f$

  </ol>

  <b>References</b>

  <ul>
  <li> R. H. Byrd, J. Nocedal, and R. B. Schnabel, Representations of
       quasi-Newton matrices and their use in limited memory methods,
       <em>Math. Prog.</em> 63 (1994): 129-156.
  <li> J. Nocedal, Updating quasi-Newton matrices with limited storage,
       <em>Mathematics of Computation</em> 35 (1980): 773-782.
  </ul>


  <b>Parameters</b>

  To use this direction, specify that the "Method" is "QuasiNewton" in the
  "Direction" sublist of the parameters that are passed to the solver
  (see NOX::Direction::Manager for more information on choosing the
  search direction).

  The following options are valid in the "Quasi-Newton" sublist of the
  "Direction" sublist of the parameters that are passed to the
  nonlinear solver:

  - "Memory" - Number of saved updates, defaults to 5.

 */

class QuasiNewton : public Generic {

public:

  //----------------------------------------------------------------------
  // MemoryUnit helper class
  //----------------------------------------------------------------------

  /*!
    \brief Utility class for NOX::Direction::QuasiNewton method

    Stores an (s,y)-pair for the limited-memory BFGS update
  */
  class MemoryUnit {

  public:

    //! Constructor
    MemoryUnit();

    //! Destuctor
    ~MemoryUnit();

    //! Reset s and y with new vectors
    /*! Sets
      \f[
      s = x_{\rm new} - x_{\rm old},
      \quad
      y = g_{\rm new} - g_{\rm old}
      \f]
    */
    void reset(const NOX::Abstract::Vector& newX, const NOX::Abstract::Vector& oldX,
           const NOX::Abstract::Vector& newG, const NOX::Abstract::Vector& oldG);

    //! Get s
    const NOX::Abstract::Vector& s() const;

    //! Get y
    const NOX::Abstract::Vector& y() const;

    //! Get \f$s^T y \f$
    double sdoty() const;

    //! Get \f$y^T y \f$
    double ydoty() const;

    //! Get \f$\rho = 1 / s^T y\f$
    double rho() const;

  private:

    //! s-vector pointer
    Teuchos::RCP<NOX::Abstract::Vector> sPtr;

    //! y-vector pointer
    Teuchos::RCP<NOX::Abstract::Vector> yPtr;

    //! \f$s^Ty\f$
    double sdotyValue;

    //! \f$y^Ty\f$
    double ydotyValue;

    //! \f$\rho = 1 / s^T y\f$
    double rhoValue;

  };

  //----------------------------------------------------------------------
  // Memory helper class
  //----------------------------------------------------------------------

  /*!
    \brief
    Utility class for NOX::Direction::QuasiNewton method
    to manage the information stored in "limited" memory.

    Store up to \f$m\f$ MemoryUnit objects where \f$m\f$ was passed to
    the constructor or to reset().  Every time add() is called, a new
    MemoryUnit is added. If there are already \f$m\f$ MemoryUnit's, the
    oldest is bumped off the list. The zero entry is always the oldest.
  */
  class Memory {

  public:

    //! Constructor. Maximum number of MemoryUnit's is \c m.
    Memory(int m = 0);

    //! Destructor
    ~Memory();

    //! Reset the memory
    void reset(int m = 0);

    //! Add new information to the memory
    void add(const NOX::Abstract::Vector& newX, const NOX::Abstract::Vector& oldX,
         const NOX::Abstract::Vector& newG, const NOX::Abstract::Vector& oldG);

    //! Returns true if the memory is empty
    bool empty() const;

    //! Number of items in memory
    int size() const;

    /*!
      \brief Return the ith MemoryUnit (unchecked access)

      The zero entry is the oldest memory.

      The m-1 entry is the newest entry (where m denotes the memory size).
    */
    const MemoryUnit& operator[](int i) const;

  private:

    //! Index Into memory vector
    std::vector<int> index;

    //! Memory
    std::vector<MemoryUnit> memory;
  };

//----------------------------------------------------------------------
// END HELPER CLASSES
//----------------------------------------------------------------------

public:

  //! Constructor
  QuasiNewton(const Teuchos::RCP<NOX::GlobalData>& gd,
          Teuchos::ParameterList& params);

  //! Destructor
  virtual ~QuasiNewton();

  // derived
  virtual bool reset(const Teuchos::RCP<NOX::GlobalData>& gd,
             Teuchos::ParameterList& params);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::Generic& solver);

  // derived
  virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
               const NOX::Solver::LineSearchBased& solver);

private:

  //! Print an error message and throw an error
  void throwError(const std::string& functionName, const std::string& errorMsg);

private:

  //! Global data pointer.  Keep this so the parameter list remains valid.
  Teuchos::RCP<NOX::GlobalData> globalDataPtr;

  //! Printing Utils
  Teuchos::RCP<NOX::Utils> utils;

  //! "Direction" sublist with parameters for the direction vector
  /*!
    \note This is pointer rather than a reference to allow for the
    reset function. It should not be deleted.
   */
  Teuchos::ParameterList* paramsPtr;

  //! Storage of \f$(s,y)\f$ pairs and related attributes
  Memory memory;

};
} // namespace Direction
} // namespace NOX

#endif
#endif
