// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MINIMALLYAUGMENTED_CONSTRAINT_H
#define LOCA_HOPF_MINIMALLYAUGMENTED_CONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "NOX_Abstract_MultiVector.H"                       // for dense matrix

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace Hopf {
    namespace MinimallyAugmented {
      class AbstractGroup;
    }
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MinimallyAugmented {

      /*!
       * \brief Implementation of
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX for computing
       * Hopf bifurcations for the minimally augmented Hopf formulation.
       */
      /*!
       * This class implements the turning point constraint equation
       * \f$\sigma(x,p,\omega) = 0\f$ for the minimally augmented Hopf
       * formulation where \f$\sigma\f$ is defined via
       * \f[
       *     \begin{bmatrix}
       *       J+i\omega M & a \\
       *       b^H         & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       v \\
       *       \sigma_1
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *       0 \\
       *       n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \begin{bmatrix}
       *       J^T-i\omega M^T & b \\
       *       a^H             & 0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *       w \\
       *       \sigma_2
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *       0 \\
       *       n
       *     \end{bmatrix},
       * \f]
       * \f[
       *     \sigma = w^H J+i\omega M v/n
       * \f]
       * for any vectors \f$ a\f$ and \f$ b\f$ in \f$\mathbb{C}^n\f$.  Using
       * these relationships, it is easy to show
       * \f[
       *     \begin{split}
       *        \sigma_x &= (w^H(J+i\omega M)v)_x/n = w^H(J+i\omega M)_x v/n \\
       *        \sigma_p &= (w^H(J+i\omega M)v)_p/n = w^H(J+i\omega M)_p v/n
       *     \end{split}
       * \f]
       *
       * The class is intialized via the \c hpfParams parameter list argument
       * to the constructor.  The parameters this class recognizes are:
       * <ul>
       * <li> "Update Null Vectors Every Continuation Step"  --
       *    [bool] (default: true) -
       *    Flag indicating whether to update \f$ a\f$ and \f$ b\f$ vectors via
       *    \f$ a = w\f$ and \f$ b = v\f$ every continuation step
       * <li> "Update Null Vectors Every Nonlinear Iteration"  --
       *    [bool] (default: false) -
       *    Flag indicating whether to update \f$ a\f$ and \f$ b\f$ vectors via
       *    \f$ a = w\f$ and \f$ b = v\f$ every nonlinear iteration
       * </ul>
       */
      class Constraint :
    public LOCA::MultiContinuation::ConstraintInterfaceMVDX {

      public:

    //! Constructor
    Constraint(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& hpfParams,
     const Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup>& g,
     bool is_symmetric,
     const NOX::Abstract::Vector& a_real,
     const NOX::Abstract::Vector& a_imag,
     const NOX::Abstract::Vector* b_real,
     const NOX::Abstract::Vector* b_imag,
     int bif_param,
     double freq);

    //! Copy constructor
    Constraint(const Constraint& source,
           NOX::CopyType type = NOX::DeepCopy);

    //! Destructor
    virtual ~Constraint();

    //! Set the group pointer
    /*!
     * This method should be called when ever the constrained group
     * is copied, since we don't explicitly copy the underlying group
     * here.
     */
    virtual void setGroup(const Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup>& g);

    //! Set Hopf frequency
    virtual void setFrequency(double freq);

    //! Returns real component of left null vector w
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getLeftNullVecReal() const;

    //! Returns imaginary component of left null vector w
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getLeftNullVecImag() const;

    //! Returns real component of right null vector v
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getRightNullVecReal() const;

    //! Returns imaginary component of right null vector v
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getRightNullVecImag() const;

    //! Returns real component of sigma
    virtual double getSigmaReal() const;

    //! Returns imaginary component of sigma
    virtual double getSigmaImag() const;

    /*!
     * @name Implementation of LOCA::MultiContinuation::ConstraintInterface
     * virtual methods
     */
    //@{

    //! Copy
    virtual void
    copy(const LOCA::MultiContinuation::ConstraintInterface& source);

    //! Cloning function
    virtual
    Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
    clone(NOX::CopyType type = NOX::DeepCopy) const;

    //! Return number of constraints
    virtual int numConstraints() const;

    //! Set the solution vector to y.
    virtual void setX(const NOX::Abstract::Vector& y);

    //! Sets parameter indexed by paramID
    virtual void setParam(int paramID, double val);

    //! Sets parameters indexed by paramIDs
    virtual void setParams(
              const std::vector<int>& paramIDs,
              const NOX::Abstract::MultiVector::DenseMatrix& vals);

    //! Compute continuation constraint equations
    virtual NOX::Abstract::Group::ReturnType
    computeConstraints();

    //! Compute derivative of constraints w.r.t. solution vector x
    virtual NOX::Abstract::Group::ReturnType
    computeDX();

    //! Compute derivative of constraints w.r.t. supplied parameters.
    /*!
     * The first column of \c dgdp should be filled with the constraint
     * residuals \f$ g\f$ if \c isValidG is \c false.  If \c isValidG is
     * \c true, then the \c dgdp contains \f$ g\f$ on input.
     */
    virtual NOX::Abstract::Group::ReturnType
    computeDP(const std::vector<int>& paramIDs,
          NOX::Abstract::MultiVector::DenseMatrix& dgdp,
          bool isValidG);

    //! Return \c true if constraint residuals are valid
    virtual bool isConstraints() const;

    //! Return \c true if derivatives of constraints w.r.t. x are valid
    virtual bool isDX() const;

    //! Return constraint residuals
    virtual const NOX::Abstract::MultiVector::DenseMatrix&
    getConstraints() const;

    //! Return solution component of constraint derivatives
    virtual const NOX::Abstract::MultiVector*
    getDX() const;

    /*!
     * \brief Return \c true if solution component of constraint
     * derivatives is zero
     */
    virtual bool isDXZero() const;

    //! Perform any postprocessing after a continuation step finishes.
    /*!
     * The \c stepStatus argument indicates whether the step was
     * successful.  Here we update the \f$ a\f$ and \f$ b\f$ vectors to
     * \f$ w\f$ and \f$ v\f$ respectively if requested.
     */
    virtual void
    postProcessContinuationStep(
                  LOCA::Abstract::Iterator::StepStatus stepStatus);

    //@}

    //! Compute derivative of sigma w.r.t. frequency omega
    virtual NOX::Abstract::Group::ReturnType
    computeDOmega(NOX::Abstract::MultiVector::DenseMatrix& domega);

      private:

    //! Prohibit generation and use of operator=()
    Constraint& operator=(const Constraint& source);

      protected:

    //! Pointer LOCA global data object
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Parsed top-level parameters
    Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

    //! Bifurcation parameter list
    Teuchos::RCP<Teuchos::ParameterList> hopfParams;

    //! Pointer to base group that defines \f$ F\f$
    Teuchos::RCP<LOCA::Hopf::MinimallyAugmented::AbstractGroup> grpPtr;

    //! Vector for \f$ a\f$
    Teuchos::RCP<NOX::Abstract::MultiVector> a_vector;

    //! Vector for \f$ b\f$
    Teuchos::RCP<NOX::Abstract::MultiVector> b_vector;

    //! Stores left null vector
    Teuchos::RCP<NOX::Abstract::MultiVector> w_vector;

    //! Stores right null vector
    Teuchos::RCP<NOX::Abstract::MultiVector> v_vector;

    //! Stores C*v
    Teuchos::RCP<NOX::Abstract::MultiVector> Cv_vector;

    //! Stores sigma_x
    Teuchos::RCP<NOX::Abstract::MultiVector> sigma_x;

    //! Constraint values
    NOX::Abstract::MultiVector::DenseMatrix constraints;

    //! Stores bordered solver strategy
    Teuchos::RCP<LOCA::BorderedSolver::AbstractStrategy> borderedSolver;

    //! Stores vector length as a double
    double dn;

    //! Stores scale factor on sigma
    double sigma_scale;

    /*!
     * \brief Flag indicating whether Jacobian is symmetric
     */
    bool isSymmetric;

    //! Flag indicating whether constraints are valid
    bool isValidConstraints;

    //! Flag indicating whether sigma_x is valid
    bool isValidDX;

    //! Stores the bifurcation parameter index
    std::vector<int> bifParamID;

    //! Stores Hopf frequency
    double omega;

    /*!
     * \brief Flag indicating whether to update \f$ a\f$ and \f$ b\f$
     * every continuation step.
     */
    bool updateVectorsEveryContinuationStep;

    /*!
     * \brief Flag indicating whether to update \f$ a\f$ and \f$ b\f$
     * every nonlinear iteration.
     */
    bool updateVectorsEveryIteration;

      }; // Class Constraint

    } // namespace MinimallyAugmented

  } // namespace Hopf

} // namespace LOCA

#endif // LOCA_HOPF_MINIMALLYAUGMENTED_CONSTRAINT_H
