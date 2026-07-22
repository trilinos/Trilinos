// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_CONSTRAINT_H
#define LOCA_TURNINGPOINT_MINIMALLYAUGMENTED_CONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterfaceMVDX.H" // base class
#include "NOX_Abstract_MultiVector.H"                       // for dense matrix

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
  namespace TurningPoint {
    namespace MinimallyAugmented {
      class AbstractGroup;
    }
  }
  namespace BorderedSolver {
    class AbstractStrategy;
  }
  namespace TimeDependent {
    class AbstractGroup;
  }
}

namespace LOCA {

  namespace TurningPoint {

    namespace MinimallyAugmented {

      /*!
       * \brief Implementation of
       * LOCA::MultiContinuation::ConstraintInterfaceMVDX for computing
       * turning points for the minimally augmented turning point formulation.
       */
      /*!
       * This class implements the turning point constraint equation
       * \f$\sigma(x,p) = 0\f$ for the minimally augmented turning point
       * formulation where \f$\sigma\f$ is defined via
       * \f[
       *     \begin{bmatrix}
       *       J   & a \\
       *       b^T & 0
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
       *       J^T & b \\
       *       a^T & 0
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
       *     \sigma = -w^T J v/n
       * \f]
       * for any vectors \f$a\f$ and \f$b\f$ in \f$\Re^n\f$.  Using these
       * relationships, it is easy to show
       * \f[
       *     \begin{split}
       *        \sigma_x &= -(w^T J v)_x/n = -w^T J_x v/n \\
       *        \sigma_p &= -(w^T J v)_p/n = -w^T J_p v/n
       *     \end{split}
       * \f]
       *
       * The class is intialized via the \c tpParams parameter list argument
       * to the constructor.  The parameters this class recognizes are:
       * <ul>
       * <li> "Symmetric Jacobian"  -- [bool] (default: false) -
       *    Flag indicating whether Jacobian matrix \f$J\f$ is symmetric,
       *    in which case we force \f$a = b\f$ and therefore the second
       *    tranpose solve for \f$w\f$ is unnecessary
       * <li> "Initial Null Vector Compuation" -- [string]
       *      (default: "User Provided") - Method to compute initial \f$a\f$
       *      and \f$b\f$ vectors.  Valid choices are:
       *      <ul>
       *      <li> "User Provided" - Initial vectors are provided in the
       *           parameter list, in which case the following parameters
       *           are relevant:
       *           <ul>
       *           <li> "Initial A Vector"  --
       *                [Teuchos::RCP<NOX::Abstract::Vector>]
       *                (Must be supplied) -
       *                 Vector storing initial value for \f$a\f$ vector
       *           <li> "Initial B Vector"  --
       *                 [Teuchos::RCP<NOX::Abstract::Vector>]
       *                 (Must be supplied for nonsymmetric Jacobians) -
       *                 Vector storing initial value for \f$b\f$ vector
       *           </ul>
       *      <li> "Solve df/dp" - Compute \f$a = J^{-T}df/dp\f$ and
       *           \f$b = J^{-1} df/dp\f$ where \f$p\f$ is the bifurcation
       *           parameter.
       *      <li> "Constant" - Entries of \f$a\f$ and \f$b\f$ are set to 1.0
       *      </ul>
       * <li> "Null Vector Scaling" -- [string] (default:  "Order N") -
       *      Method to scale \f$a\f$ and \f$b\f$.  This determines the norm of
       *      these vectors and the scaling of \f$\sigma\f$.  Valid choices are:
       *      <ul>
       *      <li> "None" -- Use initial scaling
       *      <li> "Order 1" -- Scale to unit norm
       *      <li> "Order N" -- Use vector length scaling
       *      </ul>
       * <li> "Update Null Vectors Every Continuation Step"  --
       *    [bool] (default: true) -
       *    Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors via
       *    \f$a = w\f$ and \f$b = v\f$ every continuation step
       * <li> "Update Null Vectors Every Nonlinear Iteration"  --
       *    [bool] (default: false) -
       *    Flag indicating whether to update \f$a\f$ and \f$b\f$ vectors via
       *    \f$a = w\f$ and \f$b = v\f$ every nonlinear iteration
       * <li> "Multiply Null Vectors by Mass Matrix"  --
       *    [bool] (default: false) -
       *    Flag indicating whether to multiply \f$a\f$ and \f$b\f$ vectors by
       *    the mass matrix \f$M = \partial f/\partial\dot{x}\f$ at the strart
       *    of a turning point calculation, and each time \f$a\f$ and \f$b\f$
       *    are updated.  This can improve the scaling of these vectors, and
       *    may orthogonalize them against structural null spaces (i.e.,
       *    pressure null space for incompressible Navier-Stokes).
       * </ul>
       */
      class Constraint :
    public LOCA::MultiContinuation::ConstraintInterfaceMVDX {

      public:

    //! Constructor
    Constraint(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& tpParams,
     const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g,
     int bif_param);

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
    virtual void setGroup(const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g);

    //! Returns left null vector w
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getLeftNullVec() const;

    //! Returns right null vector v
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getRightNullVec() const;

    //! Returns a vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getAVec() const;

    //! Returns b vector
    virtual Teuchos::RCP<const NOX::Abstract::Vector>
    getBVec() const;

    //! Returns sigma
    virtual double getSigma() const;

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
     * residuals \f$g\f$ if \c isValidG is \c false.  If \c isValidG is
     * \c true, then the \c dgdp contains \f$g\f$ on input.
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
     * successful.  Here we update the \f$a\f$ and \f$b\f$ vectors to
     * \f$w\f$ and \f$v\f$ respectively if requested.
     */
    virtual void
    postProcessContinuationStep(
                  LOCA::Abstract::Iterator::StepStatus stepStatus);

    //@}

      protected:

    //! Enumerated type determining type of scaling
    enum NullVectorScaling {
      NVS_None,      //! No scaling
      NVS_OrderOne,  //! Scale to O(1)
      NVS_OrderN     //! Scale to O(N) when N is the vector length
    };

    //! Scale a & b vectors
    virtual void scaleNullVectors(NOX::Abstract::Vector& a,
                      NOX::Abstract::Vector& b);

    //! Get initial a & b vectors
    virtual void getInitialVectors(NOX::Abstract::Vector& a,
                       NOX::Abstract::Vector& b);

      private:

    //! Prohibit generation and use of operator=()
    Constraint& operator=(const Constraint& source);

      protected:

    //! Pointer LOCA global data object
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Parsed top-level parameters
    Teuchos::RCP<LOCA::Parameter::SublistParser> parsedParams;

    //! Bifurcation parameter list
    Teuchos::RCP<Teuchos::ParameterList> turningPointParams;

    //! Pointer to base group that defines \f$F\f$
    Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup> grpPtr;

    //! Vector for \f$a\f$
    Teuchos::RCP<NOX::Abstract::MultiVector> a_vector;

    //! Vector for \f$b\f$
    Teuchos::RCP<NOX::Abstract::MultiVector> b_vector;

    //! Stores left null vector
    Teuchos::RCP<NOX::Abstract::MultiVector> w_vector;

    //! Stores right null vector
    Teuchos::RCP<NOX::Abstract::MultiVector> v_vector;

    //! Stores J*v
    Teuchos::RCP<NOX::Abstract::MultiVector> Jv_vector;

    //! Stores J^T*w
    Teuchos::RCP<NOX::Abstract::MultiVector> Jtw_vector;

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

    /*!
     * \brief Flag indicating whether to update \f$a\f$ and \f$b\f$
     * every continuation step.
     */
    bool updateVectorsEveryContinuationStep;

    /*!
     * \brief Flag indicating whether to update \f$a\f$ and \f$b\f$
     * every nonlinear iteration.
     */
    bool updateVectorsEveryIteration;

    //! Null vector scaling method
    NullVectorScaling nullVecScaling;

    //! Multiply null vectors by mass matrix
    bool multiplyMass;

    //! Time dependent group interface for multiplying by mass matrix
    Teuchos::RCP<LOCA::TimeDependent::AbstractGroup> tdGrp;

    //! Temporary vector for mulplying null vectors by mass matrix
    Teuchos::RCP<NOX::Abstract::Vector> tmp_mass;

      }; // Class Constraint

    } // namespace MinimallyAugmented

  } // namespace TurningPoint

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
