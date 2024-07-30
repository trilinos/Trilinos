// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_PITCHFORK_MOORESPENCE_SALINGERBORDERING_H
#define LOCA_PITCHFORK_MOORESPENCE_SALINGERBORDERING_H

#include "LOCA_Pitchfork_MooreSpence_SolverStrategy.H"
#include "NOX_Abstract_MultiVector.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Pitchfork {

    namespace MooreSpence {

      /*!
       * \brief Moore-Spence pitchfork solver strategy based on
       * "Salinger" bordering.  This is the classic 6-solve bordering
       * method.
       */
      /*!
       * This class solves the Moore-Spence pitchfork Newton equations:
       * \f[
       *     \begin{bmatrix}
       *          J   &    0   & \psi &  f_p   \\
       *       (Jv)_x &    J   &   0  & (Jv)_p \\
       *       \psi^T &    0   &   0  &   0    \\
       *          0   & \phi^T &   0  &   0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y \\
       *        w \\
       *        z
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *        F \\
       *        G \\
       *        s \\
       *        h
       *     \end{bmatrix}
       *  \f]
       *  via the following block elimination scheme:
       *  \f[
       *    \begin{split}
       *       J [A \; b \; c] &= [F \; f_p \; \psi] \\
       *       J [D \; e \; g] &= [G \; (Jv)_p \; 0] - (Jv)_x[A \; b \; c] \\
       *       w &= \frac{\phi^T e(\langle A, \psi\rangle - s) -
                          \langle b, \psi\rangle(\phi^T D - h)}
             {\phi^T e \langle c, \psi\rangle -
              \phi^T g \langle b, \psi\rangle} \\
       *       z &= \frac{\phi^T D - h - (\phi^T g) w}{\phi^T e} \\
       *       X &= A - b z - c w \\
       *       Y &= D - e z - g w
       *    \end{split}
       *   \f]
       */
      class SalingerBordering :
    public LOCA::Pitchfork::MooreSpence::SolverStrategy {

      public:

    //! Constructor.
    /*!
     * \param global_data [in] Global data object
     * \param topParams [in] Parsed top-level parameter list
     * \param solverParams [in] Bordered solver parameters.  Currently none
     * are referenced.
     */
    SalingerBordering(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

    //! Destructor
    virtual ~SalingerBordering();

    //! Set blocks in extended linear system
    /*!
     * \param group [in] Underlying group representing J
     * \param pfGroup [in] Pitchfork group representing the pitchfork
     *        equations.
     * \param asymMultiVector [in] Multivector representing the
     *                             asymmetric vector
     * \param nullVector [in] Vector representing v
     * \param JnVector [in] Vector representing Jv
     * \param dfdp [in] Vector representing df/dp
     * \param dJndp [in] Vector representing d(Jv)/dp
     */
    virtual void setBlocks(
      const Teuchos::RCP<LOCA::Pitchfork::MooreSpence::AbstractGroup>& group,
      const Teuchos::RCP<LOCA::Pitchfork::MooreSpence::ExtendedGroup>& pfGroup,
      const Teuchos::RCP<const NOX::Abstract::MultiVector>& asymMultiVector,
      const Teuchos::RCP<const NOX::Abstract::Vector>& nullVector,
      const Teuchos::RCP<const NOX::Abstract::Vector>& JnVector,
      const Teuchos::RCP<const NOX::Abstract::Vector>& dfdp,
      const Teuchos::RCP<const NOX::Abstract::Vector>& dJndp);

    /*!
     * \brief Solves the extended system as defined above
     */
    /*!
     * The \em params argument is the linear solver parameters.
     */
    virtual NOX::Abstract::Group::ReturnType
    solve(
         Teuchos::ParameterList& params,
             const LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& input,
             LOCA::Pitchfork::MooreSpence::ExtendedMultiVector& result) const;

      protected:

    //! Solves equations with contiguous arguments
    NOX::Abstract::Group::ReturnType
    solveContiguous(
              Teuchos::ParameterList& params,
          const NOX::Abstract::MultiVector& input_x,
          const NOX::Abstract::MultiVector& input_null,
          const NOX::Abstract::MultiVector::DenseMatrix& input_slack,
              const NOX::Abstract::MultiVector::DenseMatrix& input_param,
          NOX::Abstract::MultiVector& result_x,
          NOX::Abstract::MultiVector& result_null,
          NOX::Abstract::MultiVector::DenseMatrix& result_slack,
              NOX::Abstract::MultiVector::DenseMatrix& result_param) const;

      private:

    //! Private to prohibit copying
    SalingerBordering(const SalingerBordering&);

    //! Private to prohibit copying
    SalingerBordering& operator = (const SalingerBordering&);

      protected:

    //! Global data object
    Teuchos::RCP<LOCA::GlobalData> globalData;

    //! Solver parameters
    Teuchos::RCP<Teuchos::ParameterList> solverParams;

    //! Underlying group
    Teuchos::RCP<LOCA::Pitchfork::MooreSpence::AbstractGroup> group;

    //! Turning point group
    Teuchos::RCP<LOCA::Pitchfork::MooreSpence::ExtendedGroup> pfGroup;

    //! Antisymmetric vector (psi)
    Teuchos::RCP<const NOX::Abstract::MultiVector> asymMultiVector;

    //! Antisymmetric vector (psi)
    Teuchos::RCP<const NOX::Abstract::Vector> asymVector;

    //! Null vector
    Teuchos::RCP<const NOX::Abstract::Vector> nullVector;

    //! Jacobian times null vector
    Teuchos::RCP<const NOX::Abstract::Vector> JnVector;

    //! df/dp
    Teuchos::RCP<const NOX::Abstract::Vector> dfdp;

    //! d(Jn)/dp
    Teuchos::RCP<const NOX::Abstract::Vector> dJndp;

      }; // Class SalingerBordering

    } // Namespace MooreSpence

  } // Namespace Pitchfork

} // Namespace LOCA

#endif // LOCA_PITCHFORK_MOORESPENCE_SALINGERBORDERING_H
