// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_HOPF_MOORESPENCE_SALINGERBORDERING_H
#define LOCA_HOPF_MOORESPENCE_SALINGERBORDERING_H

#include "LOCA_Hopf_MooreSpence_SolverStrategy.H"
#include "NOX_Abstract_MultiVector.H"

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace Hopf {

    namespace MooreSpence {

      /*!
       * \brief Moore-Spence Hopf solver strategy based on
       * "Salinger" bordering.  This is the classic 5-solve Hopf bordering
       * method.
       */
      /*!
       * This class solves the Moore-Spence Hopf Newton equations:
       * \f[
       *     \begin{bmatrix}
       *          J       &    0   &    0   &   0  &    f_p     \\
       *       (Jy-wBz)_x &    J   &  -wB   & -Bz  & (Jy-wBz)_p \\
       *       (Jz+wBy)_x &   wB   &    J   &  By  & (Jz+wBy)_p \\
       *          0       & \phi^T &    0   &   0  &     0      \\
       *          0       &    0   & \phi^T &   0  &     0
       *     \end{bmatrix}
       *     \begin{bmatrix}
       *        X \\
       *        Y \\
       *        Z \\
       *        \omega \\
       *        \lambda
       *     \end{bmatrix} =
       *     \begin{bmatrix}
       *        F \\
       *        G \\
       *        H \\
       *        u \\
       *        v
       *     \end{bmatrix}.
       *  \f]
       *  via the following block elimination scheme:
       *  \f[
       *    \begin{split}
       *       J [A \; b] &= [F \; f_p] \\
       *       \begin{bmatrix}
       *          J & -wB \\
       *         wB &  J
       *       \end{bmatrix}
       *       \begin{bmatrix}
       *          C & e & g \\
       *          D & f & h
       *       \end{bmatrix} &=
       *       \begin{bmatrix}
       *          G - (Jy-wBz)_x A & (Jy-wBz)_p - (Jy-wBz)_x b & -Bz \\
       *          H - (Jz+wBy)_x A & (Jz+wBy)_p - (Jy-wBz)_x b &  By
       *       \end{bmatrix} \\
       *       \lambda &= \frac{(\phi^T h)(\phi^T C-u)-(\phi^T g)(\phi^T D-v)}
       *                       {(\phi^T h)(\phi^T e)-(\phi^T g)(\phi^T f)} \\
       *       \omega &= \frac{\phi^T D - v - (\phi^T f)\lambda}{\phi^T h} \\
       *       X &= A - b \lambda \\
       *       \begin{bmatrix}
       *          Y \\
       *          Z
       *       \end{bmatrix} &=
       *       \begin{bmatrix}
       *          C \\
       *          D
       *       \end{bmatrix} -
       *       \begin{bmatrix}
       *          e \\
       *          f
       *       \end{bmatrix}\lambda -
       *       \begin{bmatrix}
       *          g \\
       *          h
       *       \end{bmatrix}\omega
       *    \end{split}
       *   \f]
       */
      class SalingerBordering :
    public LOCA::Hopf::MooreSpence::SolverStrategy {

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
     * \param hopfGroup [in] Hopf group representing the Hopf equations.
     * \param yVector [in] Vector representing y
     * \param zVector [in] Vector representing z
     * \param CeRealVector [in] Vector representing Jy-wBz
     * \param CeImagVector [in] Vector representing Jz+wBy
     * \param dfdp [in] Vector representing df/dp
     * \param dCedpReal [in] Vector representing d(Jy-wBz)/dp
     * \param dCedpImag [in] Vector representing d(Jz+wBy)/dp
     * \param ByVector [in] Vector representing By
     * \param mBzVector [in] Vector representing -Bz
     * \param w [in] Bifurcation frequency w
     */
    virtual void setBlocks(
     const Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup>& group,
     const Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedGroup>& hopfGroup,
     const Teuchos::RCP<const NOX::Abstract::Vector>& yVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& zVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& CeRealVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& CeImagVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& dfdp,
     const Teuchos::RCP<const NOX::Abstract::Vector>& dCedpReal,
     const Teuchos::RCP<const NOX::Abstract::Vector>& dCedpImag,
     const Teuchos::RCP<const NOX::Abstract::Vector>& ByVector,
     const Teuchos::RCP<const NOX::Abstract::Vector>& mBzVector,
     double w);

    /*!
     * \brief Solves the extended system as defined above
     */
    /*!
     * The \em params argument is the linear solver parameters.
     */
    virtual NOX::Abstract::Group::ReturnType
    solve(
       Teuchos::ParameterList& params,
       const LOCA::Hopf::MooreSpence::ExtendedMultiVector& input,
       LOCA::Hopf::MooreSpence::ExtendedMultiVector& result) const;

      protected:

    //! Solves equations with contiguous arguments
    NOX::Abstract::Group::ReturnType
    solveContiguous(
                  Teuchos::ParameterList& params,
              const NOX::Abstract::MultiVector& input_x,
              const NOX::Abstract::MultiVector& input_y,
              const NOX::Abstract::MultiVector& input_z,
              const NOX::Abstract::MultiVector::DenseMatrix& input_w,
              const NOX::Abstract::MultiVector::DenseMatrix& input_p,
              NOX::Abstract::MultiVector& result_x,
              NOX::Abstract::MultiVector& result_y,
              NOX::Abstract::MultiVector& result_z,
              NOX::Abstract::MultiVector::DenseMatrix& result_w,
                  NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

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
    Teuchos::RCP<LOCA::Hopf::MooreSpence::AbstractGroup> group;

    //! Hopf group
    Teuchos::RCP<LOCA::Hopf::MooreSpence::ExtendedGroup> hopfGroup;

    //! Real component of eigenvector
    Teuchos::RCP<const NOX::Abstract::Vector> yVector;

    //! Imaginary component of eigenvector
    Teuchos::RCP<const NOX::Abstract::Vector> zVector;

    //! Jy-wBz vector
    Teuchos::RCP<const NOX::Abstract::Vector> CeRealVector;

    //! Jz+wBy vector
    Teuchos::RCP<const NOX::Abstract::Vector> CeImagVector;

    //! df/dp vector
    Teuchos::RCP<const NOX::Abstract::Vector> dfdp;

    //! d(Jy-wBz)/dp vector
    Teuchos::RCP<const NOX::Abstract::Vector> dCedpReal;

    //! d(Jz+wBy)/dp vector
    Teuchos::RCP<const NOX::Abstract::Vector> dCedpImag;

    //! By vector
    Teuchos::RCP<const NOX::Abstract::Vector> ByVector;

    //! -Bz vector
    Teuchos::RCP<const NOX::Abstract::Vector> minusBzVector;

    //! Frequency w
    double w;

      }; // Class SalingerBordering

    } // Namespace MooreSpence

  } // Namespace Hopf

} // Namespace LOCA

#endif // LOCA_HOPF_MOORESPENCE_SALINGERBORDERING_H
