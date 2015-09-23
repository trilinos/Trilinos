/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

#ifndef __Teko_ALOperator_hpp__
#define __Teko_ALOperator_hpp__

#include "Teko_BlockedEpetraOperator.hpp"
#include "Teko_Utilities.hpp"

namespace Teko
{

namespace NS
{

/** \brief Sparse matrix vector multiplication
 * for augmented Lagrangian-based preconditioners.
 *
 * This class implements sparse matrix vector multiplication
 * for augmented Lagrangian-based preconditioners.
 * Details can be found in the following papers:
 *
 * [1] M. Benzi and M. A. Olshanskii,
 * An Augmented Lagrangian-Based Approach to the Oseen Problem,
 * SIAM J. Scientific Computing, 28 (2006), pp. 2095-2113.
 *
 * [2] Benzi, M. A. Olshanskii and Z. Wang,
 * Modified Augmented Lagrangian Preconditioners for the Incompressible Navier-Stokes Equations,
 * International Journal for Numerical Methods in Fluids, 66 (2011), pp. 486-508.
 *
 * Suppose we are solving the following linear system:
 *
 * \f$
 * \left[ \begin{array}{cc}
 * A & B^T \\
 * B & -C
 * \end{array} \right]
 * \left[ \begin{array}{c}
 * u \\
  * p
 * \end{array} \right]
 * =
 * \left[ \begin{array}{c}
 * f \\
 * g
 * \end{array} \right].
 * \f$
 *
 * The equivalent augmented Lagrangian formulation is:
 *
 * \f$
 * \left[ \begin{array}{cc}
 * A + \gamma B^T W^{-1} B & B^T - \gamma B^T W^{-1} C \\
 * B & -C
 * \end{array} \right]
 * \left[ \begin{array}{c}
 * u \\
 * p
 * \end{array} \right]
 * =
 * \left[ \begin{array}{c}
 * f + \gamma B^T W^{-1} g \\
 * g
 * \end{array} \right]
 * \f$
 *
 * or
 *
 * \f$
 * \widehat{\mathcal{A}} x = \hat{b}.
 * \f$
 *
 * Here \f$ W \f$ can be take as the diagonal of the pressure
 * mass matrix and \f$ \gamma \f$ is a positive number.
 *
 * This class implements the matrix vector product with
 * \f$ \widehat{\mathcal{A}} \f$.
 */

class ALOperator : public Teko::Epetra::BlockedEpetraOperator
{
public:

   /** Build an augmented Lagrangian operator based on a vector of vector
    * of global IDs.
    *
    * \param[in] vars
    *            Vector of vectors of global ids specifying
    *            how the operator is to be blocked.
    * \param[in] content
    *            Operator to be blocked
    * \param[in] pressureMassMatrix
    *            Pressure mass matrix
    * \param[in] gamma
    *            Augmentation parameter
    * \param[in] label
    *            Label for name the operator
    */
   ALOperator(const std::vector<std::vector<int> > & vars,
         const Teuchos::RCP<Epetra_Operator> & content,
         LinearOp pressureMassMatrix,
         double gamma = 0.05, const std::string & label = "<ANYM>");

   /** Build a modified augmented Lagrangian operator based on a vector of vector
    * of global IDs.
    *
    * \param[in] vars
    *            Vector of vectors of global ids specifying
    *            how the operator is to be blocked.
    * \param[in] content
    *            Operator to be blocked
    * \param[in] gamma
    *            Augmentation parameter
    * \param[in] label
    *            Name of the operator
    */
   ALOperator(const std::vector<std::vector<int> > & vars,
         const Teuchos::RCP<Epetra_Operator> & content,
         double gamma = 0.05, const std::string & label = "<ANYM>");

   // Destructor
   virtual
   ~ALOperator()
   {
   }

   /** Set the pressure mass matrix.
    *
    * \param[in] pressureMassMatrix
    *            Pressure mass matrix.
    *
    */
   void
   setPressureMassMatrix(LinearOp pressureMassMatrix);

   /**
    * \returns Pressure mass matrix that can be used to construct preconditioner.
    */
   const LinearOp &
   getPressureMassMatrix() const
   {
      return pressureMassMatrix_;
   }

   /** Set gamma.
    *
    * \param[in] gamma
    *            Augmentation parameter.
    */
   void
   setGamma(double gamma);

   /**
    * \returns Gamma
    *          Augmentation parameter.
    */
   const double &
   getGamma() const
   {
      return gamma_;
   }

   /**
    * \param[in] b
    *            Right-hand side.
    * \param[out] bAugmented
    *             Augmented right-hand side.
    */
   void
   augmentRHS(const Epetra_MultiVector & b, Epetra_MultiVector & bAugmented);

   /**
    * \returns Number of rows.
    */
   int
   getNumberOfBlockRows() const
   {
      return numBlockRows_;
   }

   /**
    * Force a rebuild of the blocked operator from the stored
    * content operator.
    */
   virtual void
   RebuildOps()
   {
      BuildALOperator();
   }

   /** Get the (i,j) block of the original (non-augmented) operator.
    *
    * \param[in] i
    *            Row index.
    * \param[in] j
    *            Column index.
    */
   const Teuchos::RCP<const Epetra_Operator>
   GetBlock(int i, int j) const;

protected:

   /**
    * AL operator.
    */
   Teuchos::RCP<Thyra::LinearOpBase<double> > alOperator_;

   /**
    * Operator for augmenting the right-hand side.
    */
   Teuchos::RCP<Thyra::LinearOpBase<double> > alOperatorRhs_;

   /**
    * Pressure mass matrix and its inverse.
    */
   LinearOp pressureMassMatrix_;

   /**
    * Inverse of the pressure mass matrix.
    */
   LinearOp invPressureMassMatrix_;

   /**
    * Augmentation parameter.
    */
   double gamma_;

   /**
    * Dimension of the problem.
    */
   int dim_;

   /**
    * Number of block rows.
    */
   int numBlockRows_;

   /**
    * Check dimension. Only implemente for 2D and 3D problems.
    */
   void
   checkDim(const std::vector<std::vector<int> > & vars);

   /**
    * Build AL operator.
    */
   void
   BuildALOperator();
};

} // end namespace NS

} // end namespace Teko

#endif /* __Teko_ALOperator_hpp__ */
