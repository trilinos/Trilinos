#ifndef __PB_Helpers_hpp__
#define __PB_Helpers_hpp__

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace PB {

Teuchos::RCP<Thyra::LinearOpBase<double> > inverse(const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > & lowsFact,
                                                   const Teuchos::RCP<const Thyra::PreconditionerFactoryBase<double> > & precFact,
                                                   const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A);


/** \brief Compute the spectral radius of a matrix
  *
  * Compute the spectral radius of matrix A.  This utilizes the 
  * Trilinos-Anasazi BlockKrylovShcur method for computing eigenvalues.
  * It attempts to compute the largest (in magnitude) eigenvalue to a given
  * level of tolerance.
  *
  * \param[in] A   matrix whose spectral radius is needed
  * \param[in] tol The <em>most</em> accuracy needed (the algorithm will run until
  *            it reaches this level of accuracy and then it will quit).
  *            If this level is not reached it will return something to indicate
  *            it has not converged.
  * \param[in] isHermitian Is the matrix Hermitian
  * \param[in] numBlocks The size of the orthogonal basis built (like in GMRES) before
  *                  restarting.  Increase the memory usage by O(restart*n). At least
  *                  restart=3 is required.
  * \param[in] restart How many restarts are permitted
  * \param[in] verbosity See the Anasazi documentation
  *
  * \return The spectral radius of the matrix.  If the algorithm didn't converge the
  *         number is the negative of the ritz-values. If a <code>NaN</code> is returned
  *         there was a problem constructing the Anasazi problem
  */
double computeSpectralRad(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A,double tol,
                          bool isHermitian=false,int numBlocks=5,int restart=0,int verbosity=0);

} // end namespace PB

#endif
