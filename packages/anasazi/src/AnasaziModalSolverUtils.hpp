// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_MODAL_SOLVER_UTILS_HPP
#define ANASAZI_MODAL_SOLVER_UTILS_HPP

/*!     \file AnasaziModalSolverUtils.hpp
        \brief Class which provides internal utilities for the Anasazi modal solvers.
*/

/*!    \class Anasazi::ModalSolverUtils
       \brief Anasazi's templated, static class providing utilities for
       the modal solvers.

       This class provides concrete, templated implementations of utilities necessary
       for the modal solvers (Davidson, LOBPCG, ...).  These utilities include
       sorting, orthogonalization, projecting/solving local eigensystems, and sanity
       checking.  These are internal utilties, so the user should not alter this class.

       \author Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"
#include "AnasaziOutputManager.hpp"
#include "Teuchos_BLAS.hpp"

namespace Anasazi {

  template<class ScalarType, class MV, class OP>
  class ModalSolverUtils 
  {  
  public:
    
    //@{ \name Constructor/Destructor

    ModalSolverUtils( const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om ); 
    
    virtual ~ModalSolverUtils() {};
    
    //@}
    
    //@{ \name Sorting Methods
    
    int sortScalars(int n, ScalarType *y, int *perm = 0) const;
    
    int sortScalars_Vectors(int n, ScalarType* lambda, MV* Q, std::vector<ScalarType>* resids = 0) const;

    //@} 

    //@{ \name Eigensolver Projection Methods

    int massOrthonormalize(MV &X, MV &MX, const OP *M, const MV &Q, int howMany,
			   int orthoType = 0, ScalarType kappa = 1.5625) const;
    
    int directSolver(int size, const Teuchos::SerialDenseMatrix<int,ScalarType> &KK, 
		     const Teuchos::SerialDenseMatrix<int,ScalarType> *MM,
		     Teuchos::SerialDenseMatrix<int,ScalarType> *EV,
		     std::vector<ScalarType>* theta,
		     int nev, int esType = 0) const;

    //@}

    //@{ \name Sanity Checking Methods

    ScalarType errorOrthogonality(const MV *X, const MV *R, const OP *M = 0) const;
    
    ScalarType errorOrthonormality(const MV *X, const OP *M = 0) const;
    
    ScalarType errorEquality(const MV *X, const MV *MX, const OP *M = 0) const;
    
    //@}
    
  private:

    // Reference counted pointer to output manager used by eigensolver.
    Teuchos::RefCountPtr<OutputManager<ScalarType> > _om;

    //@{ \name Internal Typedefs

    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef OperatorTraits<ScalarType,MV,OP> OPT;

    //@}
  };

  //-----------------------------------------------------------------------------
  // 
  //  CONSTRUCTOR
  //
  //-----------------------------------------------------------------------------  

  template<class ScalarType, class MV, class OP>
  ModalSolverUtils<ScalarType, MV, OP>::ModalSolverUtils( const Teuchos::RefCountPtr<OutputManager<ScalarType> > &om ) 
    : _om(om)
  {}

  //-----------------------------------------------------------------------------
  // 
  //  SORTING METHODS
  //
  //-----------------------------------------------------------------------------
  
  template<class ScalarType, class MV, class OP>
  int ModalSolverUtils<ScalarType, MV, OP>::sortScalars( int n, ScalarType *y, int *perm) const 
  {
    // Sort a vector into increasing order of algebraic values
    //
    // Input:
    //
    // n    (integer ) = Size of the array (input)
    // y    (double* ) = Array of length n to be sorted (input/output)
    // perm (integer*) = Array of length n with the permutation (input/output)
    //                   Optional argument
    
    int i, j;
    int igap = n / 2;
    
    if (igap == 0) {
      if ((n > 0) && (perm != 0)) {
	perm[0] = 0;
      }
      return 0;
    }
    
    if (perm) {
      for (i = 0; i < n; ++i)
	perm[i] = i;
    }
    
    while (igap > 0) {
      for (i=igap; i<n; ++i) {
	for (j=i-igap; j>=0; j-=igap) {
	  if (y[j] > y[j+igap]) {
	    double tmpD = y[j];
	    y[j] = y[j+igap];
	    y[j+igap] = tmpD;
	    if (perm) {
	      int tmpI = perm[j];
	      perm[j] = perm[j+igap];
	      perm[j+igap] = tmpI;
	    }
	  }
	  else {
	    break;
	  }
	}
      }
      igap = igap / 2;
    }
    
    return 0;
  }
  
  template<class ScalarType, class MV, class OP>
  int ModalSolverUtils<ScalarType, MV, OP>::sortScalars_Vectors( int n, ScalarType *lambda, MV *Q, 
								 std::vector<ScalarType> *resids ) const
  {
    // This routines sorts the scalars (stored in lambda) in ascending order.
    // The associated vectors (stored in Q) are accordingly ordered.
    
    int info = 0;
    int i, j;
    std::vector<int> index(1);
    ScalarType tmp;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    
    if ( n > MVT::GetNumberVecs( *Q ) ) { return -1; }

    int igap = n / 2;
    
    if ( Q ) {
      while (igap > 0) {
	for (i=igap; i < n; ++i) {
	  for (j=i-igap; j>=0; j-=igap) {
	    if (lambda[j] > lambda[j+igap]) {

	      // Swap two scalars
	      tmp = lambda[j];
	      lambda[j] = lambda[j+igap];
	      lambda[j+igap] = tmp;

	      // Swap residuals (if they exist)
	      if (resids) {
		tmp = (*resids)[j];
		(*resids)[j] = (*resids)[j+igap];
		(*resids)[j+igap] = tmp;
	      }

	      // Swap corresponding vectors
	      index[0] = j;
	      Teuchos::RefCountPtr<MV> tmpQ = MVT::CloneCopy( *Q, index );
	      Teuchos::RefCountPtr<MV> tmpQj = MVT::CloneView( *Q, index );
	      index[0] = j + igap;
	      Teuchos::RefCountPtr<MV> tmpQgap = MVT::CloneView( *Q, index );
	      MVT::MvAddMv( one, *tmpQgap, zero, *tmpQgap, *tmpQj );
	      MVT::MvAddMv( one, *tmpQ, zero, *tmpQ, *tmpQgap );
	    } 
	    else {
	      break;
	    }
	  }
	}
	igap = igap / 2;
      } // while (igap > 0)
    } // if ( Q )
    else {
      while (igap > 0) {
	for (i=igap; i < n; ++i) {
	  for (j=i-igap; j>=0; j-=igap) {
	    if (lambda[j] > lambda[j+igap]) {
	      // Swap two scalars
	      tmp = lambda[j];
	      lambda[j] = lambda[j+igap];
	      lambda[j+igap] = tmp;
	      
	      // Swap residuals (if they exist)
	      if (resids) {
		tmp = (*resids)[j];
		(*resids)[j] = (*resids)[j+igap];
		(*resids)[j+igap] = tmp;
	      }	      
	    } 
	    else {
	      break;
	    }
	  }
	}
	igap = igap / 2;
      } // while (igap > 0)
    } // if ( Q )
    
    return info;  
  }
  
  //-----------------------------------------------------------------------------
  // 
  //  EIGENSOLVER PROJECTION METHODS
  //
  //-----------------------------------------------------------------------------
  
  template<class ScalarType, class MV, class OP>
  int ModalSolverUtils<ScalarType, MV, OP>::massOrthonormalize(MV &X, MV &MX, const OP *M, const MV &Q, 
							       int howMany, int orthoType, ScalarType kappa) const
  {
    // For the inner product defined by the operator M or the identity (M = 0)
    //   -> Orthogonalize X against Q
    //   -> Orthonormalize X 
    // Modify MX accordingly
    //
    // Note that when M is 0, MX is not referenced
    //
    // Parameter variables
    //
    // X  : Vectors to be transformed
    //
    // MX : Image of the block vector X by the mass matrix
    //
    // Q  : Vectors to orthogonalize against
    //
    // howMany : Number of vectors of X to orthogonalized
    //           If this number is smaller than the total number of vectors in X,
    //           then it is assumed that the last "howMany" vectors are not orthogonal
    //           while the other vectors in X are othogonal to Q and orthonormal.
    //
    // orthoType = 0 (default) > Performs both operations
    // orthoType = 1           > Performs Q^T M X = 0
    // orthoType = 2           > Performs X^T M X = I
    //
    // kappa= Coefficient determining when to perform a second Gram-Schmidt step
    //        Default value = 1.5625 = (1.25)^2 (as suggested in Parlett's book)
    //
    // Output the status of the computation
    //
    // info =   0 >> Success.
    //
    // info >   0 >> Indicate how many vectors have been tried to avoid rank deficiency for X 
    //
    // info =  -1 >> Failure >> X has zero columns
    //                       >> It happens when # col of X     > # rows of X 
    //                       >> It happens when # col of [Q X] > # rows of X 
    //                       >> It happens when no good random vectors could be found
    
    int i;
    int info = 0;
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    ScalarType eps = Teuchos::ScalarTraits<ScalarType>::eps();
    
    // Orthogonalize X against Q
    //timeProj -= MyWatch.WallTime();
    if (orthoType != 2) {
      
      int xc = howMany;
      int xr = MVT::GetNumberVecs( X );
      int qc = MVT::GetNumberVecs( Q );
      
      std::vector<int> index(howMany);
      for (i=0; i<howMany; i++)
	index[i] = xr - howMany + i;

      Teuchos::RefCountPtr<MV> XX = MVT::CloneView( X, index );
      
      Teuchos::RefCountPtr<MV> MXX;

      if (M) {
	int mxr = MVT::GetNumberVecs( MX );
	for (i=0; i<howMany; i++)
	  index[i] = mxr - howMany + i;
	MXX = MVT::CloneView( MX, index );
      } 
      else {
	MXX = MVT::CloneView( X, index );
      }

      // Perform the Gram-Schmidt transformation for a block of vectors
      
      // Compute the initial M-norms
      std::vector<ScalarType> oldDot( xc );
      MVT::MvDot( *XX, *MXX, &oldDot );
      
      // Define the product Q^T * (M*X)
      // Multiply Q' with MX
      Teuchos::SerialDenseMatrix<int,ScalarType> qTmx( qc, xc );
      MVT::MvTransMv( one, Q, MX, qTmx );
      
      // Multiply by Q and substract the result in X
      MVT::MvTimesMatAddMv( -one, Q, qTmx, one, *XX );
      
      // Update MX
      if (M) {
	if (qc >= xc) {
	  //timeProj_MassMult -= MyWatch.WallTime();
	  OPT::Apply( *M, *XX, *MXX);
	  //timeProj_MassMult += MyWatch.WallTime();
	  //numProj_MassMult += xc;
	}
	else {
	  Teuchos::RefCountPtr<MV> MQ = MVT::Clone( Q, qc );
	  //timeProj_MassMult -= MyWatch.WallTime();
	  OPT::Apply( *M, Q, *MQ );
	  //timeProj_MassMult += MyWatch.WallTime();
	  //numProj_MassMult += qc;
	  MVT::MvTimesMatAddMv( -one, *MQ, qTmx, one, *MXX );
	}  // if (qc >= xc)
      } // if (M)
      
      // Compute new M-norms
      std::vector<ScalarType> newDot(xc);
      MVT::MvDot( *XX, *MXX, &newDot );
      
      int j;
      for (j = 0; j < xc; ++j) {
	
	if (kappa*newDot[j] < oldDot[j]) {
	  
	  // Apply another step of classical Gram-Schmidt
	  //timeQtMult -= MyWatch.WallTime();
	  MVT::MvTransMv( one, Q, *MXX, qTmx );
	  //timeQtMult += MyWatch.WallTime();
	  
	  //timeQMult -= MyWatch.WallTime();
	  MVT::MvTimesMatAddMv( -one, Q, qTmx, one, *XX );
	  //timeQMult += MyWatch.WallTime();
	  
	  // Update MX
	  if (M) {
	    if (qc >= xc) {
	      //timeProj_MassMult -= MyWatch.WallTime();
	      OPT::Apply( *M, *XX, *MXX);
	      //timeProj_MassMult += MyWatch.WallTime();
	      //numProj_MassMult += xc;
	    }
	    else {
	      Teuchos::RefCountPtr<MV> MQ = MVT::Clone( Q, qc );
	      //timeProj_MassMult -= MyWatch.WallTime();
	      OPT::Apply( *M, Q, *MQ);
	      //timeProj_MassMult += MyWatch.WallTime();
	      //numProj_MassMult += qc;
	      MVT::MvTimesMatAddMv( -one, *MQ, qTmx, one, *MXX );
	    } // if (qc >= xc)
	  } // if (M)
	  
	  break;
	} // if (kappa*newDot[j] < oldDot[j])
      } // for (j = 0; j < xc; ++j)
      
    } // if (orthoType != 2)

    //timeProj += MyWatch.WallTime();
    
    // Orthonormalize X 
    //  timeNorm -= MyWatch.WallTime();
    if (orthoType != 1) {
      
      int j;
      int xc = MVT::GetNumberVecs( X );
      int xr = MVT::GetVecLength( X );
      int shift = (orthoType == 2) ? 0 : MVT::GetNumberVecs( Q );
      int mxc = (M) ? MVT::GetNumberVecs( MX ) : xc;

      std::vector<int> index( 1 );      
      std::vector<ScalarType> oldDot( 1 ), newDot( 1 );

      for (j = 0; j < howMany; ++j) {
	
	int numX = xc - howMany + j;
	int numMX = mxc - howMany + j;
	
	// Put zero vectors in X when we are exceeding the space dimension
	if (numX + shift >= xr) {
	  index[0] = numX;
	  Teuchos::RefCountPtr<MV> XXj = MVT::CloneView( X, index );
	  MVT::MvInit( *XXj, zero );
	  if (M) {
	    index[0] = numMX;
	    Teuchos::RefCountPtr<MV> MXXj = MVT::CloneView( MX, index );
	    MVT::MvInit( *MXXj, zero );
	  }
	  info = -1;
	}
	
	// Get a view of the vectors currently being worked on.
	index[0] = numX;
	Teuchos::RefCountPtr<MV> Xj = MVT::CloneView( X, index );
	index[0] = numMX;
	Teuchos::RefCountPtr<MV> MXj;
	if (M)
	  MXj = MVT::CloneView( MX, index );
	else
	  MXj = MVT::CloneView( X, index );

	// Get a view of the previous vectors.
	std::vector<int> prev_idx( numX );
	Teuchos::RefCountPtr<MV> prevXj;

	if (numX > 0) {
	  for (i=0; i<numX; i++)
	    prev_idx[i] = i;
	  prevXj = MVT::CloneView( X, prev_idx );
	} 

	// Make storage for these Gram-Schmidt iterations.
	Teuchos::SerialDenseMatrix<int,ScalarType> product( numX, 1 );

	int numTrials;
	bool rankDef = true;
	for (numTrials = 0; numTrials < 10; ++numTrials) {

	  //
	  // Compute M-norm
	  //
	  MVT::MvDot( *Xj, *MXj, &oldDot );
	  //
	  // Save old MXj vector.
	  //
	  Teuchos::RefCountPtr<MV> oldMXj = MVT::CloneCopy( *MXj );

	  if (numX > 0) {
	    //
	    // Apply the first step of Gram-Schmidt
	    //
	    MVT::MvTransMv( one, *prevXj, *MXj, product );
	    //
	    MVT::MvTimesMatAddMv( -one, *prevXj, product, one, *Xj );
	    //
	    if (M) {
	      if (xc == mxc) {
		Teuchos::RefCountPtr<MV> prevMXj = MVT::CloneView( MX, prev_idx );
		MVT::MvTimesMatAddMv( -one, *prevMXj, product, one, *MXj );
	      }
	      else {
		//timeNorm_MassMult -= MyWatch.WallTime();
		OPT::Apply( *M, *Xj, *MXj );
		//timeNorm_MassMult += MyWatch.WallTime();
		//numNorm_MassMult += 1;
	      }
	    }
	    //
	    // Compute new M-norm
	    //
	    MVT::MvDot( *Xj, *MXj, &newDot );
	    //
	    // Check if a correction is needed.
	    //
	    if (kappa*newDot[0] < oldDot[0]) {
	      //
	      // Apply the second step of Gram-Schmidt
	      //
	      MVT::MvTransMv( one, *prevXj, *MXj, product );
	      //
	      MVT::MvTimesMatAddMv( -one, *prevXj, product, one, *Xj );
	      //
	      if (M) {
		if (xc == mxc) {
		  Teuchos::RefCountPtr<MV> prevMXj = MVT::CloneView( MX, prev_idx );
		  MVT::MvTimesMatAddMv( -one, *prevMXj, product, one, *MXj );
		}
		else {
		  //timeNorm_MassMult -= MyWatch.WallTime();
		  OPT::Apply( *M, *Xj, *MXj );
		  //timeNorm_MassMult += MyWatch.WallTime();
		  //numNorm_MassMult += 1;
		}
	      }
	    } // if (kappa*newDot[0] < oldDot[0])
	    
	  } // if (numX > 0)
	    
	  //
	  // Compute M-norm with old MXj
	  // NOTE:  Re-using newDot vector
	  //
	  MVT::MvDot( *Xj, *oldMXj, &newDot );
	  
	  if (newDot[0] > oldDot[0]*eps*eps) {
	    MVT::MvAddMv( one/Teuchos::ScalarTraits<ScalarType>::squareroot(newDot[0]), *Xj, zero, *Xj, *Xj );
	    if (M)
	      MVT::MvAddMv( one/Teuchos::ScalarTraits<ScalarType>::squareroot(newDot[0]), *MXj, zero, *MXj, *MXj );
	    rankDef = false;
	    break;
	  }
	  else {
	    info += 1;
	    MVT::MvRandom( *Xj );
	    if (M) {
	      //timeNorm_MassMult -= MyWatch.WallTime();
	      OPT::Apply( *M, *Xj, *MXj );
	      //timeNorm_MassMult += MyWatch.WallTime();
	      //numNorm_MassMult += 1;
	    }
	    if (orthoType == 0)
	      massOrthonormalize(*Xj, *MXj, M, Q, 1, 1, kappa);
	  } // if (norm > oldDot*eps*eps)
	  
	}  // for (numTrials = 0; numTrials < 10; ++numTrials)
	
	if (rankDef == true) {
	  MVT::MvInit( *Xj, zero );
	  if (M)
	    MVT::MvInit( *MXj, zero );
	  info = -1;
	  break;
	}
	
      } // for (j = 0; j < howMany; ++j)
	
    } // if (orthoType != 1)
    //timeNorm += MyWatch.WallTime();

    return info;
    
  }
  
  template<class ScalarType, class MV, class OP>
  int ModalSolverUtils<ScalarType, MV, OP>::directSolver(int size, const Teuchos::SerialDenseMatrix<int,ScalarType> &KK, 
							 const Teuchos::SerialDenseMatrix<int,ScalarType> *MM,
							 Teuchos::SerialDenseMatrix<int,ScalarType>* EV,
							 std::vector<ScalarType>* theta,
							 int nev, int esType) const
  {
    // Routine for computing the first NEV generalized eigenpairs of the symmetric pencil (KK, MM)
    //
    // Parameter variables:
    //
    // size : Dimension of the eigenproblem (KK, MM)
    //
    // KK : Symmetric "stiffness" matrix 
    //
    // MM : Symmetric Positive "mass" matrix
    //
    // EV : Array to store the nev eigenvectors 
    //
    // theta : Array to store the eigenvalues (Size = nev )
    //
    // nev : Number of the smallest eigenvalues requested (input)
    //       Number of the smallest computed eigenvalues (output)
    //
    // esType : Flag to select the algorithm
    //
    // esType =  0 (default) Uses LAPACK routine (Cholesky factorization of MM)
    //                       with deflation of MM to get orthonormality of 
    //                       eigenvectors (S^T MM S = I)
    //
    // esType =  1           Uses LAPACK routine (Cholesky factorization of MM)
    //                       (no check of orthonormality)
    //
    // esType = 10           Uses LAPACK routine for simple eigenproblem on KK
    //                       (MM is not referenced in this case)
    //
    // Note: The code accesses only the upper triangular part of KK and MM.
    //
    // Return the integer info on the status of the computation
    //
    // info = 0 >> Success
    //
    // info = - 20 >> Failure in LAPACK routine
    
    // Define local arrays

    // Create blas/lapack objects.
    Teuchos::LAPACK<int,ScalarType> lapack;
    Teuchos::BLAS<int,ScalarType> blas;
    
    int i, j;
    int rank = 0;
    int info = 0;
   
    std::string lapack_name = "sytrd";
    std::string lapack_opts = "u";
    int NB = 5 + lapack.ILAENV(1, lapack_name, lapack_opts, size, -1, -1, -1);
    int lwork = size*NB;
    std::vector<ScalarType> work(lwork);
    std::vector<ScalarType> tt( size );
    
    //  ScalarType tol = sqrt(eps);
    ScalarType tol = 1e-12;
    ScalarType zero = Teuchos::ScalarTraits<ScalarType>::zero();
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();

    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > KKcopy, MMcopy;
    Teuchos::RefCountPtr<Teuchos::SerialDenseMatrix<int,ScalarType> > U;
 
    switch (esType) {
      
    default:
    case 0:
      //
      // Use the Cholesky factorization of MM to compute the generalized eigenvectors
      //
      for (rank = size; rank > 0; --rank) {

	U = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>(size,size) );      
	//
	// Copy KK & MM
	//
	KKcopy = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, KK, rank, rank ) );
	MMcopy = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, *MM, rank, rank ) );
	//
	// Solve the generalized eigenproblem with LAPACK
	//
	info = 0;
	lapack.SYGV(1, 'V', 'U', rank, KKcopy->values(), KKcopy->stride(), 
		    MMcopy->values(), MMcopy->stride(), &tt[0], &work[0], lwork, &info);
	//
	// Treat error messages
	//
	if (info < 0) {
	  //	    if (verbose > 0) {
	  cerr << endl;
	  cerr << " In DSYGV, argument " << -info << "has an illegal value.\n";
	  cerr << endl;
	  //}
	  return -20;
	}
	if (info > 0) {
	  if (info > rank)
	    rank = info - rank;
	  continue;
	}
	//
	// Check the quality of eigenvectors
	// ( using mass-orthonormality )
	//
	Teuchos::SerialDenseMatrix<int,ScalarType> MMcopy2( Teuchos::Copy, *MM, size, size );	  
	for (i = 0; i < size; ++i) {
	  for (j = 0; j < i; ++j)
	    MMcopy2(i,j) = (*MM)(j,i);
	}
	blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, size, rank, size, one, MMcopy2.values(), MMcopy2.stride(), 
		  KKcopy->values(), KKcopy->stride(), zero, U->values(), U->stride());
	blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, rank, rank, size, one, KKcopy->values(), KKcopy->stride(), 
		  U->values(), U->stride(), zero, MMcopy2.values(), MMcopy2.stride());
	ScalarType maxNorm = zero;
	ScalarType maxOrth = zero;
	for (i = 0; i < rank; ++i) {
	  for (j = i; j < rank; ++j) {
	    if (j == i)
	      maxNorm = (Teuchos::ScalarTraits<ScalarType>::magnitude(MMcopy2(i,j)-one) > maxNorm) 
		? Teuchos::ScalarTraits<ScalarType>::magnitude(MMcopy2(i,j)-one) : maxNorm;	    
	    else 
	      maxOrth = (Teuchos::ScalarTraits<ScalarType>::magnitude(MMcopy2(i,j)) > maxOrth)
		? Teuchos::ScalarTraits<ScalarType>::magnitude(MMcopy2(i,j)) : maxOrth;
	  }
	}
	/*        if (verbose > 4) {
	cout << " >> Local eigensolve >> Size: " << rank;
	cout.precision(2);
	cout.setf(ios::scientific, ios::floatfield);
	cout << " Normalization error: " << maxNorm;
	cout << " Orthogonality error: " << maxOrth;
	cout << endl;
	}*/
	if ((maxNorm <= tol) && (maxOrth <= tol))
	  break;
      } // for (rank = size; rank > 0; --rank)
      //
      // Copy the computed eigenvectors and eigenvalues
      // ( they may be less than the number requested because of deflation )
      //
      nev = (rank < nev) ? rank : nev;
      EV->putScalar( zero );
      blas.COPY( nev, &tt[0], 1, &(*theta)[0], 1 );
      for (i = 0; i < nev; ++i) {
	blas.COPY( rank, (*KKcopy)[i], 1, (*EV)[i], 1 );
      }
      
      break;
      
    case 1:
      //
      // Use the Cholesky factorization of MM to compute the generalized eigenvectors
      //
      // Copy KK & MM
      //
      KKcopy = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, KK, size, size ) );
      MMcopy = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, *MM, size, size ) );
      //
      // Solve the generalized eigenproblem with LAPACK
      //
      info = 0;
      lapack.SYGV(1, 'V', 'U', size, KKcopy->values(), KKcopy->stride(), 
		  MMcopy->values(), MMcopy->stride(), &tt[0], &work[0], lwork, &info);
      //
      // Treat error messages
      //
      if (info < 0) {
	//if (verbose > 0) {
	cerr << endl;
	cerr << " In DSYGV, argument " << -info << "has an illegal value.\n";
	cerr << endl;
	//      }
	return -20;
      }
      if (info > 0) {
	if (info > size)
	  nev = 0;
	else {
	  //	if (verbose > 0) {
	  cerr << endl;
	  cerr << " In DSYGV, DPOTRF or DSYEV returned an error code (" << info << ").\n";
	  cerr << endl;
	  //          }
	  return -20; 
	}
      }
      //
      // Copy the eigenvectors and eigenvalues
      //
      blas.COPY( nev, &tt[0], 1, &(*theta)[0], 1 );
      for (i = 0; i < nev; ++i) {
	blas.COPY( size, (*KKcopy)[i], 1, (*EV)[i], 1 );
      }
      
      break;
      
    case 10:
      //
      // Simple eigenproblem
      //
      // Copy KK
      //
      KKcopy = Teuchos::rcp( new Teuchos::SerialDenseMatrix<int,ScalarType>( Teuchos::Copy, KK, size, size ) );
      //
      // Solve the generalized eigenproblem with LAPACK
      //
      lapack.SYEV('V', 'U', size, KKcopy->values(), KKcopy->stride(), &tt[0], &work[0], lwork, &info);
      //
      // Treat error messages
      if (info != 0) {
	//      if (verbose > 0) {
	cerr << endl;
	if (info < 0) 
	  cerr << " In DSYEV, argument " << -info << " has an illegal value\n";
	else
	  cerr << " In DSYEV, the algorithm failed to converge (" << info << ").\n";
	cerr << endl;
	//}
	info = -20;
	break;
      }
      //
      // Copy the eigenvectors
      //
      blas.COPY( nev, &tt[0], 1, &(*theta)[0], 1 );
      for (i = 0; i < nev; ++i) {
	blas.COPY( size, (*KKcopy)[i], 1, (*EV)[i], 1 );
      }
      
      break;
      
    }
    
    // Clear the memory
    
    return info;
  }
  
  //-----------------------------------------------------------------------------
  // 
  //  SANITY CHECKING METHODS
  //
  //-----------------------------------------------------------------------------

  template<class ScalarType, class MV, class OP>
  ScalarType ModalSolverUtils<ScalarType, MV, OP>::errorOrthogonality(const MV *X, const MV *R, 
								      const OP *M) const
  {
    // Return the maximum value of R_i^T * M * X_j / || MR_i || || X_j ||
    // When M is not specified, the identity is used.
    ScalarType maxDot = Teuchos::ScalarTraits<ScalarType>::zero();
    
    int xc = (X) ? MVT::GetNumberVecs( *X ) : 0;
    int rc = (R) ? MVT::GetNumberVecs( *R ) : 0;
    
    if (xc*rc == 0)
      return maxDot;
    
    int i, j;
    Teuchos::RefCountPtr<MV> MR;
    std::vector<ScalarType> normMR( rc );
    std::vector<ScalarType> normX( xc );
    if (M) {
      MR = MVT::Clone( *R, rc );
      OPT::Apply( *M, *R, *MR );
    }
    else {
      MR = MVT::CloneCopy( *R );
    }
    MVT::MvNorm( *MR, &normMR );
    MVT::MvNorm( *X, &normX );

    ScalarType dot = Teuchos::ScalarTraits<ScalarType>::zero();
    Teuchos::SerialDenseMatrix<int, ScalarType> xTMr( xc, rc );
    MVT::MvTransMv( 1.0, *X, *MR, xTMr );    
    for (i = 0; i < xc; ++i) {
      for (j = 0; j < rc; ++j) {
	dot = Teuchos::ScalarTraits<ScalarType>::magnitude(xTMr(i,j))/(normMR[j]*normX[i]);
	maxDot = (dot > maxDot) ? dot : maxDot;
      }
    }
    
    return maxDot;
    
  }
  
  template<class ScalarType, class MV, class OP>
  ScalarType ModalSolverUtils<ScalarType, MV, OP>::errorOrthonormality(const MV *X, const OP *M) const
  {
    // Return the maximum coefficient of the matrix X^T * M * X - I
    // When M is not specified, the identity is used.
    ScalarType maxDot = Teuchos::ScalarTraits<ScalarType>::zero();
    ScalarType one = Teuchos::ScalarTraits<ScalarType>::one();
    
    int xc = (X) ? MVT::GetNumberVecs( *X ) : 0;
    if (xc == 0)
      return maxDot;
    
    int i, j;
    std::vector<int> index( 1 );
    std::vector<ScalarType> dot( 1 );
    Teuchos::RefCountPtr<MV> MXi;
    Teuchos::RefCountPtr<const MV> Xi;

    // Create space if there is a M matrix specified.
    if (M)
      MXi = MVT::Clone( *X, 1 );

    for (i = 0; i < xc; ++i) {
      index[0] = i;
      if (M) {
	Xi = MVT::CloneView( *X, index );
	OPT::Apply( *M, *Xi, *MXi );
      }
      else {
	MXi = MVT::CloneView( *(const_cast<MV *>(X)), index );
      }
      for (j = 0; j < xc; ++j) {
	index[0] = j;
	Xi = MVT::CloneView( *X, index );
	MVT::MvDot( *Xi, *MXi, &dot );
	dot[0] = (i == j) ? fabs(dot[0] - one) : fabs(dot[0]);
	maxDot = (dot[0] > maxDot) ? dot[0] : maxDot;
      }
    }
    
    return maxDot;    
  }
  
  template<class ScalarType, class MV, class OP>
  ScalarType ModalSolverUtils<ScalarType, MV, OP>::errorEquality(const MV *X, const MV *MX, 
								 const OP *M) const
  {
    // Return the maximum coefficient of the matrix M * X - MX
    // scaled by the maximum coefficient of MX.
    // When M is not specified, the identity is used.
    
    ScalarType maxDiff = Teuchos::ScalarTraits<ScalarType>::zero();
    
    int xc = (X) ? MVT::GetNumberVecs( *X ) : 0;
    int mxc = (MX) ? MVT::GetNumberVecs( *MX ) : 0;
    
    if ((xc != mxc) || (xc*mxc == 0))
      return maxDiff;
    
    int i;
    ScalarType maxCoeffX = Teuchos::ScalarTraits<ScalarType>::zero();
    std::vector<ScalarType> tmp( xc );
    MVT::MvNorm( *MX, &tmp );

    for (i = 0; i < xc; ++i) {
      maxCoeffX = (tmp[i] > maxCoeffX) ? tmp[i] : maxCoeffX;
    }

    std::vector<int> index( 1 );
    Teuchos::RefCountPtr<MV> MtimesX; 
    if (M) {
      MtimesX = MVT::Clone( *X, xc );
      OPT::Apply( *M, *X, *MtimesX );
    }
    else {
      MtimesX = MVT::CloneCopy( *(const_cast<MV *>(X)) );
    }
    MVT::MvAddMv( -1.0, *MX, 1.0, *MtimesX, *MtimesX );
    MVT::MvNorm( *MtimesX, &tmp );
   
    for (i = 0; i < xc; ++i) {
      maxDiff = (tmp[i] > maxDiff) ? tmp[i] : maxDiff;
    }
    
    return (maxCoeffX == 0.0) ? maxDiff : maxDiff/maxCoeffX;
    
  }    
  
} // end namespace Anasazi

#endif // ANASAZI_MODAL_SOLVER_UTILS_HPP

