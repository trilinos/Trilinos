
#ifndef THYRA_CREATE_EXAMPLE_TRIDIAG_TPETRA_LINEAR_OP_HPP
#define THYRA_CREATE_EXAMPLE_TRIDIAG_TPETRA_LINEAR_OP_HPP

#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#ifdef TPETRA_MPI
#  include "Tpetra_MpiPlatform.hpp"
#else
#  include "Tpetra_SerialPlatform.hpp"
#endif

/** \brief \brief This function generates a tridiagonal linear operator using Tpetra.
 *
 * Specifically, this function returns a smart pointer to the matrix:
\f[

A=
\left[\begin{array}{rrrrrrrrrr}
2 a    & -1 \\
-1     &  2 a    & -1 \\
       & \ddots  & \ddots  & \ddots \\
       &         & -1      & 2 a     & -1 \\
       &         &         &  -1     & 2 a
\end{array}\right]
\f]
 *
 * where <tt>diagScale</tt> is \f$a\f$ and <tt>globalDim</tt> is the
 * glboal dimension of the matrix.
 */
template<class Ordinal, class Scalar>
Teuchos::RefCountPtr<Thyra::LinearOpBase<Scalar> >
createExampleTridiagTpetraLinearOp(
  const Ordinal      globalDim
#ifdef HAVE_MPI
  ,MPI_Comm          mpiComm
#endif
  ,const double      diagScale
  ,const bool        verbose
  ,std::ostream      &out
  )
{
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;

  //
  // (A) Create Tpetra::VectorSpace
  //

#ifdef HAVE_MPI
  if(verbose) out << "\nCreating Tpetra::MpiPlatform ...\n";
  const Tpetra::MpiPlatform<Ordinal,Ordinal>  ordinalPlatform(mpiComm);
  const Tpetra::MpiPlatform<Ordinal,Scalar>   scalarPlatform(mpiComm);
#else
  if(verbose) out << "\nCreating Tpetra::SerialPlatform ...\n";
  const Tpetra::SerialPlatform<Ordinal,Ordinal>  ordinalPlatform;
  const Tpetra::SerialPlatform<Ordinal,Scalar>   scalarPlatform;
#endif
  const Tpetra::ElementSpace<Ordinal> elementSpace(globalDim,0,ordinalPlatform);
  // Above: platform gets cloned!
  const Tpetra::VectorSpace<Ordinal,Scalar> vectorSpace(elementSpace,scalarPlatform);
  // Above: platform gets cloned again!
  // Note: The above Tpetra::VectorSpace object is really just a handle to
  // the real underlying object!

  //
  // (B) Create the tridiagonal Tpetra object
  //
  //       [  2  -1             ]
  //       [ -1   2  -1         ]
  //  A =  [      .   .   .     ]
  //       [          -1  2  -1 ]
  //       [             -1   2 ]
  //
  
  // (B.1) Allocate the Tpetra::CisMatrix object.
  RefCountPtr<Tpetra::CisMatrix<Ordinal,Scalar> >
    A_tpetra = rcp(new Tpetra::CisMatrix<Ordinal,Scalar>(vectorSpace));
  // Note that Tpetra::CisMatrix is a handle object but we still use a
  // RefCountPtr to wrap since there can be no confusion when using an RCP to
  // manage an object.

  // (B.2) Get the indexes of the rows on this processor
  const int numMyElements = vectorSpace.getNumMyEntries();
  const std::vector<int> &myGlobalElements = vectorSpace.elementSpace().getMyGlobalElements();

  // (B.3) Fill the local matrix entries one row at a time.
  // Note, we set up Tpetra_Map above to use zero-based indexing and that is what we must use below:
  const Scalar offDiag = -1.0, diag = 2.0*diagScale;
  int numEntries; Scalar values[3]; int indexes[3];
  for( int k = 0; k < numMyElements; ++k ) {
    const int rowIndex = myGlobalElements[k];
    if( rowIndex == 0 ) {                     // First row
      numEntries = 2;
      values[0]  = diag;             values[1]  = offDiag;
      indexes[0] = 0;                indexes[1] = 1; 
    }
    else if( rowIndex == globalDim - 1 ) {    // Last row
      numEntries = 2;
      values[0]  = offDiag;         values[1]  = diag;
      indexes[0] = globalDim-2;     indexes[1] = globalDim-1; 
    }
    else {                                    // Middle rows
      numEntries = 3;
      values[0]  = offDiag;         values[1]  = diag;          values[2]  = offDiag;
      indexes[0] = rowIndex-1;      indexes[1] = rowIndex;      indexes[2] = rowIndex+1; 
    }
    A_tpetra->submitEntries(Tpetra::Insert,rowIndex,numEntries,values,indexes);
  }

  // (B.4) Finish the construction of the Tpetra::CisMatrix
  A_tpetra->fillComplete();

  //
  // (C) Wrap the above created (Tpetra::CisMatrix) Tpetra_Operator object created above
  // into a Thyra::TpetraLinearOp object to turn it into a Thyra::LinearOpBase object
  //

  RefCountPtr<Thyra::LinearOpBase<Scalar> >
    A = rcp(
      new Thyra::TpetraLinearOp<Ordinal,Scalar>(
        Teuchos::rcp_implicit_cast<Tpetra::Operator<Ordinal,Scalar> >(A_tpetra)
        )
      );

  //
  // (D) Finally return the Thyra-wrapped Tpetra matrix object
  //
  
  return A;

  // Note that when this function returns the returned
  // RefCountPtr-wrapped Thyra::LinearOpBase object will own all of the
  // Tpetra objects that went into its construction and these objects
  // will stay around until all of the RefCountPtr objects to the
  // allocated Thyra::LinearOpBase object are removed and destruction
  // occurs!

} // end createExampleTridiagTpetraLinearOp()

#endif // THYRA_CREATE_EXAMPLE_TRIDIAG_TPETRA_LINEAR_OP_HPP
