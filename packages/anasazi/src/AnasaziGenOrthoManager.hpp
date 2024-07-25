// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziGenOrthoManager.hpp
  \brief  Templated virtual class for providing orthogonalization/orthonormalization methods with matrix-based
          inner products.
*/

#ifndef ANASAZI_GENORTHOMANAGER_HPP
#define ANASAZI_GENORTHOMANAGER_HPP

/*! \class Anasazi::GenOrthoManager
  
  This class provides an interface for orthogonalization managers to provide
  oblique projectors of the form:
  \f[
    P_{X,Y} S = S - X \langle Y, X \rangle^{-1} \langle Y, S \rangle\ .
  \f]
  Such a projector modifies the input in the range on \f$X\f$ in order to 
  make the output orthogonal to the range of \f$Y\f$.
  
  \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "AnasaziMatOrthoManager.hpp"
#include "AnasaziMultiVecTraits.hpp"
#include "AnasaziOperatorTraits.hpp"

namespace Anasazi {

  template <class ScalarType, class MV, class OP>
  class GenOrthoManager : public MatOrthoManager<ScalarType,MV,OP> {
  public:
    //! @name Constructor/Destructor
    //@{ 
    //! Default constructor.
    GenOrthoManager(Teuchos::RCP<const OP> Op = Teuchos::null);

    //! Destructor.
    virtual ~GenOrthoManager() {};
    //@}


    //! @name Orthogonalization methods
    //@{ 

    /*! \brief Applies a series of generic projectors.
     *
     * Given a list of bases <tt>X[i]</tt> and <tt>Y[i]</tt> (a projection pair), this method
     * takes a multivector \c S and applies the projectors
     * \f[
     *    P_{X[i],Y[i]} S = S - X[i] \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     * This operation projects \c S onto the space orthogonal to the <tt>Y[i]</tt>, 
     * along the range of the <tt>X[i]</tt>. The inner product specified by \f$\langle \cdot, 
     * \cdot \rangle\f$ is given by innerProd(). 
     *
     * \note The call 
     * \code 
     * projectGen(S, tuple(X1,X2), tuple(Y1,Y2)) 
     * \endcode 
     * is equivalent to the call
     * \code 
     * projectGen(S, tuple(X2,X1), tuple(Y2,Y1))
     * \endcode
     *
     * The method also returns the coefficients <tt>C[i]</tt> associated with each projection pair, so that 
     * \f[
     *    S_{in} = S_{out} + \sum_i X[i] C[i]
     * \f]
     * and therefore
     * \f[
     *    C[i] = \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     *
     * Lastly, for reasons of efficiency, the user must specify whether the projection pairs are bi-orthonormal with
     * respect to innerProd(), i.e., whether \f$\langle Y[i], X[i] \rangle = I\f$. In the case that the bases are specified
     * to be biorthogonal, the inverse \f$\langle Y, X \rangle^{-1}\f$ will not be computed. Furthermore, the user may optionally 
     * specifiy the image of \c S and the projection pairs under the inner product operator getOp().
     @param S [in/out] The multivector to be modified.<br>
       On output, the columns of \c S will be orthogonal to each <tt>Y[i]</tt>, satisfying
       \f[
          \langle Y[i], S_{out} \rangle = 0
       \f]
       Also, 
       \f[
          S_{in} = S_{out} + \sum_i X[i] C[i]
       \f]

     @param X [in] Multivectors for bases under which \f$S_{in}\f$ is modified.
        
     @param Y [in] Multivectors for bases to which \f$S_{out}\f$ should be orthogonal.

     @param isBiortho [in] A flag specifying whether the bases <tt>X[i]</tt>
            and <tt>Y[i]</tt> are biorthonormal, i.e,. whether \f$\langle Y[i],
            X[i]\rangle == I\f$.

     @param C [out] Coefficients for reconstructing \f$S_{in}\f$ via the bases <tt>X[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c S and <tt>X[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>.<br>
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c S and \c <tt>X[i]</tt>, then a std::invalid_argument
       exception will be thrown.<br>
       Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null pointer,
       the caller will not have access to the computed coefficients <tt>C[i]</tt>.

     @param MS [in/out] If specified by the user, on input \c MS is required to be the image of \c S under the operator getOp(). 
     On output, \c MS will be updated to reflect the changes in \c S.
     
     @param MX [in] If specified by the user, <tt>MX[i]</tt> is required to be the image of <tt>X[i]</tt> under the operator getOp().
     @param MY [in] If specified by the user, <tt>MY[i]</tt> is required to be the image of <tt>Y[i]</tt> under the operator getOp().

     \pre 
     <ul>
     <li>If <tt>X[i] != Teuchos::null</tt> or <tt>Y[i] != Teuchos::null</tt>, then <tt>X[i]</tt> and <tt>Y[i]</tt> are required to 
         have the same number of columns, and each should have the same number of rows as \c S.
     <li>For any <tt>i != j</tt>, \f$\langle Y[i], X[j] \rangle == 0\f$.
     <li>If <tt>biOrtho == true</tt>, \f$\langle Y[i], X[i]\rangle == I\f$
     <li>Otherwise, if <tt>biOrtho == false</tt>, then \f$\langle Y[i], X[i]\rangle\f$ should be Hermitian positive-definite.
     <li>If <tt>X[i]</tt> and <tt>Y[i]</tt> have \f$xc_i\f$ columns and \c S has \f$sc\f$ columns, then <tt>C[i]</tt> if specified must
         be \f$xc_i \times sc\f$.
     </ul>
     
     */
    virtual void projectGen( 
           MV &S,
           Teuchos::Array<Teuchos::RCP<const MV> > X,
           Teuchos::Array<Teuchos::RCP<const MV> > Y,
           bool isBiOrtho,
           Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
                 = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
           Teuchos::RCP<MV> MS                        = Teuchos::null,
           Teuchos::Array<Teuchos::RCP<const MV> > MX = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null)),
           Teuchos::Array<Teuchos::RCP<const MV> > MY = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const = 0;


    /*! \brief Applies a series of generic projectors and returns an orthonormal basis for the residual data.
     *
     * Given a list of bases <tt>X[i]</tt> and <tt>Y[i]</tt> (a projection pair), this method
     * takes a multivector \c S and applies the projectors
     * \f[
     *    P_{X[i],Y[i]} S = S - X[i] \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     * These operation project \c S onto the space orthogonal to the range of the <tt>Y[i]</tt>, 
     * along the range of \c X[i]. The inner product specified by \f$\langle \cdot, \cdot \rangle\f$
     * is given by innerProd(). 
     *
     * The method returns in \c S an orthonormal basis for the residual
     * \f[
     *    \left( \prod_{i} P_{X[i],Y[i]} \right) S_{in} = S_{out} B\ ,
     * \f]
     * where \c B contains the (not necessarily triangular) coefficients of the residual with respect to the 
     * new basis.
     *
     * The method also returns the coefficients <tt>C[i]</tt> and \c B associated with each projection pair, so that 
     * \f[
     *    S_{in} = S_{out} B + \sum_i X[i] C[i]
     * \f]
     * and 
     * \f[
     *    C[i] = \langle Y[i], X[i] \rangle^{-1} \langle Y[i], S \rangle\ .
     * \f]
     *
     * Lastly, for reasons of efficiency, the user must specify whether the projection pairs are bi-orthonormal with
     * respect to innerProd(), i.e., whether \f$\langle Y[i], X[i] \rangle = I\f$. Furthermore, the user may optionally 
     * specifiy the image of \c S and the projection pairs under the inner product operator getOp().
     @param S [in/out] The multivector to be modified.<br>
       On output, the columns of \c S will be orthogonal to each <tt>Y[i]</tt>, satisfying
       \f[
          \langle Y[i], S_{out} \rangle = 0
       \f]
       Also, 
       \f[
          S_{in}(1:m,1:n) = S_{out}(1:m,1:rank) B(1:rank,1:n) + \sum_i X[i] C[i]\ ,
       \f]
       where \c m is the number of rows in \c S, \c n is the number of
       columns in \c S, and \c rank is the value returned from the method.

     @param X [in] Multivectors for bases under which \f$S_{in}\f$ is modified.
        
     @param Y [in] Multivectors for bases to which \f$S_{out}\f$ should be orthogonal.

     @param isBiortho [in] A flag specifying whether the bases <tt>X[i]</tt>
            and <tt>Y[i]</tt> are biorthonormal, i.e,. whether \f$\langle Y[i],
            X[i]\rangle == I\f$.

     @param C [out] Coefficients for reconstructing \f$S_{in}\f$ via the bases <tt>X[i]</tt>. If <tt>C[i]</tt> is a non-null pointer 
       and <tt>C[i]</tt> matches the dimensions of \c X and <tt>Q[i]</tt>, then the coefficients computed during the orthogonalization
       routine will be stored in the matrix <tt>C[i]</tt>.<br>
       If <tt>C[i]</tt> points to a Teuchos::SerialDenseMatrix with size
       inconsistent with \c S and \c <tt>X[i]</tt>, then a std::invalid_argument
       exception will be thrown.<br>
       Otherwise, if <tt>C.size() < i</tt> or <tt>C[i]</tt> is a null pointer,
       the caller will not have access to the computed coefficients <tt>C[i]</tt>.

     @param B [out] The coefficients of the original \c S with respect to the computed basis. If \c B is a non-null pointer and 
     \c B matches the dimensions of \c B, then the
     coefficients computed during the orthogonalization routine will be stored in \c B, similar to calling 
       \code
          innerProd( Sout, Sin, B );
       \endcode
     If \c B points to a Teuchos::SerialDenseMatrix with size inconsistent with
     \c S, then a std::invalid_argument exception will be thrown.<br>
     Otherwise, if \c B is null, the caller will not have access to the computed
     coefficients.<br>

     @param MS [in/out] If specified by the user, on input \c MS is required to be the image of \c S under the operator getOp(). 
     On output, \c MS will be updated to reflect the changes in \c S.
     
     @param MX [in] If specified by the user, <tt>MX[i]</tt> is required to be the image of <tt>X[i]</tt> under the operator getOp().
     @param MY [in] If specified by the user, <tt>MY[i]</tt> is required to be the image of <tt>Y[i]</tt> under the operator getOp().

     \note The matrix \c B is not necessarily triangular (as in a QR
         factorization); see the documentation of specific orthogonalization managers.

     \pre 
     <ul>
     <li>If <tt>X[i] != Teuchos::null</tt> or <tt>Y[i] != Teuchos::null</tt>, then <tt>X[i]</tt> and <tt>Y[i]</tt> are required to 
         have the same number of columns, and each should have the same number of rows as \c S.
     <li>For any <tt>i != j</tt>, \f$\langle Y[i], X[j] \rangle == 0\f$.
     <li>If <tt>biOrtho == true</tt>, \f$\langle Y[i], X[i]\rangle == I\f$
     <li>Otherwise, if <tt>biOrtho == false</tt>, then \f$\langle Y[i], X[i]\rangle\f$ should be Hermitian positive-definite.
     <li>If <tt>X[i]</tt> and <tt>Y[i]</tt> have \f$xc_i\f$ columns and \c S has \f$sc\f$ columns, then <tt>C[i]</tt> if specified must
         be \f$xc_i \times sc\f$.
     <li>If <tt>S</tt> has \f$sc\f$ columns, then \c B if specified must be \f$sc \times sc \f$.
     </ul>

     @return Rank of the basis computed by this method.
     */
    virtual int projectAndNormalizeGen (
           MV &S,
           Teuchos::Array<Teuchos::RCP<const MV> > X,
           Teuchos::Array<Teuchos::RCP<const MV> > Y,
           bool isBiOrtho,
           Teuchos::Array<Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > > C
                 = Teuchos::tuple(Teuchos::RCP< Teuchos::SerialDenseMatrix<int,ScalarType> >(Teuchos::null)),
           Teuchos::RCP<Teuchos::SerialDenseMatrix<int,ScalarType> > B = Teuchos::null,
           Teuchos::RCP<MV> MS                                         = Teuchos::null,
           Teuchos::Array<Teuchos::RCP<const MV> > MX                  = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null)),
           Teuchos::Array<Teuchos::RCP<const MV> > MY                  = Teuchos::tuple(Teuchos::RCP<const MV>(Teuchos::null))
        ) const = 0;

    //@}

  };

  template <class ScalarType,class MV,class OP>
  GenOrthoManager<ScalarType,MV,OP>::GenOrthoManager(Teuchos::RCP<const OP> Op)
    : MatOrthoManager<ScalarType,MV,OP>(Op) {}

} // end of Anasazi namespace


#endif

// end of file AnasaziGenOrthoManager.hpp
