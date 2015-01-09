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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_HELPER_TRAITS_HPP
#define ANASAZI_HELPER_TRAITS_HPP

/*!     \file AnasaziOperatorTraits.hpp
        \brief Virtual base class which defines basic traits for the operator type
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziTypes.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Anasazi {

    /*!  \brief Class which defines basic traits for working with different scalar types.

      An adapter for this traits class must exist for the <tt>ScalarType</tt>.
      If not, this class will produce a compile-time error.

      \ingroup anasazi_opvec_interfaces
    */
    template <class ScalarType>
    class HelperTraits 
    {
        public:

            //! Helper function for correctly storing the Ritz values when the eigenproblem is non-Hermitian
            /*! This allows us to use template specialization to compute the right index vector and correctly
             *  handle complex-conjugate pairs.  
             */
            static void sortRitzValues( 
                    const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& rRV,
                    const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& iRV,
                    std::vector<Value<ScalarType> >* RV, std::vector<int>* RO, std::vector<int>* RI );

            //! Helper function for correctly scaling the eigenvectors of the projected eigenproblem.
            /*! This allows us to use template specialization to compute the right scaling so the
             * Ritz residuals are correct.
             */
            static void scaleRitzVectors( 
                    const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& iRV,
                    Teuchos::SerialDenseMatrix<int, ScalarType>* S );

            //! Helper function for correctly computing the Ritz residuals of the projected eigenproblem.
            /*! This allows us to use template specialization to ensure the Ritz residuals are correct.
            */
            static void computeRitzResiduals( 
                    const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& iRV,
                    const Teuchos::SerialDenseMatrix<int, ScalarType>& S,
                    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>* RR);

    };


    template<class ScalarType>
    void HelperTraits<ScalarType>::sortRitzValues( 
            const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& rRV,
            const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& iRV,
            std::vector<Value<ScalarType> >* RV, std::vector<int>* RO, std::vector<int>* RI )
    {
        typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
        MagnitudeType MT_ZERO = Teuchos::ScalarTraits<MagnitudeType>::zero();

        int curDim = (int)rRV.size();
        int i = 0;

        // Clear the current index.
        RI->clear();

        // Place the Ritz values from rRV and iRV into the RV container.
        while( i < curDim ) {
            if ( iRV[i] != MT_ZERO ) {
                //
                // We will have this situation for real-valued, non-Hermitian matrices.
                (*RV)[i].set(rRV[i], iRV[i]);
                (*RV)[i+1].set(rRV[i+1], iRV[i+1]);

                // Make sure that complex conjugate pairs have their positive imaginary part first.
                if ( (*RV)[i].imagpart < MT_ZERO ) {
                    // The negative imaginary part is first, so swap the order of the ritzValues and ritzOrders.
                    Anasazi::Value<ScalarType> tmp_ritz( (*RV)[i] );
                    (*RV)[i] = (*RV)[i+1];
                    (*RV)[i+1] = tmp_ritz;

                    int tmp_order = (*RO)[i];
                    (*RO)[i] = (*RO)[i+1];
                    (*RO)[i+1] = tmp_order;

                }
                RI->push_back(1); RI->push_back(-1);
                i = i+2;
            } else {
                //
                // The Ritz value is not complex.
                (*RV)[i].set(rRV[i], MT_ZERO);
                RI->push_back(0);
                i++;
            }
        }
    }


    template<class ScalarType>
    void HelperTraits<ScalarType>::scaleRitzVectors( 
            const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& iRV,
            Teuchos::SerialDenseMatrix<int, ScalarType>* S )
    {
        ScalarType ST_ONE = Teuchos::ScalarTraits<ScalarType>::one();

        typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
        MagnitudeType MT_ZERO = Teuchos::ScalarTraits<MagnitudeType>::zero();

        Teuchos::LAPACK<int,MagnitudeType> lapack_mag;
        Teuchos::BLAS<int,ScalarType> blas;

        int i = 0, curDim = S->numRows();
        ScalarType temp;
        ScalarType* s_ptr = S->values();
        while( i < curDim ) {
            if ( iRV[i] != MT_ZERO ) {
                temp = lapack_mag.LAPY2( blas.NRM2( curDim, s_ptr+i*curDim, 1 ), 
                        blas.NRM2( curDim, s_ptr+(i+1)*curDim, 1 ) );
                blas.SCAL( curDim, ST_ONE/temp, s_ptr+i*curDim, 1 );
                blas.SCAL( curDim, ST_ONE/temp, s_ptr+(i+1)*curDim, 1 );
                i = i+2;
            } else {
                temp = blas.NRM2( curDim, s_ptr+i*curDim, 1 );
                blas.SCAL( curDim, ST_ONE/temp, s_ptr+i*curDim, 1 );
                i++;
            }
        }
    }

    template<class ScalarType>
    void HelperTraits<ScalarType>::computeRitzResiduals( 
            const std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& iRV,
            const Teuchos::SerialDenseMatrix<int, ScalarType>& S,
            std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>* RR )
    {
        typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
        MagnitudeType MT_ZERO = Teuchos::ScalarTraits<MagnitudeType>::zero();

        Teuchos::LAPACK<int,MagnitudeType> lapack_mag;
        Teuchos::BLAS<int,ScalarType> blas;

        int i = 0;
        int s_stride = S.stride();
        int s_rows = S.numRows();
        int s_cols = S.numCols();
        ScalarType* s_ptr = S.values();

        while( i < s_cols ) {
            if ( iRV[i] != MT_ZERO ) {
                (*RR)[i] = lapack_mag.LAPY2( blas.NRM2(s_rows, s_ptr + i*s_stride, 1),
                                             blas.NRM2(s_rows, s_ptr + (i+1)*s_stride, 1) );
                (*RR)[i+1] = (*RR)[i];
                i = i+2;
            } else {
                (*RR)[i] = blas.NRM2(s_rows, s_ptr + i*s_stride, 1);
                i++;
            }
        }          
    }

#ifdef HAVE_TEUCHOS_COMPLEX
    // Partial template specializations for the complex scalar type.

    /*!  \brief Class which defines basic traits for working with different scalar types.

      An adapter for this traits class must exist for the <tt>ScalarType</tt>.
      If not, this class will produce a compile-time error.

      \ingroup anasazi_opvec_interfaces
    */
    template <class T>
    class HelperTraits<ANSZI_CPLX_CLASS<T> >
    {
        public:
          static void sortRitzValues( 
              const std::vector<T>& rRV, 
              const std::vector<T>& iRV,
              std::vector<Value<ANSZI_CPLX_CLASS<T> > >* RV, 
              std::vector<int>* RO, std::vector<int>* RI );

            static void scaleRitzVectors( 
                const std::vector<T>& iRV,
                Teuchos::SerialDenseMatrix<int, ANSZI_CPLX_CLASS<T> >* S );

            static void computeRitzResiduals( 
                const std::vector<T>& iRV,
                const Teuchos::SerialDenseMatrix<int, ANSZI_CPLX_CLASS<T> >& S,
                std::vector<T>* RR );
    };

    template<class T>
    void HelperTraits<ANSZI_CPLX_CLASS<T> >::sortRitzValues( 
            const std::vector<T>& rRV, 
            const std::vector<T>& iRV,
            std::vector<Value<ANSZI_CPLX_CLASS<T> > >* RV, 
            std::vector<int>* RO, std::vector<int>* RI )
    {
        (void)RO;
        int curDim = (int)rRV.size();
        int i = 0;

        // Clear the current index.
        RI->clear();

        // Place the Ritz values from rRV and iRV into the RV container.
        while( i < curDim ) {
            (*RV)[i].set(rRV[i], iRV[i]);
            RI->push_back(0);
            i++;
        }    
    }

    template<class T>
    void HelperTraits<ANSZI_CPLX_CLASS<T> >::scaleRitzVectors( 
            const std::vector<T>& iRV,
            Teuchos::SerialDenseMatrix<int, ANSZI_CPLX_CLASS<T> >* S )
    {
      (void)iRV;
      typedef ANSZI_CPLX_CLASS<T> ST;
      ST ST_ONE = Teuchos::ScalarTraits<ST>::one();

      Teuchos::BLAS<int,ST> blas;

      int i = 0, curDim = S->numRows();
      ST temp;
      ST* s_ptr = S->values();
      while( i < curDim ) {
        temp = blas.NRM2( curDim, s_ptr+i*curDim, 1 );
        blas.SCAL( curDim, ST_ONE/temp, s_ptr+i*curDim, 1 );
        i++;
      }
    }

    template<class T>
    void HelperTraits<ANSZI_CPLX_CLASS<T> >::computeRitzResiduals( 
            const std::vector<T>& iRV,
            const Teuchos::SerialDenseMatrix<int, ANSZI_CPLX_CLASS<T> >& S,
            std::vector<T>* RR )
    {
        (void)iRV;
        Teuchos::BLAS<int,ANSZI_CPLX_CLASS<T> > blas;

        int s_stride = S.stride();
        int s_rows = S.numRows();
        int s_cols = S.numCols();
        ANSZI_CPLX_CLASS<T>* s_ptr = S.values();

        for (int i=0; i<s_cols; ++i ) {
            (*RR)[i] = blas.NRM2(s_rows, s_ptr + i*s_stride, 1);
        }
    }          
#endif

} // end Anasazi namespace


#endif // ANASAZI_HELPER_TRAITS_HPP
