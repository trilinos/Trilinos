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
//
#ifndef ANASAZI_MULTI_VEC_TRAITS_HPP
#define ANASAZI_MULTI_VEC_TRAITS_HPP

#include "AnasaziMultiVec.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Anasazi {

  template< class TYPE, class MV >
  struct UndefinedMultiVecTraits
  {
    //! This function should not compile if there is an attempt to instantiate!
    static inline TYPE notDefined() { return MV::this_type_is_missing_a_specialization(); };
  };
  
  template<class TYPE, class MV>
  class MultiVecTraits 
  {
  public:
    ///
    static Teuchos::RefCountPtr<MV> Clone( const MV& mv, const int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<MV> CloneCopy( const MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<MV> CloneCopy( const MV& mv, int index[], int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<MV> CloneView( MV& mv, int index[], int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static Teuchos::RefCountPtr<const MV> CloneView( const MV& mv, int index[], int numvecs )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return Teuchos::null; }     
    ///
    static int GetVecLength( const MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return 0; }     
    ///
    static int GetNumberVecs( const MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); return 0; }     
    ///
    static void MvTimesMatAddMv( TYPE alpha, const MV& A, 
				 const Teuchos::SerialDenseMatrix<int,TYPE>& B, TYPE beta, MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvAddMv( TYPE alpha, const MV& A, TYPE beta, const MV& B, MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvTransMv( const MV& mv, TYPE alpha, const MV& A, Teuchos::SerialDenseMatrix<int,TYPE>& B )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvNorm( const MV& mv, TYPE *normvec )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void SetBlock( const MV& A, int index[], int numvecs, MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvRandom( MV& mv )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvInit( MV& mv, TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
    ///
    static void MvPrint( const MV& mv, ostream& os )
    { UndefinedMultiVecTraits<TYPE, MV>::notDefined(); }     
  };


  template<class TYPE>
  class MultiVecTraits<TYPE,MultiVec<TYPE> >
  {
  public:
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > Clone( const MultiVec<TYPE>& mv, const int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).Clone(numvecs) ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > CloneCopy( const MultiVec<TYPE>& mv )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).CloneCopy() ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > CloneCopy( const MultiVec<TYPE>& mv, int index[], int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).CloneCopy(index,numvecs) ); }
    ///
    static Teuchos::RefCountPtr<MultiVec<TYPE> > CloneView( MultiVec<TYPE>& mv, int index[], int numvecs )
    { return Teuchos::rcp( mv.CloneView(index,numvecs) ); }
    ///
    static Teuchos::RefCountPtr<const MultiVec<TYPE> > CloneView( const MultiVec<TYPE>& mv, int index[], int numvecs )
    { return Teuchos::rcp( const_cast<MultiVec<TYPE>&>(mv).CloneView(index,numvecs) ); }
    ///
    static int GetVecLength( const MultiVec<TYPE>& mv )
    { return mv.GetVecLength(); }
    ///
    static int GetNumberVecs( const MultiVec<TYPE>& mv )
    { return mv.GetNumberVecs(); }
    ///
    static void MvTimesMatAddMv( TYPE alpha, const MultiVec<TYPE>& A, 
				 const Teuchos::SerialDenseMatrix<int,TYPE>& B, 
				 TYPE beta, MultiVec<TYPE>& mv )
    { mv.MvTimesMatAddMv(alpha,const_cast<MultiVec<TYPE>&>(A),const_cast<Teuchos::SerialDenseMatrix<int,TYPE>&>(B),beta); }
    ///
    static void MvAddMv( TYPE alpha, const MultiVec<TYPE>& A, TYPE beta, const MultiVec<TYPE>& B, MultiVec<TYPE>& mv )
    { mv.MvAddMv(alpha,const_cast<MultiVec<TYPE>&>(A),beta,const_cast<MultiVec<TYPE>&>(B)); }
    ///
    static void MvTransMv( const MultiVec<TYPE>& mv, TYPE alpha, const MultiVec<TYPE>& A, Teuchos::SerialDenseMatrix<int,TYPE>& B )
    { const_cast<MultiVec<TYPE>&>(mv).MvTransMv(alpha,const_cast<MultiVec<TYPE>&>(A),B); }
    ///
    static void MvNorm( const MultiVec<TYPE>& mv, TYPE *normvec )
    { const_cast<MultiVec<TYPE>&>(mv).MvNorm(normvec); }
    ///
    static void SetBlock( const MultiVec<TYPE>& A, int index[], int numvecs, MultiVec<TYPE>& mv )
    { mv.SetBlock(const_cast<MultiVec<TYPE>&>(A),index,numvecs); }
    ///
    static void MvRandom( MultiVec<TYPE>& mv )
    { mv.MvRandom(); }
    ///
    static void MvInit( MultiVec<TYPE>& mv, TYPE alpha = Teuchos::ScalarTraits<TYPE>::zero() )
    { mv.MvInit(alpha); }
    ///
    static void MvPrint( const MultiVec<TYPE>& mv, ostream& os )
    { const_cast<MultiVec<TYPE>&>(mv).MvPrint(os); }
    
  };
  
} // namespace Anasazi

#endif // ANASAZI_MULTI_VEC_TRAITS_HPP
