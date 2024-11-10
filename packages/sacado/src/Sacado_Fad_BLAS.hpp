// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_BLAS_HPP
#define SACADO_FAD_BLAS_HPP

#include "Teuchos_BLAS.hpp"
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_dummy_arg.hpp"

namespace Sacado {

  namespace Fad {

    template <typename OrdinalType, typename FadType>
    class ArrayTraits {

      typedef typename Sacado::ValueType<FadType>::type ValueType;
      typedef typename Sacado::ScalarType<FadType>::type scalar_type;
      typedef typename Sacado::dummy<ValueType,scalar_type>::type ScalarType;
      
    public:
      
      ArrayTraits(bool use_dynamic = true,
		  OrdinalType workspace_size = 0);

      ArrayTraits(const ArrayTraits& a);

      ~ArrayTraits();
      
      void unpack(const FadType& a, OrdinalType& n_dot, ValueType& val, 
		  const ValueType*& dot) const;
      
      void unpack(const FadType* a, OrdinalType n, OrdinalType inc,
		  OrdinalType& n_dot, OrdinalType& inc_val, 
		  OrdinalType& inc_dot,
		  const ValueType*& val, const ValueType*& dot) const;
      
      void unpack(const FadType* A, OrdinalType m, OrdinalType n, 
		  OrdinalType lda, OrdinalType& n_dot, 
		  OrdinalType& lda_val, OrdinalType& lda_dot,
		  const ValueType*& val, const ValueType*& dot) const;

      void unpack(const ValueType& a, OrdinalType& n_dot, ValueType& val, 
		  const ValueType*& dot) const;
      
      void unpack(const ValueType* a, OrdinalType n, OrdinalType inc,
		  OrdinalType& n_dot, OrdinalType& inc_val, 
		  OrdinalType& inc_dot,
		  const ValueType*& val, const ValueType*& dot) const;
      
      void unpack(const ValueType* A, OrdinalType m, OrdinalType n, 
		  OrdinalType lda, OrdinalType& n_dot, 
		  OrdinalType& lda_val, OrdinalType& lda_dot,
		  const ValueType*& val, const ValueType*& dot) const;

      void unpack(const ScalarType& a, OrdinalType& n_dot, ScalarType& val, 
		  const ScalarType*& dot) const;
      
      void unpack(const ScalarType* a, OrdinalType n, OrdinalType inc,
		  OrdinalType& n_dot, OrdinalType& inc_val, 
		  OrdinalType& inc_dot,
		  const ScalarType*& val, const ScalarType*& dot) const;
      
      void unpack(const ScalarType* A, OrdinalType m, OrdinalType n, 
		  OrdinalType lda, OrdinalType& n_dot, 
		  OrdinalType& lda_val, OrdinalType& lda_dot,
		  const ScalarType*& val, const ScalarType*& dot) const;

      void unpack(FadType& a, OrdinalType& n_dot, OrdinalType& final_n_dot, 
		  ValueType& val, ValueType*& dot) const;
      
      void unpack(FadType* a, OrdinalType n, OrdinalType inc,
		  OrdinalType& n_dot, OrdinalType& final_n_dot, 
		  OrdinalType& inc_val, OrdinalType& inc_dot,
		  ValueType*& val, ValueType*& dot) const;
      
      void unpack(FadType* A, OrdinalType m, OrdinalType n, OrdinalType lda, 
		  OrdinalType& n_dot, OrdinalType& final_n_dot, 
		  OrdinalType& lda_val, OrdinalType& lda_dot,
		  ValueType*& val, ValueType*& dot) const;

      void pack(FadType& a, OrdinalType n_dot, const ValueType& val, 
		const ValueType* dot) const;
      
      void pack(FadType* a, OrdinalType n, OrdinalType inc,
		OrdinalType n_dot, OrdinalType inc_val, OrdinalType inc_dot,
		const ValueType* val, const ValueType* dot) const;
      
      void pack(FadType* A, OrdinalType m, OrdinalType n, 
		OrdinalType lda, OrdinalType n_dot, 
		OrdinalType lda_val, OrdinalType lda_dot,
		const ValueType* val, const ValueType* dot) const;

      void free(const FadType& a, OrdinalType n_dot, 
		const ValueType* dot) const;
      
      void free(const FadType* a, OrdinalType n, OrdinalType n_dot,
		OrdinalType inc_val, OrdinalType inc_dot,
		const ValueType* val, const ValueType* dot) const;

      void free(const FadType* A, OrdinalType m, OrdinalType n, 
		OrdinalType n_dot, OrdinalType lda_val, OrdinalType lda_dot,
		const ValueType* val, const ValueType* dot) const;

      void free(const ValueType& a, OrdinalType n_dot, 
		const ValueType* dot) const {}
      
      void free(const ValueType* a, OrdinalType n, OrdinalType n_dot,
		OrdinalType inc_val, OrdinalType inc_dot,
		const ValueType* val, const ValueType* dot) const {}

      void free(const ValueType* A, OrdinalType m, OrdinalType n, 
		OrdinalType n_dot, OrdinalType lda_val, OrdinalType lda_dot,
		const ValueType* val, const ValueType* dot) const {}

      void free(const ScalarType& a, OrdinalType n_dot, 
		const ScalarType* dot) const {}
      
      void free(const ScalarType* a, OrdinalType n, OrdinalType n_dot,
		OrdinalType inc_val, OrdinalType inc_dot,
		const ScalarType* val, const ScalarType* dot) const {}

      void free(const ScalarType* A, OrdinalType m, OrdinalType n, 
		OrdinalType n_dot, OrdinalType lda_val, OrdinalType lda_dot,
		const ScalarType* val, const ScalarType* dot) const {}

      ValueType* allocate_array(OrdinalType size) const;

      void free_array(const ValueType* ptr, OrdinalType size) const;

      bool is_array_contiguous(const FadType* a, OrdinalType n, 
			       OrdinalType n_dot) const;

    protected:

      //! Use dynamic memory allocation
      bool use_dynamic;

      //! Size of static workspace
      OrdinalType workspace_size;

      //! Workspace for holding contiguous values/derivatives
      mutable ValueType *workspace;

      //! Pointer to current free entry in workspace
      mutable ValueType *workspace_pointer;
		
    };

    template <typename T> struct ArrayValueType { typedef T type; };

    //! Fad specializations for Teuchos::BLAS wrappers
    template <typename OrdinalType, typename FadType>
    class BLAS : public Teuchos::DefaultBLASImpl<OrdinalType,FadType> {    
      
      typedef typename Teuchos::ScalarTraits<FadType>::magnitudeType MagnitudeType;
      typedef typename Sacado::ValueType<FadType>::type ValueType;
      typedef typename Sacado::ScalarType<FadType>::type scalar_type;
      typedef typename Sacado::dummy<ValueType,scalar_type>::type ScalarType;
      typedef Teuchos::DefaultBLASImpl<OrdinalType,FadType> BLASType;
    
    public:
      //! @name Constructor/Destructor.
      //@{ 
    
      //! Default constructor.
      BLAS(bool use_default_impl = true,
	   bool use_dynamic = true, OrdinalType static_workspace_size = 0);

      //! Copy constructor.

      BLAS(const BLAS& x);

      //! Destructor.
      virtual ~BLAS();

      //@}
      
      //! @name Level 1 BLAS Routines.
      //@{ 
      
      //! Computes a Givens plane rotation.
      void ROTG(FadType* da, FadType* db, MagnitudeType* c, FadType* s) const { 
	BLASType::ROTG(da,db,c,s); 
      }
      
      //! Applies a Givens plane rotation.
      void ROT(const OrdinalType n, FadType* dx, const OrdinalType incx, 
	       FadType* dy, const OrdinalType incy, MagnitudeType* c, 
	       FadType* s) const { 
	BLASType::ROT(n,dx,incx,dy,incy,c,s); 
      }
      
      //! Scale the std::vector \c x by the constant \c alpha.
      void SCAL(const OrdinalType n, const FadType& alpha, FadType* x, 
		const OrdinalType incx) const;

      //! Copy the std::vector \c x to the std::vector \c y.
      void COPY(const OrdinalType n, const FadType* x, 
		const OrdinalType incx, FadType* y, 
		const OrdinalType incy) const;

      //! Perform the operation: \c y \c <- \c y+alpha*x.
      template <typename alpha_type, typename x_type>
      void AXPY(const OrdinalType n, const alpha_type& alpha, 
		const x_type* x, const OrdinalType incx, FadType* y, 
		const OrdinalType incy) const;

      //! Sum the absolute values of the entries of \c x.
      typename Teuchos::ScalarTraits<FadType>::magnitudeType 
      ASUM(const OrdinalType n, const FadType* x, 
	   const OrdinalType incx) const {
	return BLASType::ASUM(n,x,incx);
      }

      //! Form the dot product of the vectors \c x and \c y.
      template <typename x_type, typename y_type>
      FadType DOT(const OrdinalType n, const x_type* x, 
		  const OrdinalType incx, const y_type* y, 
		  const OrdinalType incy) const;

      //! Compute the 2-norm of the std::vector \c x.
      MagnitudeType NRM2(const OrdinalType n, const FadType* x, 
			 const OrdinalType incx) const;

      //! Return the index of the element of \c x with the maximum magnitude.
      OrdinalType IAMAX(const OrdinalType n, const FadType* x, 
			const OrdinalType incx) const {
	return BLASType::IAMAX(n,x,incx); 
      }

      //@}
      
      //! @name Level 2 BLAS Routines.
      //@{ 
      
      /*! 
       * \brief Performs the matrix-std::vector operation:  
       * \c y \c <- \c alpha*A*x+beta*y or \c y \c <- \c alpha*A'*x+beta*y 
       * where \c A is a general \c m by \c n matrix.
       */
      template <typename alpha_type, typename A_type, typename x_type,
		typename beta_type>
      void GEMV(Teuchos::ETransp trans, const OrdinalType m, 
		const OrdinalType n, 
		const alpha_type& alpha, const A_type* A, 
		const OrdinalType lda, const x_type* x, 
		const OrdinalType incx, const beta_type& beta, 
		FadType* y, const OrdinalType incy) const;

      /*!
       * \brief Performs the matrix-std::vector operation:  
       * \c x \c <- \c A*x or \c x \c <- \c A'*x where \c A is a unit/non-unit 
       * \c n by \c n upper/lower triangular matrix.
       */
      template <typename A_type>
      void TRMV(Teuchos::EUplo uplo, Teuchos::ETransp trans, 
		Teuchos::EDiag diag, const OrdinalType n, 
		const A_type* A, const OrdinalType lda, FadType* x, 
		const OrdinalType incx) const;

      //! Performs the rank 1 operation:  \c A \c <- \c alpha*x*y'+A. 
      template <typename alpha_type, typename x_type, typename y_type>
      void GER(const OrdinalType m, const OrdinalType n, 
	       const alpha_type& alpha, 
	       const x_type* x, const OrdinalType incx, 
	       const y_type* y, const OrdinalType incy, 
	       FadType* A, const OrdinalType lda) const;
      
      //@}
      
      //! @name Level 3 BLAS Routines. 
      //@{ 
      
      /*! 
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*op(B)+beta*C where \c op(A) is either \c A 
       * or \c A', \c op(B) is either \c B or \c B', and C is an \c m by \c k 
       * matrix.
       */
      template <typename alpha_type, typename A_type, typename B_type,
		typename beta_type>
      void GEMM(Teuchos::ETransp transa, Teuchos::ETransp transb, 
		const OrdinalType m, const OrdinalType n, const OrdinalType k, 
		const alpha_type& alpha, const A_type* A, const OrdinalType lda,
		const B_type* B, const OrdinalType ldb, const beta_type& beta, 
		FadType* C, const OrdinalType ldc) const;
      
      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*A*B+beta*C or \c C \c <- \c alpha*B*A+beta*C where 
       * \c A is an \c m by \c m or \c n by \c n symmetric matrix and \c B is a 
       * general matrix.
       */
      template <typename alpha_type, typename A_type, typename B_type,
		typename beta_type>
      void SYMM(Teuchos::ESide side, Teuchos::EUplo uplo, const OrdinalType m, 
		const OrdinalType n, 
		const alpha_type& alpha, const A_type* A, 
		const OrdinalType lda, const B_type* B, 
		const OrdinalType ldb,
		const beta_type& beta, FadType* C, 
		const OrdinalType ldc) const;
      
      /*!
       * \brief Performs the matrix-matrix operation: 
       * \c C \c <- \c alpha*op(A)*B+beta*C or 
       * \c C \c <- \c alpha*B*op(A)+beta*C where \c op(A) is an unit/non-unit, 
       * upper/lower triangular matrix and \c B is a general matrix.
       */
      template <typename alpha_type, typename A_type>
      void TRMM(Teuchos::ESide side, Teuchos::EUplo uplo, 
		Teuchos::ETransp transa, Teuchos::EDiag diag, 
		const OrdinalType m, const OrdinalType n, 
		const alpha_type& alpha, 
		const A_type* A, const OrdinalType lda, 
		FadType* B, const OrdinalType ldb) const;

      /*! 
       * \brief Solves the matrix equations:  
       * \c op(A)*X=alpha*B or \c X*op(A)=alpha*B where \c X and \c B are \c m 
       * by \c n matrices, \c A is a unit/non-unit, upper/lower triangular 
       * matrix and \c op(A) is \c A or \c A'.  The matrix \c X is overwritten 
       * on \c B.
       */
      template <typename alpha_type, typename A_type>
      void TRSM(Teuchos::ESide side, Teuchos::EUplo uplo, 
		Teuchos::ETransp transa, Teuchos::EDiag diag, 
		const OrdinalType m, const OrdinalType n, 
		const alpha_type& alpha, 
		const A_type* A, const OrdinalType lda, 
		FadType* B, const OrdinalType ldb) const;

      //@}

    protected:

      //! ArrayTraits for packing/unpacking value/derivative arrays
      ArrayTraits<OrdinalType,FadType> arrayTraits;

      //! BLAS for values
      Teuchos::BLAS<OrdinalType, ValueType> blas;

      //! Use custom or default implementation
      bool use_default_impl;
      
      //! Temporary array for GEMV
      mutable std::vector<ValueType> gemv_Ax;

      //! Temporary array for GEMM
      mutable std::vector<ValueType> gemm_AB;

    protected:

      //! Implementation of DOT
      template <typename x_type, typename y_type>
      void Fad_DOT(const OrdinalType n,
		   const x_type* x, 
		   const OrdinalType incx, 
		   const OrdinalType n_x_dot,  
		   const x_type* x_dot, 
		   const OrdinalType incx_dot,
		   const y_type* y, 
		   const OrdinalType incy, 
		   const OrdinalType n_y_dot, 
		   const y_type* y_dot,
		   const OrdinalType incy_dot,
		   ValueType& z,
		   const OrdinalType n_z_dot,
		   ValueType* zdot) const;

      //! Implementation of GEMV
      template <typename alpha_type, typename A_type, typename x_type,
		typename beta_type>
      void Fad_GEMV(Teuchos::ETransp trans, 
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const alpha_type& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const alpha_type* alpha_dot,
		    const A_type* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const A_type* A_dot,
		    const OrdinalType lda_dot,
		    const x_type* x, 
		    const OrdinalType incx, 
		    const OrdinalType n_x_dot, 
		    const x_type* x_dot, 
		    const OrdinalType incx_dot, 
		    const beta_type& beta, 
		    const OrdinalType n_beta_dot, 
		    const beta_type* beta_dot,
		    ValueType* y, 
		    const OrdinalType incy, 
		    const OrdinalType n_y_dot, 
		    ValueType* y_dot,
		    const OrdinalType incy_dot,
		    const OrdinalType n_dot) const;

      //! Implementation of GER
      template <typename alpha_type, typename x_type, typename y_type>
      void Fad_GER(const OrdinalType m, 
		   const OrdinalType n, 
		   const alpha_type& alpha, 
		   const OrdinalType n_alpha_dot, 
		   const alpha_type* alpha_dot,
		   const x_type* x, 
		   const OrdinalType incx, 
		   const OrdinalType n_x_dot, 
		   const x_type* x_dot, 
		   const OrdinalType incx_dot, 
		   const y_type* y, 
		   const OrdinalType incy, 
		   const OrdinalType n_y_dot, 
		   const y_type* y_dot,
		   const OrdinalType incy_dot,
		   ValueType* A, 
		   const OrdinalType lda, 
		   const OrdinalType n_A_dot, 
		   ValueType* A_dot,
		   const OrdinalType lda_dot,
		   const OrdinalType n_dot) const;

      //! Implementation of GEMM
      template <typename alpha_type, typename A_type, typename B_type,
		typename beta_type>
      void Fad_GEMM(Teuchos::ETransp transa,
		    Teuchos::ETransp transb,
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const OrdinalType k,
		    const alpha_type& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const alpha_type* alpha_dot,
		    const A_type* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const A_type* A_dot,
		    const OrdinalType lda_dot,
		    const B_type* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    const B_type* B_dot, 
		    const OrdinalType ldb_dot, 
		    const beta_type& beta, 
		    const OrdinalType n_beta_dot, 
		    const beta_type* beta_dot,
		    ValueType* C, 
		    const OrdinalType ldc, 
		    const OrdinalType n_C_dot, 
		    ValueType* C_dot,
		    const OrdinalType ldc_dot,
		    const OrdinalType n_dot) const;

      //! Implementation of SYMM
      template <typename alpha_type, typename A_type, typename B_type,
		typename beta_type>
      void Fad_SYMM(Teuchos::ESide side, 
		    Teuchos::EUplo uplo,
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const alpha_type& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const alpha_type* alpha_dot,
		    const A_type* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const A_type* A_dot,
		    const OrdinalType lda_dot,
		    const B_type* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    const B_type* B_dot, 
		    const OrdinalType ldb_dot, 
		    const beta_type& beta, 
		    const OrdinalType n_beta_dot, 
		    const beta_type* beta_dot,
		    ValueType* C, 
		    const OrdinalType ldc, 
		    const OrdinalType n_C_dot, 
		    ValueType* C_dot,
		    const OrdinalType ldc_dot,
		    const OrdinalType n_dot) const;

      //! Implementation of TRMM
      template <typename alpha_type, typename A_type>
      void Fad_TRMM(Teuchos::ESide side, 
		    Teuchos::EUplo uplo,
		    Teuchos::ETransp transa, 
		    Teuchos::EDiag diag, 
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const alpha_type& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const alpha_type* alpha_dot,
		    const A_type* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const A_type* A_dot,
		    const OrdinalType lda_dot,
		    ValueType* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    ValueType* B_dot, 
		    const OrdinalType ldb_dot, 
		    const OrdinalType n_dot) const;

      //! Implementation of TRMM
      template <typename alpha_type, typename A_type>
      void Fad_TRSM(Teuchos::ESide side, 
		    Teuchos::EUplo uplo,
		    Teuchos::ETransp transa, 
		    Teuchos::EDiag diag, 
		    const OrdinalType m, 
		    const OrdinalType n, 
		    const alpha_type& alpha, 
		    const OrdinalType n_alpha_dot, 
		    const alpha_type* alpha_dot,
		    const A_type* A, 
		    const OrdinalType lda, 
		    const OrdinalType n_A_dot, 
		    const A_type* A_dot,
		    const OrdinalType lda_dot,
		    ValueType* B, 
		    const OrdinalType ldb, 
		    const OrdinalType n_B_dot, 
		    ValueType* B_dot, 
		    const OrdinalType ldb_dot, 
		    const OrdinalType n_dot) const;

    }; // class FadBLAS

  }  // namespace Fad

  // template <typename FadType> ArrayValueType<FadType> { typedef ValueType type; };
  // template <> ArrayValueType<ValueType> { typedef ValueType type; };
  // template <> ArrayValueType<ScalarType> { typedef ScalarType type; };

} // namespace Sacado

// Here we provide partial specializations for Teuchos::BLAS for each Fad type
#define TEUCHOS_BLAS_FAD_SPEC(FADTYPE)					\
namespace Teuchos {							\
  template <typename OrdinalType, typename ValueT>			\
  class BLAS< OrdinalType, FADTYPE<ValueT> > :				\
    public Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT> > {		\
  public:								\
    BLAS(bool use_default_impl = true,	bool use_dynamic = true,	\
	 OrdinalType static_workspace_size = 0) :			\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT> >(		\
	use_default_impl, use_dynamic,static_workspace_size) {}		\
    BLAS(const BLAS& x) :						\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT> >(x) {}		\
    virtual ~BLAS() {}							\
  };									\
}									\
namespace Sacado {							\
  namespace Fad {							\
    template <typename ValueT>						\
    struct ArrayValueType< FADTYPE<ValueT> > {				\
      typedef ValueT type;						\
    };									\
  }									\
}
#define TEUCHOS_BLAS_SFAD_SPEC(FADTYPE)					\
namespace Teuchos {							\
  template <typename OrdinalType, typename ValueT, int Num>		\
  class BLAS< OrdinalType, FADTYPE<ValueT,Num> > :			\
    public Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,Num> > {	\
  public:								\
    BLAS(bool use_default_impl = true,	bool use_dynamic = true,	\
	 OrdinalType static_workspace_size = 0) :			\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,Num> >(		\
	use_default_impl, use_dynamic, static_workspace_size) {}	\
    BLAS(const BLAS& x) :						\
      Sacado::Fad::BLAS< OrdinalType, FADTYPE<ValueT,Num> >(x) {}	\
    virtual ~BLAS() {}							\
  };									\
}									\
namespace Sacado {							\
  namespace Fad {							\
    template <typename ValueT, int Num>					\
    struct ArrayValueType< FADTYPE<ValueT,Num> > {			\
      typedef ValueT type;						\
    };									\
  }									\
}
TEUCHOS_BLAS_FAD_SPEC(Sacado::Fad::DFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::Fad::SFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::Fad::SLFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::Fad::DVFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::ELRFad::DFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::ELRFad::SFad)
TEUCHOS_BLAS_SFAD_SPEC(Sacado::ELRFad::SLFad)
TEUCHOS_BLAS_FAD_SPEC(Sacado::CacheFad::DFad)

#undef TEUCHOS_BLAS_FAD_SPEC
#undef TEUCHOS_BLAS_SFAD_SPEC

#include "Sacado_Fad_BLASImp.hpp"

#endif // SACADO_FAD_BLAS_HPP 
