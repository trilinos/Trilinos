// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_Assert.hpp"

template <typename OrdinalType, typename FadType>
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
ArrayTraits(bool use_dynamic_,
            OrdinalType workspace_size_) :
  use_dynamic(use_dynamic_),
  workspace_size(workspace_size_),
  workspace(NULL),
  workspace_pointer(NULL)
{
  if (workspace_size > 0) {
    workspace = new ValueType[workspace_size];
    workspace_pointer = workspace;
  }
}

template <typename OrdinalType, typename FadType>
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
ArrayTraits(const ArrayTraits& a) :
  use_dynamic(a.use_dynamic),
  workspace_size(a.workspace_size),
  workspace(NULL),
  workspace_pointer(NULL)
{
  if (workspace_size > 0) {
    workspace = new ValueType*[workspace_size];
    workspace_pointer = workspace;
  }
}


template <typename OrdinalType, typename FadType>
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
~ArrayTraits()
{
// #ifdef SACADO_DEBUG
//   TEUCHOS_TEST_FOR_EXCEPTION(workspace_pointer != workspace,
//                   std::logic_error,
//                   "ArrayTraits::~ArrayTraits(): " <<
//                   "Destructor called with non-zero used workspace. " <<
//                   "Currently used size is " << workspace_pointer-workspace <<
//                   ".");

// #endif

  if (workspace_size > 0)
    delete [] workspace;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const FadType& a, OrdinalType& n_dot, ValueType& val,
       const ValueType*& dot) const
{
  n_dot = a.size();
  val = a.val();
  if (n_dot > 0)
    dot = &a.fastAccessDx(0);
  else
    dot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const FadType* a, OrdinalType n, OrdinalType inc,
       OrdinalType& n_dot, OrdinalType& inc_val, OrdinalType& inc_dot,
       const ValueType*& cval, const ValueType*& cdot) const
{
  cdot = NULL; // Ensure dot is always initialized

  if (n == 0) {
    n_dot = 0;
    inc_val = 0;
    inc_dot = 0;
    cval = NULL;
    cdot = NULL;
    return;
  }

  n_dot = a[0].size();
  bool is_contiguous = is_array_contiguous(a, n, n_dot);
  if (is_contiguous) {
    inc_val = inc;
    inc_dot = inc;
    cval = &a[0].val();
    if (n_dot > 0)
      cdot = &a[0].fastAccessDx(0);
  }
  else {
    inc_val = 1;
    inc_dot = 0;
    ValueType *val = allocate_array(n);
    ValueType *dot = NULL;
    if (n_dot > 0) {
      inc_dot = 1;
      dot = allocate_array(n*n_dot);
    }
    for (OrdinalType i=0; i<n; i++) {
      val[i] = a[i*inc].val();
      for (OrdinalType j=0; j<n_dot; j++)
        dot[j*n+i] = a[i*inc].fastAccessDx(j);
    }

    cval = val;
    cdot = dot;
  }
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const FadType* A, OrdinalType m, OrdinalType n, OrdinalType lda,
       OrdinalType& n_dot, OrdinalType& lda_val, OrdinalType& lda_dot,
       const ValueType*& cval, const ValueType*& cdot) const
{
  cdot = NULL; // Ensure dot is always initialized

  if (m*n == 0) {
    n_dot = 0;
    lda_val = 0;
    lda_dot = 0;
    cval = NULL;
    cdot = NULL;
    return;
  }

  n_dot = A[0].size();
  bool is_contiguous = is_array_contiguous(A, m*n, n_dot);
  if (is_contiguous) {
    lda_val = lda;
    lda_dot = lda;
    cval = &A[0].val();
    if (n_dot > 0)
      cdot = &A[0].fastAccessDx(0);
  }
  else {
    lda_val = m;
    lda_dot = 0;
    ValueType *val = allocate_array(m*n);
    ValueType *dot = NULL;
    if (n_dot > 0) {
      lda_dot = m;
      dot = allocate_array(m*n*n_dot);
    }
    for (OrdinalType j=0; j<n; j++) {
      for (OrdinalType i=0; i<m; i++) {
        val[j*m+i] = A[j*lda+i].val();
        for (OrdinalType k=0; k<n_dot; k++)
          dot[(k*n+j)*m+i] = A[j*lda+i].fastAccessDx(k);
      }
    }

    cval = val;
    cdot = dot;
  }
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const ValueType& a, OrdinalType& n_dot, ValueType& val,
       const ValueType*& dot) const
{
  n_dot = 0;
  val = a;
  dot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const ValueType* a, OrdinalType n, OrdinalType inc,
       OrdinalType& n_dot, OrdinalType& inc_val, OrdinalType& inc_dot,
       const ValueType*& cval, const ValueType*& cdot) const
{
  n_dot = 0;
  inc_val = inc;
  inc_dot = 0;
  cval = a;
  cdot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const ValueType* A, OrdinalType m, OrdinalType n, OrdinalType lda,
       OrdinalType& n_dot, OrdinalType& lda_val, OrdinalType& lda_dot,
       const ValueType*& cval, const ValueType*& cdot) const
{
  n_dot = 0;
  lda_val = lda;
  lda_dot = 0;
  cval = A;
  cdot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const ScalarType& a, OrdinalType& n_dot, ScalarType& val,
       const ScalarType*& dot) const
{
  n_dot = 0;
  val = a;
  dot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const ScalarType* a, OrdinalType n, OrdinalType inc,
       OrdinalType& n_dot, OrdinalType& inc_val, OrdinalType& inc_dot,
       const ScalarType*& cval, const ScalarType*& cdot) const
{
  n_dot = 0;
  inc_val = inc;
  inc_dot = 0;
  cval = a;
  cdot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(const ScalarType* A, OrdinalType m, OrdinalType n, OrdinalType lda,
       OrdinalType& n_dot, OrdinalType& lda_val, OrdinalType& lda_dot,
       const ScalarType*& cval, const ScalarType*& cdot) const
{
  n_dot = 0;
  lda_val = lda;
  lda_dot = 0;
  cval = A;
  cdot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(FadType& a, OrdinalType& n_dot, OrdinalType& final_n_dot, ValueType& val,
       ValueType*& dot) const
{
  n_dot = a.size();
  val = a.val();
#ifdef SACADO_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(n_dot > 0 && final_n_dot > 0 && final_n_dot != n_dot,
                     std::logic_error,
                     "ArrayTraits::unpack(): FadType has wrong number of " <<
                     "derivative components.  Got " << n_dot <<
                     ", expected " << final_n_dot << ".");
#endif
  if (n_dot > final_n_dot)
    final_n_dot = n_dot;

  OrdinalType n_avail = a.availableSize();
  if (n_avail < final_n_dot) {
    dot = alloate(final_n_dot);
    for (OrdinalType i=0; i<final_n_dot; i++)
        dot[i] = 0.0;
  }
  else if (n_avail > 0)
    dot = &a.fastAccessDx(0);
  else
    dot = NULL;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(FadType* a, OrdinalType n, OrdinalType inc, OrdinalType& n_dot,
       OrdinalType& final_n_dot, OrdinalType& inc_val, OrdinalType& inc_dot,
       ValueType*& val, ValueType*& dot) const
{
  if (n == 0) {
    inc_val = 0;
    inc_dot = 0;
    val = NULL;
    dot = NULL;
    return;
  }

  n_dot = a[0].size();
  bool is_contiguous = is_array_contiguous(a, n, n_dot);
#ifdef SACADO_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(n_dot > 0 && final_n_dot > 0 && final_n_dot != n_dot,
                     std::logic_error,
                     "ArrayTraits::unpack(): FadType has wrong number of " <<
                     "derivative components.  Got " << n_dot <<
                     ", expected " << final_n_dot << ".");
#endif
  if (n_dot > final_n_dot)
    final_n_dot = n_dot;

  if (is_contiguous) {
    inc_val = inc;
    val = &a[0].val();
  }
  else {
    inc_val = 1;
    val = allocate_array(n);
    for (OrdinalType i=0; i<n; i++)
      val[i] = a[i*inc].val();
  }

  OrdinalType n_avail = a[0].availableSize();
  if (is_contiguous && n_avail >= final_n_dot && final_n_dot > 0) {
    inc_dot = inc;
    dot = &a[0].fastAccessDx(0);
  }
  else if (final_n_dot > 0) {
    inc_dot = 1;
    dot = allocate_array(n*final_n_dot);
    for (OrdinalType i=0; i<n; i++) {
      if (n_dot > 0)
        for (OrdinalType j=0; j<n_dot; j++)
          dot[j*n+i] = a[i*inc].fastAccessDx(j);
      else
        for (OrdinalType j=0; j<final_n_dot; j++)
          dot[j*n+i] = 0.0;
    }
  }
  else {
    inc_dot = 0;
    dot = NULL;
  }
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
unpack(FadType* A, OrdinalType m, OrdinalType n, OrdinalType lda,
       OrdinalType& n_dot, OrdinalType& final_n_dot,
       OrdinalType& lda_val, OrdinalType& lda_dot,
       ValueType*& val, ValueType*& dot) const
{
  if (m*n == 0) {
    lda_val = 0;
    lda_dot = 0;
    val = NULL;
    dot = NULL;
    return;
  }

  n_dot = A[0].size();
  bool is_contiguous = is_array_contiguous(A, m*n, n_dot);
#ifdef SACADO_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(n_dot > 0 && final_n_dot > 0 && final_n_dot != n_dot,
                     std::logic_error,
                     "ArrayTraits::unpack(): FadType has wrong number of " <<
                     "derivative components.  Got " << n_dot <<
                     ", expected " << final_n_dot << ".");
#endif
  if (n_dot > final_n_dot)
    final_n_dot = n_dot;

  if (is_contiguous) {
    lda_val = lda;
    val = &A[0].val();
  }
  else {
    lda_val = m;
    val = allocate_array(m*n);
    for (OrdinalType j=0; j<n; j++)
      for (OrdinalType i=0; i<m; i++)
        val[j*m+i] = A[j*lda+i].val();
  }

  OrdinalType n_avail = A[0].availableSize();
  if (is_contiguous && n_avail >= final_n_dot && final_n_dot > 0) {
    lda_dot = lda;
    dot = &A[0].fastAccessDx(0);
  }
  else if (final_n_dot > 0) {
    lda_dot = m;
    dot = allocate_array(m*n*final_n_dot);
    for (OrdinalType j=0; j<n; j++) {
      for (OrdinalType i=0; i<m; i++) {
        if (n_dot > 0)
          for (OrdinalType k=0; k<n_dot; k++)
            dot[(k*n+j)*m+i] = A[j*lda+i].fastAccessDx(k);
        else
          for (OrdinalType k=0; k<final_n_dot; k++)
            dot[(k*n+j)*m+i] = 0.0;
      }
    }
  }
  else {
    lda_dot = 0;
    dot = NULL;
  }
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
pack(FadType& a, OrdinalType n_dot, const ValueType& val,
     const ValueType* dot) const
{
  a.val() = val;

  if (n_dot == 0)
    return;

  if (a.size() != n_dot)
    a.resize(n_dot);
  if (a.dx() != dot)
    for (OrdinalType i=0; i<n_dot; i++)
      a.fastAccessDx(i) = dot[i];
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
pack(FadType* a, OrdinalType n, OrdinalType inc,
     OrdinalType n_dot, OrdinalType inc_val, OrdinalType inc_dot,
     const ValueType* val, const ValueType* dot) const
{
  if (n == 0)
    return;

  // Copy values
  if (&a[0].val() != val)
    for (OrdinalType i=0; i<n; i++)
      a[i*inc].val() = val[i*inc_val];

  if (n_dot == 0)
    return;

  // Resize derivative arrays
  if (a[0].size() != n_dot)
    for (OrdinalType i=0; i<n; i++)
      a[i*inc].resize(n_dot);

  // Copy derivatives
  if (a[0].dx() != dot)
    for (OrdinalType i=0; i<n; i++)
      for (OrdinalType j=0; j<n_dot; j++)
        a[i*inc].fastAccessDx(j) = dot[(i+j*n)*inc_dot];
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
pack(FadType* A, OrdinalType m, OrdinalType n, OrdinalType lda,
     OrdinalType n_dot, OrdinalType lda_val, OrdinalType lda_dot,
     const ValueType* val, const ValueType* dot) const
{
  if (m*n == 0)
    return;

  // Copy values
  if (&A[0].val() != val)
    for (OrdinalType j=0; j<n; j++)
      for (OrdinalType i=0; i<m; i++)
        A[i+j*lda].val() = val[i+j*lda_val];

  if (n_dot == 0)
    return;

  // Resize derivative arrays
  if (A[0].size() != n_dot)
    for (OrdinalType j=0; j<n; j++)
      for (OrdinalType i=0; i<m; i++)
        A[i+j*lda].resize(n_dot);

  // Copy derivatives
  if (A[0].dx() != dot)
    for (OrdinalType j=0; j<n; j++)
      for (OrdinalType i=0; i<m; i++)
        for (OrdinalType k=0; k<n_dot; k++)
          A[i+j*lda].fastAccessDx(k) = dot[i+(j+k*n)*lda_dot];
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
free(const FadType& a, OrdinalType n_dot, const ValueType* dot) const
{
  if (n_dot > 0 && a.dx() != dot) {
    free_array(dot, n_dot);
  }
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
free(const FadType* a, OrdinalType n, OrdinalType n_dot,
     OrdinalType inc_val, OrdinalType inc_dot,
     const ValueType* val, const ValueType* dot) const
{
  if (n == 0)
    return;

  if (val != &a[0].val())
    free_array(val, n*inc_val);

  if (n_dot > 0  && a[0].dx() != dot)
    free_array(dot, n*inc_dot*n_dot);
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
free(const FadType* A, OrdinalType m, OrdinalType n, OrdinalType n_dot,
     OrdinalType lda_val, OrdinalType lda_dot,
     const ValueType* val, const ValueType* dot) const
{
  if (m*n == 0)
    return;

  if (val != &A[0].val())
    free_array(val, lda_val*n);

  if (n_dot > 0  && A[0].dx() != dot)
    free_array(dot, lda_dot*n*n_dot);
}

template <typename OrdinalType, typename FadType>
typename Sacado::Fad::ArrayTraits<OrdinalType,FadType>::ValueType*
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
allocate_array(OrdinalType size) const
{
  if (use_dynamic)
    return new ValueType[size];

#ifdef SACADO_DEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(workspace_pointer + size - workspace > workspace_size,
                     std::logic_error,
                     "ArrayTraits::allocate_array(): " <<
                     "Requested workspace memory beyond size allocated. " <<
                     "Workspace size is " << workspace_size <<
                     ", currently used is " << workspace_pointer-workspace <<
                     ", requested size is " << size << ".");

#endif

  ValueType *v = workspace_pointer;
  workspace_pointer += size;
  return v;
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
free_array(const ValueType* ptr, OrdinalType size) const
{
  if (use_dynamic && ptr != NULL)
    delete [] ptr;
  else
    workspace_pointer -= size;
}

template <typename OrdinalType, typename FadType>
bool
Sacado::Fad::ArrayTraits<OrdinalType,FadType>::
is_array_contiguous(const FadType* a, OrdinalType n, OrdinalType n_dot) const
{
  return (n > 0) &&
    (&(a[n-1].val())-&(a[0].val()) == n-1) &&
    (a[n-1].dx()-a[0].dx() == n-1);
}

template <typename OrdinalType, typename FadType>
Sacado::Fad::BLAS<OrdinalType,FadType>::
BLAS(bool use_default_impl_,
     bool use_dynamic_, OrdinalType static_workspace_size_) :
  BLASType(),
  arrayTraits(use_dynamic_, static_workspace_size_),
  blas(),
  use_default_impl(use_default_impl_)
{
}

template <typename OrdinalType, typename FadType>
Sacado::Fad::BLAS<OrdinalType,FadType>::
BLAS(const BLAS& x) :
  BLASType(x),
  arrayTraits(x.arrayTraits),
  blas(x.blas),
  use_default_impl(x.use_default_impl)
{
}

template <typename OrdinalType, typename FadType>
Sacado::Fad::BLAS<OrdinalType,FadType>::
~BLAS()
{
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
SCAL(const OrdinalType n, const FadType& alpha, FadType* x,
     const OrdinalType incx) const
{
  if (use_default_impl) {
    BLASType::SCAL(n,alpha,x,incx);
    return;
  }

  // Unpack input values & derivatives
  ValueType alpha_val;
  const ValueType *alpha_dot;
  ValueType *x_val, *x_dot;
  OrdinalType n_alpha_dot = 0, n_x_dot = 0, n_dot = 0;
  OrdinalType incx_val, incx_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  n_dot = n_alpha_dot;
  arrayTraits.unpack(x, n, incx, n_x_dot, n_dot, incx_val, incx_dot,
                     x_val, x_dot);

#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_x_dot != n_dot && n_x_dot != 0),
                     std::logic_error,
                     "BLAS::SCAL(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif

  // Call differentiated routine
  if (n_x_dot > 0)
    blas.SCAL(n*n_x_dot, alpha_val, x_dot, incx_dot);
  for (OrdinalType i=0; i<n_alpha_dot; i++)
    blas.AXPY(n, alpha_dot[i], x_val, incx_val, x_dot+i*n*incx_dot, incx_dot);
  blas.SCAL(n, alpha_val, x_val, incx_val);

  // Pack values and derivatives for result
  arrayTraits.pack(x, n, incx, n_dot, incx_val, incx_dot, x_val, x_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(x, n, n_dot, incx_val, incx_dot, x_val, x_dot);
}

template <typename OrdinalType, typename FadType>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
COPY(const OrdinalType n, const FadType* x, const OrdinalType incx,
     FadType* y, const OrdinalType incy) const
{
  if (use_default_impl) {
    BLASType::COPY(n,x,incx,y,incy);
    return;
  }

  if (n == 0)
    return;

  OrdinalType n_x_dot = x[0].size();
  OrdinalType n_y_dot = y[0].size();
  if (n_x_dot == 0 || n_y_dot == 0 || n_x_dot != n_y_dot ||
      !arrayTraits.is_array_contiguous(x, n, n_x_dot) ||
      !arrayTraits.is_array_contiguous(y, n, n_y_dot))
    BLASType::COPY(n,x,incx,y,incy);
  else {
    blas.COPY(n, &x[0].val(), incx, &y[0].val(), incy);
    blas.COPY(n*n_x_dot, &x[0].fastAccessDx(0), incx, &y[0].fastAccessDx(0),
              incy);
  }
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename x_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
AXPY(const OrdinalType n, const alpha_type& alpha, const x_type* x,
     const OrdinalType incx, FadType* y, const OrdinalType incy) const
{
  if (use_default_impl) {
    BLASType::AXPY(n,alpha,x,incx,y,incy);
    return;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  const typename ArrayValueType<x_type>::type *x_val, *x_dot;
  ValueType *y_val, *y_dot;
  OrdinalType n_alpha_dot, n_x_dot, n_y_dot, n_dot;
  OrdinalType incx_val, incy_val, incx_dot, incy_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(x, n, incx, n_x_dot, incx_val, incx_dot, x_val, x_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_x_dot > 0)
    n_dot = n_x_dot;

  // Unpack and allocate y
  arrayTraits.unpack(y, n, incy, n_y_dot, n_dot, incy_val, incy_dot, y_val,
                     y_dot);

#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_x_dot != n_dot && n_x_dot != 0) ||
                     (n_y_dot != n_dot && n_y_dot != 0),
                     std::logic_error,
                     "BLAS::AXPY(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif

  // Call differentiated routine
  if (n_x_dot > 0)
    blas.AXPY(n*n_x_dot, alpha_val, x_dot, incx_dot, y_dot, incy_dot);
  for (OrdinalType i=0; i<n_alpha_dot; i++)
    blas.AXPY(n, alpha_dot[i], x_val, incx_val, y_dot+i*n*incy_dot, incy_dot);
  blas.AXPY(n, alpha_val, x_val, incx_val, y_val, incy_val);

  // Pack values and derivatives for result
  arrayTraits.pack(y, n, incy, n_dot, incy_val, incy_dot, y_val, y_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(x, n, n_x_dot, incx_val, incx_dot, x_val, x_dot);
  arrayTraits.free(y, n, n_dot, incy_val, incy_dot, y_val, y_dot);
}

template <typename OrdinalType, typename FadType>
template <typename x_type, typename y_type>
FadType
Sacado::Fad::BLAS<OrdinalType,FadType>::
DOT(const OrdinalType n, const x_type* x, const OrdinalType incx,
    const y_type* y, const OrdinalType incy) const
{
  if (use_default_impl)
    return BLASType::DOT(n,x,incx,y,incy);

  // Unpack input values & derivatives
  const typename ArrayValueType<x_type>::type *x_val, *x_dot;
  const typename ArrayValueType<y_type>::type *y_val, *y_dot;
  OrdinalType n_x_dot, n_y_dot;
  OrdinalType incx_val, incy_val, incx_dot, incy_dot;
  arrayTraits.unpack(x, n, incx, n_x_dot, incx_val, incx_dot, x_val, x_dot);
  arrayTraits.unpack(y, n, incy, n_y_dot, incy_val, incy_dot, y_val, y_dot);

  // Compute size
  OrdinalType n_z_dot = 0;
  if (n_x_dot > 0)
    n_z_dot = n_x_dot;
  else if (n_y_dot > 0)
    n_z_dot = n_y_dot;

  // Unpack and allocate z
  FadType z(n_z_dot, 0.0);
  ValueType& z_val = z.val();
  ValueType *z_dot = &z.fastAccessDx(0);

  // Call differentiated routine
  Fad_DOT(n, x_val, incx_val, n_x_dot, x_dot, incx_dot,
          y_val, incy_val, n_y_dot, y_dot, incy_dot,
          z_val, n_z_dot, z_dot);

  // Free temporary arrays
  arrayTraits.free(x, n, n_x_dot, incx_val, incx_dot, x_val, x_dot);
  arrayTraits.free(y, n, n_y_dot, incy_val, incy_dot, y_val, y_dot);

  return z;
}

template <typename OrdinalType, typename FadType>
typename Sacado::Fad::BLAS<OrdinalType,FadType>::MagnitudeType
Sacado::Fad::BLAS<OrdinalType,FadType>::
NRM2(const OrdinalType n, const FadType* x, const OrdinalType incx) const
{
  if (use_default_impl)
    return BLASType::NRM2(n,x,incx);

  // Unpack input values & derivatives
  const ValueType *x_val, *x_dot;
  OrdinalType n_x_dot, incx_val, incx_dot;
  arrayTraits.unpack(x, n, incx, n_x_dot, incx_val, incx_dot, x_val, x_dot);

  // Unpack and allocate z
  MagnitudeType z(n_x_dot, 0.0);

  // Call differentiated routine
  z.val() = blas.NRM2(n, x_val, incx_val);
  // if (!Teuchos::ScalarTraits<FadType>::isComplex && incx_dot == 1)
  //   blas.GEMV(Teuchos::TRANS, n, n_x_dot, 1.0/z.val(), x_dot, n, x_val,
  //          incx_val, 1.0, &z.fastAccessDx(0), OrdinalType(1));
  // else
    for (OrdinalType i=0; i<n_x_dot; i++)
      z.fastAccessDx(i) =
        Teuchos::ScalarTraits<ValueType>::magnitude(blas.DOT(n, x_dot+i*n*incx_dot, incx_dot, x_val, incx_val)) / z.val();

  // Free temporary arrays
  arrayTraits.free(x, n, n_x_dot, incx_val, incx_dot, x_val, x_dot);

  return z;
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type, typename x_type,
          typename beta_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
GEMV(Teuchos::ETransp trans, const OrdinalType m, const OrdinalType n,
     const alpha_type& alpha, const A_type* A,
     const OrdinalType lda, const x_type* x,
     const OrdinalType incx, const beta_type& beta,
     FadType* y, const OrdinalType incy) const
{
  if (use_default_impl) {
    BLASType::GEMV(trans,m,n,alpha,A,lda,x,incx,beta,y,incy);
    return;
  }

  OrdinalType n_x_rows = n;
  OrdinalType n_y_rows = m;
  if (trans != Teuchos::NO_TRANS) {
    n_x_rows = m;
    n_y_rows = n;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  typename ArrayValueType<beta_type>::type beta_val;
  const typename ArrayValueType<beta_type>::type *beta_dot;
  const typename ArrayValueType<A_type>::type *A_val = 0, *A_dot = 0;
  const typename ArrayValueType<x_type>::type *x_val = 0, *x_dot = 0;
  ValueType *y_val, *y_dot;
  OrdinalType n_alpha_dot, n_A_dot, n_x_dot, n_beta_dot, n_y_dot = 0, n_dot;
  OrdinalType lda_val, incx_val, incy_val, lda_dot, incx_dot, incy_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(A, m, n, lda, n_A_dot, lda_val, lda_dot, A_val, A_dot);
  arrayTraits.unpack(x, n_x_rows, incx, n_x_dot, incx_val, incx_dot, x_val,
                     x_dot);
  arrayTraits.unpack(beta, n_beta_dot, beta_val, beta_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_A_dot > 0)
    n_dot = n_A_dot;
  else if (n_x_dot > 0)
    n_dot = n_x_dot;
  else if (n_beta_dot > 0)
    n_dot = n_beta_dot;

  // Unpack and allocate y
  arrayTraits.unpack(y, n_y_rows, incy, n_y_dot, n_dot, incy_val, incy_dot,
                     y_val, y_dot);

  // Call differentiated routine
  Fad_GEMV(trans, m, n, alpha_val, n_alpha_dot, alpha_dot, A_val, lda_val,
           n_A_dot, A_dot, lda_dot, x_val, incx_val, n_x_dot, x_dot, incx_dot,
           beta_val, n_beta_dot, beta_dot, y_val, incy_val, n_y_dot, y_dot,
           incy_dot, n_dot);

  // Pack values and derivatives for result
  arrayTraits.pack(y, n_y_rows, incy, n_dot, incy_val, incy_dot, y_val, y_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(A, m, n, n_A_dot, lda_val, lda_dot, A_val, A_dot);
  arrayTraits.free(x, n_x_rows, n_x_dot, incx_val, incx_dot, x_val, x_dot);
  arrayTraits.free(beta, n_beta_dot, beta_dot);
  arrayTraits.free(y, n_y_rows, n_dot, incy_val, incy_dot, y_val, y_dot);
}

template <typename OrdinalType, typename FadType>
template <typename A_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
TRMV(Teuchos::EUplo uplo, Teuchos::ETransp trans, Teuchos::EDiag diag,
     const OrdinalType n, const A_type* A, const OrdinalType lda,
     FadType* x, const OrdinalType incx) const
{
  if (use_default_impl) {
    BLASType::TRMV(uplo,trans,diag,n,A,lda,x,incx);
    return;
  }

  // Unpack input values & derivatives
  const typename ArrayValueType<A_type>::type *A_val, *A_dot;
  ValueType *x_val, *x_dot;
  OrdinalType n_A_dot = 0, n_x_dot = 0, n_dot = 0;
  OrdinalType lda_val, incx_val, lda_dot, incx_dot;
  arrayTraits.unpack(A, n, n, lda, n_A_dot, lda_val, lda_dot, A_val, A_dot);
  n_dot = n_A_dot;
  arrayTraits.unpack(x, n, incx, n_x_dot, n_dot, incx_val, incx_dot, x_val,
                     x_dot);

#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_x_dot != n_dot && n_x_dot != 0),
                     std::logic_error,
                     "BLAS::TRMV(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif

  // Compute [xd_1 .. xd_n] = A*[xd_1 .. xd_n]
  if (n_x_dot > 0) {
    if (incx_dot == 1)
      blas.TRMM(Teuchos::LEFT_SIDE, uplo, trans, diag, n, n_x_dot, 1.0, A_val,
                lda_val, x_dot, n);
    else
      for (OrdinalType i=0; i<n_x_dot; i++)
        blas.TRMV(uplo, trans, diag, n, A_val, lda_val, x_dot+i*incx_dot*n,
                  incx_dot);
  }

  // Compute [xd_1 .. xd_n] = [Ad_1*x .. Ad_n*x]
  if (gemv_Ax.size() != std::size_t(n))
      gemv_Ax.resize(n);
  for (OrdinalType i=0; i<n_A_dot; i++) {
    blas.COPY(n, x_val, incx_val, &gemv_Ax[0], OrdinalType(1));
    blas.TRMV(uplo, trans, Teuchos::NON_UNIT_DIAG, n, A_dot+i*lda_dot*n,
              lda_dot, &gemv_Ax[0], OrdinalType(1));
    blas.AXPY(n, 1.0, &gemv_Ax[0], OrdinalType(1), x_dot+i*incx_dot*n,
              incx_dot);
  }

  // Compute x = A*x
  blas.TRMV(uplo, trans, diag, n, A_val, lda_val, x_val, incx_val);

  // Pack values and derivatives for result
  arrayTraits.pack(x, n, incx, n_dot, incx_val, incx_dot, x_val, x_dot);

  // Free temporary arrays
  arrayTraits.free(A, n, n, n_A_dot, lda_val, lda_dot, A_val, A_dot);
  arrayTraits.free(x, n, n_dot, incx_val, incx_dot, x_val, x_dot);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename x_type, typename y_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
GER(const OrdinalType m, const OrdinalType n, const alpha_type& alpha,
    const x_type* x, const OrdinalType incx,
    const y_type* y, const OrdinalType incy,
    FadType* A, const OrdinalType lda) const
{
  if (use_default_impl) {
    BLASType::GER(m,n,alpha,x,incx,y,incy,A,lda);
    return;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  const typename ArrayValueType<x_type>::type *x_val = 0, *x_dot = 0;
  const typename ArrayValueType<y_type>::type *y_val = 0, *y_dot = 0;
  ValueType *A_val, *A_dot;
  OrdinalType n_alpha_dot, n_x_dot, n_y_dot, n_A_dot, n_dot;
  OrdinalType lda_val, incx_val, incy_val, lda_dot, incx_dot, incy_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(x, m, incx, n_x_dot, incx_val, incx_dot, x_val, x_dot);
  arrayTraits.unpack(y, n, incy, n_y_dot, incy_val, incy_dot, y_val, y_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_x_dot > 0)
    n_dot = n_x_dot;
  else if (n_y_dot > 0)
    n_dot = n_y_dot;

  // Unpack and allocate A
  arrayTraits.unpack(A, m, n, lda, n_A_dot, n_dot, lda_val, lda_dot, A_val,
                     A_dot);

  // Call differentiated routine
  Fad_GER(m, n, alpha_val, n_alpha_dot, alpha_dot, x_val, incx_val,
          n_x_dot, x_dot, incx_dot, y_val, incy_val, n_y_dot, y_dot,
          incy_dot, A_val, lda_val, n_A_dot, A_dot, lda_dot, n_dot);

  // Pack values and derivatives for result
  arrayTraits.pack(A, m, n, lda, n_dot, lda_val, lda_dot, A_val, A_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(x, m, n_x_dot, incx_val, incx_dot, x_val, x_dot);
  arrayTraits.free(y, n, n_y_dot, incy_val, incy_dot, y_val, y_dot);
  arrayTraits.free(A, m, n, n_dot, lda_val, lda_dot, A_val, A_dot);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type, typename B_type,
          typename beta_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
GEMM(Teuchos::ETransp transa, Teuchos::ETransp transb,
     const OrdinalType m, const OrdinalType n, const OrdinalType k,
     const alpha_type& alpha, const A_type* A, const OrdinalType lda,
     const B_type* B, const OrdinalType ldb, const beta_type& beta,
     FadType* C, const OrdinalType ldc) const
{
  if (use_default_impl) {
    BLASType::GEMM(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
    return;
  }

  OrdinalType n_A_rows = m;
  OrdinalType n_A_cols = k;
  if (transa != Teuchos::NO_TRANS) {
    n_A_rows = k;
    n_A_cols = m;
  }

  OrdinalType n_B_rows = k;
  OrdinalType n_B_cols = n;
  if (transb != Teuchos::NO_TRANS) {
    n_B_rows = n;
    n_B_cols = k;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  typename ArrayValueType<beta_type>::type beta_val;
  const typename ArrayValueType<beta_type>::type *beta_dot;
  const typename ArrayValueType<A_type>::type *A_val, *A_dot;
  const typename ArrayValueType<B_type>::type *B_val, *B_dot;
  ValueType *C_val, *C_dot;
  OrdinalType n_alpha_dot, n_A_dot, n_B_dot, n_beta_dot, n_C_dot = 0, n_dot;
  OrdinalType lda_val, ldb_val, ldc_val, lda_dot, ldb_dot, ldc_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(A, n_A_rows, n_A_cols, lda, n_A_dot, lda_val, lda_dot,
                     A_val, A_dot);
  arrayTraits.unpack(B, n_B_rows, n_B_cols, ldb, n_B_dot, ldb_val, ldb_dot,
                     B_val, B_dot);
  arrayTraits.unpack(beta, n_beta_dot, beta_val, beta_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_A_dot > 0)
    n_dot = n_A_dot;
  else if (n_B_dot > 0)
    n_dot = n_B_dot;
  else if (n_beta_dot > 0)
    n_dot = n_beta_dot;

  // Unpack and allocate C
  arrayTraits.unpack(C, m, n, ldc, n_C_dot, n_dot, ldc_val, ldc_dot, C_val,
                     C_dot);

  // Call differentiated routine
  Fad_GEMM(transa, transb, m, n, k,
           alpha_val, n_alpha_dot, alpha_dot,
           A_val, lda_val, n_A_dot, A_dot, lda_dot,
           B_val, ldb_val, n_B_dot, B_dot, ldb_dot,
           beta_val, n_beta_dot, beta_dot,
           C_val, ldc_val, n_C_dot, C_dot, ldc_dot, n_dot);

  // Pack values and derivatives for result
  arrayTraits.pack(C, m, n, ldc, n_dot, ldc_val, ldc_dot, C_val, C_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(A, n_A_rows, n_A_cols, n_A_dot, lda_val, lda_dot, A_val,
                   A_dot);
  arrayTraits.free(B, n_B_rows, n_B_cols, n_B_dot, ldb_val, ldb_dot, B_val,
                   B_dot);
  arrayTraits.free(beta, n_beta_dot, beta_dot);
  arrayTraits.free(C, m, n, n_dot, ldc_val, ldc_dot, C_val, C_dot);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type, typename B_type,
          typename beta_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
SYMM(Teuchos::ESide side, Teuchos::EUplo uplo,
     const OrdinalType m, const OrdinalType n,
     const alpha_type& alpha, const A_type* A, const OrdinalType lda,
     const B_type* B, const OrdinalType ldb, const beta_type& beta,
     FadType* C, const OrdinalType ldc) const
{
  if (use_default_impl) {
    BLASType::SYMM(side,uplo,m,n,alpha,A,lda,B,ldb,beta,C,ldc);
    return;
  }

  OrdinalType n_A_rows = m;
  OrdinalType n_A_cols = m;
  if (side == Teuchos::RIGHT_SIDE) {
    n_A_rows = n;
    n_A_cols = n;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  typename ArrayValueType<beta_type>::type beta_val;
  const typename ArrayValueType<beta_type>::type *beta_dot;
  const typename ArrayValueType<A_type>::type *A_val, *A_dot;
  const typename ArrayValueType<B_type>::type *B_val, *B_dot;
  ValueType *C_val, *C_dot;
  OrdinalType n_alpha_dot, n_A_dot, n_B_dot, n_beta_dot, n_C_dot, n_dot;
  OrdinalType lda_val, ldb_val, ldc_val, lda_dot, ldb_dot, ldc_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(A, n_A_rows, n_A_cols, lda, n_A_dot, lda_val, lda_dot,
                     A_val, A_dot);
  arrayTraits.unpack(B, m, n, ldb, n_B_dot, ldb_val, ldb_dot, B_val, B_dot);
  arrayTraits.unpack(beta, n_beta_dot, beta_val, beta_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_A_dot > 0)
    n_dot = n_A_dot;
  else if (n_B_dot > 0)
    n_dot = n_B_dot;
  else if (n_beta_dot > 0)
    n_dot = n_beta_dot;

  // Unpack and allocate C
  arrayTraits.unpack(C, m, n, ldc, n_C_dot, n_dot, ldc_val, ldc_dot, C_val,
                     C_dot);

  // Call differentiated routine
  Fad_SYMM(side, uplo, m, n,
           alpha_val, n_alpha_dot, alpha_dot,
           A_val, lda_val, n_A_dot, A_dot, lda_dot,
           B_val, ldb_val, n_B_dot, B_dot, ldb_dot,
           beta_val, n_beta_dot, beta_dot,
           C_val, ldc_val, n_C_dot, C_dot, ldc_dot, n_dot);

  // Pack values and derivatives for result
  arrayTraits.pack(C, m, n, ldc, n_dot, ldc_val, ldc_dot, C_val, C_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(A, n_A_rows, n_A_cols, n_A_dot, lda_val, lda_dot, A_val,
                   A_dot);
  arrayTraits.free(B, m, n, n_B_dot, ldb_val, ldb_dot, B_val, B_dot);
  arrayTraits.free(beta, n_beta_dot, beta_dot);
  arrayTraits.free(C, m, n, n_dot, ldc_val, ldc_dot, C_val, C_dot);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
TRMM(Teuchos::ESide side, Teuchos::EUplo uplo,
     Teuchos::ETransp transa, Teuchos::EDiag diag,
     const OrdinalType m, const OrdinalType n,
     const alpha_type& alpha, const A_type* A, const OrdinalType lda,
     FadType* B, const OrdinalType ldb) const
{
  if (use_default_impl) {
    BLASType::TRMM(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
    return;
  }

  OrdinalType n_A_rows = m;
  OrdinalType n_A_cols = m;
  if (side == Teuchos::RIGHT_SIDE) {
    n_A_rows = n;
    n_A_cols = n;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  const typename ArrayValueType<A_type>::type *A_val, *A_dot;
  ValueType *B_val, *B_dot;
  OrdinalType n_alpha_dot, n_A_dot, n_B_dot, n_dot;
  OrdinalType lda_val, ldb_val, lda_dot, ldb_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(A, n_A_rows, n_A_cols, lda, n_A_dot, lda_val, lda_dot,
                     A_val, A_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_A_dot > 0)
    n_dot = n_A_dot;

  // Unpack and allocate B
  arrayTraits.unpack(B, m, n, ldb, n_B_dot, n_dot, ldb_val, ldb_dot, B_val,
                     B_dot);

  // Call differentiated routine
  Fad_TRMM(side, uplo, transa, diag, m, n,
           alpha_val, n_alpha_dot, alpha_dot,
           A_val, lda_val, n_A_dot, A_dot, lda_dot,
           B_val, ldb_val, n_B_dot, B_dot, ldb_dot, n_dot);

  // Pack values and derivatives for result
  arrayTraits.pack(B, m, n, ldb, n_dot, ldb_val, ldb_dot, B_val, B_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(A, n_A_rows, n_A_cols, n_A_dot, lda_val, lda_dot, A_val,
                   A_dot);
  arrayTraits.free(B, m, n, n_dot, ldb_val, ldb_dot, B_val, B_dot);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
TRSM(Teuchos::ESide side, Teuchos::EUplo uplo,
     Teuchos::ETransp transa, Teuchos::EDiag diag,
     const OrdinalType m, const OrdinalType n,
     const alpha_type& alpha, const A_type* A, const OrdinalType lda,
     FadType* B, const OrdinalType ldb) const
{
  if (use_default_impl) {
    BLASType::TRSM(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb);
    return;
  }

  OrdinalType n_A_rows = m;
  OrdinalType n_A_cols = m;
  if (side == Teuchos::RIGHT_SIDE) {
    n_A_rows = n;
    n_A_cols = n;
  }

  // Unpack input values & derivatives
  typename ArrayValueType<alpha_type>::type alpha_val;
  const typename ArrayValueType<alpha_type>::type *alpha_dot;
  const typename ArrayValueType<A_type>::type *A_val, *A_dot;
  ValueType *B_val, *B_dot;
  OrdinalType n_alpha_dot, n_A_dot, n_B_dot, n_dot;
  OrdinalType lda_val, ldb_val, lda_dot, ldb_dot;
  arrayTraits.unpack(alpha, n_alpha_dot, alpha_val, alpha_dot);
  arrayTraits.unpack(A, n_A_rows, n_A_cols, lda, n_A_dot, lda_val, lda_dot,
                     A_val, A_dot);

  // Compute size
  n_dot = 0;
  if (n_alpha_dot > 0)
    n_dot = n_alpha_dot;
  else if (n_A_dot > 0)
    n_dot = n_A_dot;

  // Unpack and allocate B
  arrayTraits.unpack(B, m, n, ldb, n_B_dot, n_dot, ldb_val, ldb_dot, B_val,
                     B_dot);

  // Call differentiated routine
  Fad_TRSM(side, uplo, transa, diag, m, n,
           alpha_val, n_alpha_dot, alpha_dot,
           A_val, lda_val, n_A_dot, A_dot, lda_dot,
           B_val, ldb_val, n_B_dot, B_dot, ldb_dot, n_dot);

  // Pack values and derivatives for result
  arrayTraits.pack(B, m, n, ldb, n_dot, ldb_val, ldb_dot, B_val, B_dot);

  // Free temporary arrays
  arrayTraits.free(alpha, n_alpha_dot, alpha_dot);
  arrayTraits.free(A, n_A_rows, n_A_cols, n_A_dot, lda_val, lda_dot, A_val,
                   A_dot);
  arrayTraits.free(B, m, n, n_dot, ldb_val, ldb_dot, B_val, B_dot);
}

template <typename OrdinalType, typename FadType>
template <typename x_type, typename y_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_DOT(const OrdinalType n,
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
        ValueType* z_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_x_dot != n_z_dot && n_x_dot != 0) ||
                     (n_y_dot != n_z_dot && n_y_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_DOT(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif

  // Compute [zd_1 .. zd_n] = [xd_1 .. xd_n]^T*y
  if (n_x_dot > 0) {
    if (incx_dot == OrdinalType(1))
      blas.GEMV(Teuchos::TRANS, n, n_x_dot, 1.0, x_dot, n, y, incy, 0.0, z_dot,
                OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_z_dot; i++)
        z_dot[i] = blas.DOT(n, x_dot+i*incx_dot*n, incx_dot, y, incy);
  }

  // Compute [zd_1 .. zd_n] += [yd_1 .. yd_n]^T*x
  if (n_y_dot > 0) {
    if (incy_dot == OrdinalType(1) &&
        !Teuchos::ScalarTraits<ValueType>::isComplex)
      blas.GEMV(Teuchos::TRANS, n, n_y_dot, 1.0, y_dot, n, x, incx, 1.0, z_dot,
                OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_z_dot; i++)
        z_dot[i] += blas.DOT(n, x, incx, y_dot+i*incy_dot*n, incy_dot);
  }

  // Compute z = x^T*y
  z = blas.DOT(n, x, incx, y, incy);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type, typename x_type,
          typename beta_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_GEMV(Teuchos::ETransp trans,
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
         const OrdinalType n_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_x_dot != n_dot && n_x_dot != 0) ||
                     (n_beta_dot != n_dot && n_beta_dot != 0) ||
                     (n_y_dot != n_dot && n_y_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_GEMV(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif
  OrdinalType n_A_rows = m;
  OrdinalType n_A_cols = n;
  OrdinalType n_x_rows = n;
  OrdinalType n_y_rows = m;
  if (trans == Teuchos::TRANS) {
    n_A_rows = n;
    n_A_cols = m;
    n_x_rows = m;
    n_y_rows = n;
  }

  // Compute [yd_1 .. yd_n] = beta*[yd_1 .. yd_n]
  if (n_y_dot > 0)
    blas.SCAL(n_y_rows*n_y_dot, beta, y_dot, incy_dot);

  // Compute [yd_1 .. yd_n] = alpha*A*[xd_1 .. xd_n]
  if (n_x_dot > 0) {
    if (incx_dot == 1)
      blas.GEMM(trans, Teuchos::NO_TRANS, n_A_rows, n_dot, n_A_cols,
                alpha, A, lda, x_dot, n_x_rows, 1.0, y_dot, n_y_rows);
    else
      for (OrdinalType i=0; i<n_x_dot; i++)
        blas.GEMV(trans, m, n, alpha, A, lda, x_dot+i*incx_dot*n_x_rows,
                  incx_dot, 1.0, y_dot+i*incy_dot*n_y_rows, incy_dot);
  }

  // Compute [yd_1 .. yd_n] += diag([alphad_1 .. alphad_n])*A*x
  if (n_alpha_dot > 0) {
    if (gemv_Ax.size() != std::size_t(n))
      gemv_Ax.resize(n);
    blas.GEMV(trans, m, n, 1.0, A, lda, x, incx, 0.0, &gemv_Ax[0],
              OrdinalType(1));
    for (OrdinalType i=0; i<n_alpha_dot; i++)
      blas.AXPY(n_y_rows, alpha_dot[i], &gemv_Ax[0], OrdinalType(1),
                y_dot+i*incy_dot*n_y_rows, incy_dot);
  }

  // Compute [yd_1 .. yd_n] += alpha*[Ad_1*x .. Ad_n*x]
  for (OrdinalType i=0; i<n_A_dot; i++)
    blas.GEMV(trans, m, n, alpha, A_dot+i*lda_dot*n, lda_dot, x, incx, 1.0,
              y_dot+i*incy_dot*n_y_rows, incy_dot);

  // Compute [yd_1 .. yd_n] += diag([betad_1 .. betad_n])*y
  for (OrdinalType i=0; i<n_beta_dot; i++)
    blas.AXPY(n_y_rows, beta_dot[i], y, incy, y_dot+i*incy_dot*n_y_rows,
              incy_dot);

  // Compute y = alpha*A*x + beta*y
  if (n_alpha_dot > 0) {
    blas.SCAL(n_y_rows, beta, y, incy);
    blas.AXPY(n_y_rows, alpha, &gemv_Ax[0], OrdinalType(1), y, incy);
  }
  else
    blas.GEMV(trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename x_type, typename y_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_GER(const OrdinalType m,
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
        const OrdinalType n_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_x_dot != n_dot && n_x_dot != 0) ||
                     (n_y_dot != n_dot && n_y_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_GER(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif

  // Compute [Ad_1 .. Ad_n] += [alphad_1*x*y^T .. alphad_n*x*y^T]
  for (OrdinalType i=0; i<n_alpha_dot; i++)
    blas.GER(m, n, alpha_dot[i], x, incx, y, incy, A_dot+i*lda_dot*n, lda_dot);

  // Compute [Ad_1 .. Ad_n] += alpha*[xd_1*y^T .. xd_n*y^T]
  for (OrdinalType i=0; i<n_x_dot; i++)
    blas.GER(m, n, alpha, x_dot+i*incx_dot*m, incx_dot, y, incy,
             A_dot+i*lda_dot*n, lda_dot);

  // Compute [Ad_1 .. Ad_n] += alpha*x*[yd_1 .. yd_n]
  if (n_y_dot > 0)
    blas.GER(m, n*n_y_dot, alpha, x, incx, y_dot, incy_dot, A_dot, lda_dot);

  // Compute A = alpha*x*y^T + A
  blas.GER(m, n, alpha, x, incx, y, incy, A, lda);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type, typename B_type,
          typename beta_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_GEMM(Teuchos::ETransp transa,
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
         const OrdinalType n_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_B_dot != n_dot && n_B_dot != 0) ||
                     (n_beta_dot != n_dot && n_beta_dot != 0) ||
                     (n_C_dot != n_dot && n_C_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_GEMM(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif
  OrdinalType n_A_cols = k;
  if (transa != Teuchos::NO_TRANS) {
    n_A_cols = m;
  }

  OrdinalType n_B_cols = n;
  if (transb != Teuchos::NO_TRANS) {
    n_B_cols = k;
  }

  // Compute [Cd_1 .. Cd_n] = beta*[Cd_1 .. Cd_n]
  if (n_C_dot > 0) {
    if (ldc_dot == m)
      blas.SCAL(m*n*n_C_dot, beta, C_dot, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_C_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.SCAL(m, beta, C_dot+i*ldc_dot*n+j*ldc_dot, OrdinalType(1));
  }

  // Compute [Cd_1 .. Cd_n] += alpha*A*[Bd_1 .. Bd_n]
  for (OrdinalType i=0; i<n_B_dot; i++)
    blas.GEMM(transa, transb, m, n, k, alpha, A, lda, B_dot+i*ldb_dot*n_B_cols,
              ldb_dot, 1.0, C_dot+i*ldc_dot*n, ldc_dot);

  // Compute [Cd_1 .. Cd_n] += [alphad_1*A*B .. alphad_n*A*B]
  if (n_alpha_dot > 0) {
    if (gemm_AB.size() != std::size_t(m*n))
      gemm_AB.resize(m*n);
    blas.GEMM(transa, transb, m, n, k, 1.0, A, lda, B, ldb, 0.0, &gemm_AB[0],
              OrdinalType(m));
    if (ldc_dot == m)
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        blas.AXPY(m*n, alpha_dot[i], &gemm_AB[0], OrdinalType(1),
                  C_dot+i*ldc_dot*n, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.AXPY(m, alpha_dot[i], &gemm_AB[j*m], OrdinalType(1),
                    C_dot+i*ldc_dot*n+j*ldc_dot, OrdinalType(1));
  }

  // Compute [Cd_1 .. Cd_n] += alpha*[Ad_1*B .. Ad_n*B]
  for (OrdinalType i=0; i<n_A_dot; i++)
    blas.GEMM(transa, transb, m, n, k, alpha, A_dot+i*lda_dot*n_A_cols,
              lda_dot, B, ldb, 1.0, C_dot+i*ldc_dot*n, ldc_dot);

  // Compute [Cd_1 .. Cd_n] += [betad_1*C .. betad_n*C]
  if (ldc == m && ldc_dot == m)
    for (OrdinalType i=0; i<n_beta_dot; i++)
      blas.AXPY(m*n, beta_dot[i], C, OrdinalType(1), C_dot+i*ldc_dot*n,
                OrdinalType(1));
  else
    for (OrdinalType i=0; i<n_beta_dot; i++)
      for (OrdinalType j=0; j<n; j++)
        blas.AXPY(m, beta_dot[i], C+j*ldc, OrdinalType(1),
                  C_dot+i*ldc_dot*n+j*ldc_dot, OrdinalType(1));

  // Compute C = alpha*A*B + beta*C
  if (n_alpha_dot > 0) {
    if (ldc == m) {
      blas.SCAL(m*n, beta, C, OrdinalType(1));
      blas.AXPY(m*n, alpha, &gemm_AB[0], OrdinalType(1), C, OrdinalType(1));
    }
    else
      for (OrdinalType j=0; j<n; j++) {
        blas.SCAL(m, beta, C+j*ldc, OrdinalType(1));
        blas.AXPY(m, alpha, &gemm_AB[j*m], OrdinalType(1), C+j*ldc,
                  OrdinalType(1));
      }
  }
  else
    blas.GEMM(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type, typename B_type,
          typename beta_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_SYMM(Teuchos::ESide side, Teuchos::EUplo uplo,
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
         const OrdinalType n_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_B_dot != n_dot && n_B_dot != 0) ||
                     (n_beta_dot != n_dot && n_beta_dot != 0) ||
                     (n_C_dot != n_dot && n_C_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_SYMM(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif
  OrdinalType n_A_cols = m;
  if (side == Teuchos::RIGHT_SIDE) {
    n_A_cols = n;
  }

  // Compute [Cd_1 .. Cd_n] = beta*[Cd_1 .. Cd_n]
  if (n_C_dot > 0) {
    if (ldc_dot == m)
      blas.SCAL(m*n*n_C_dot, beta, C_dot, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_C_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.SCAL(m, beta, C_dot+i*ldc_dot*n+j*ldc_dot, OrdinalType(1));
  }

  // Compute [Cd_1 .. Cd_n] += alpha*A*[Bd_1 .. Bd_n]
  for (OrdinalType i=0; i<n_B_dot; i++)
    blas.SYMM(side, uplo, m, n, alpha, A, lda, B_dot+i*ldb_dot*n,
              ldb_dot, 1.0, C_dot+i*ldc_dot*n, ldc_dot);

  // Compute [Cd_1 .. Cd_n] += [alphad_1*A*B .. alphad_n*A*B]
  if (n_alpha_dot > 0) {
    if (gemm_AB.size() != std::size_t(m*n))
      gemm_AB.resize(m*n);
    blas.SYMM(side, uplo, m, n, 1.0, A, lda, B, ldb, 0.0, &gemm_AB[0],
              OrdinalType(m));
    if (ldc_dot == m)
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        blas.AXPY(m*n, alpha_dot[i], &gemm_AB[0], OrdinalType(1),
                  C_dot+i*ldc_dot*n, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.AXPY(m, alpha_dot[i], &gemm_AB[j*m], OrdinalType(1),
                    C_dot+i*ldc_dot*n+j*ldc_dot, OrdinalType(1));
  }

  // Compute [Cd_1 .. Cd_n] += alpha*[Ad_1*B .. Ad_n*B]
  for (OrdinalType i=0; i<n_A_dot; i++)
    blas.SYMM(side, uplo, m, n, alpha, A_dot+i*lda_dot*n_A_cols, lda_dot, B,
              ldb, 1.0, C_dot+i*ldc_dot*n, ldc_dot);

  // Compute [Cd_1 .. Cd_n] += [betad_1*C .. betad_n*C]
  if (ldc == m && ldc_dot == m)
    for (OrdinalType i=0; i<n_beta_dot; i++)
      blas.AXPY(m*n, beta_dot[i], C, OrdinalType(1), C_dot+i*ldc_dot*n,
                OrdinalType(1));
  else
    for (OrdinalType i=0; i<n_beta_dot; i++)
      for (OrdinalType j=0; j<n; j++)
        blas.AXPY(m, beta_dot[i], C+j*ldc, OrdinalType(1),
                  C_dot+i*ldc_dot*n+j*ldc_dot, OrdinalType(1));

  // Compute C = alpha*A*B + beta*C
  if (n_alpha_dot > 0) {
    if (ldc == m) {
      blas.SCAL(m*n, beta, C, OrdinalType(1));
      blas.AXPY(m*n, alpha, &gemm_AB[0], OrdinalType(1), C, OrdinalType(1));
    }
    else
      for (OrdinalType j=0; j<n; j++) {
        blas.SCAL(m, beta, C+j*ldc, OrdinalType(1));
        blas.AXPY(m, alpha, &gemm_AB[j*m], OrdinalType(1), C+j*ldc,
                  OrdinalType(1));
      }
  }
  else
    blas.SYMM(side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_TRMM(Teuchos::ESide side,
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
         const OrdinalType n_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_B_dot != n_dot && n_B_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_TRMM(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif
  OrdinalType n_A_cols = m;
  if (side == Teuchos::RIGHT_SIDE) {
    n_A_cols = n;
  }

  // Compute [Bd_1 .. Bd_n] = alpha*A*[Bd_1 .. Bd_n]
  for (OrdinalType i=0; i<n_B_dot; i++)
    blas.TRMM(side, uplo, transa, diag, m, n, alpha, A, lda, B_dot+i*ldb_dot*n,
              ldb_dot);

  // Compute [Bd_1 .. Bd_n] += [alphad_1*A*B .. alphad_n*A*B]
  if (n_alpha_dot > 0) {
    if (gemm_AB.size() != std::size_t(m*n))
      gemm_AB.resize(m*n);
    if (ldb == m)
      blas.COPY(m*n, B, OrdinalType(1), &gemm_AB[0], OrdinalType(1));
    else
      for (OrdinalType j=0; j<n; j++)
        blas.COPY(m, B+j*ldb, OrdinalType(1), &gemm_AB[j*m], OrdinalType(1));
    blas.TRMM(side, uplo, transa, diag, m, n, 1.0, A, lda, &gemm_AB[0],
              OrdinalType(m));
    if (ldb_dot == m)
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        blas.AXPY(m*n, alpha_dot[i], &gemm_AB[0], OrdinalType(1),
                  B_dot+i*ldb_dot*n, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.AXPY(m, alpha_dot[i], &gemm_AB[j*m], OrdinalType(1),
                    B_dot+i*ldb_dot*n+j*ldb_dot, OrdinalType(1));
  }

  // Compute [Bd_1 .. Bd_n] += alpha*[Ad_1*B .. Ad_n*B]
  if (n_A_dot > 0) {
    if (gemm_AB.size() != std::size_t(m*n))
      gemm_AB.resize(m*n);
    for (OrdinalType i=0; i<n_A_dot; i++) {
      if (ldb == m)
        blas.COPY(m*n, B, OrdinalType(1), &gemm_AB[0], OrdinalType(1));
      else
        for (OrdinalType j=0; j<n; j++)
          blas.COPY(m, B+j*ldb, OrdinalType(1), &gemm_AB[j*m], OrdinalType(1));
      blas.TRMM(side, uplo, transa, Teuchos::NON_UNIT_DIAG, m, n, alpha,
                A_dot+i*lda_dot*n_A_cols, lda_dot, &gemm_AB[0],
                OrdinalType(m));
      if (ldb_dot == m)
        blas.AXPY(m*n, 1.0, &gemm_AB[0], OrdinalType(1),
                  B_dot+i*ldb_dot*n, OrdinalType(1));
      else
        for (OrdinalType j=0; j<n; j++)
          blas.AXPY(m, 1.0, &gemm_AB[j*m], OrdinalType(1),
                    B_dot+i*ldb_dot*n+j*ldb_dot, OrdinalType(1));
    }
  }

  // Compute B = alpha*A*B
  if (n_alpha_dot > 0 && n_A_dot == 0) {
    if (ldb == m) {
      blas.SCAL(m*n, 0.0, B, OrdinalType(1));
      blas.AXPY(m*n, alpha, &gemm_AB[0], OrdinalType(1), B, OrdinalType(1));
    }
    else
      for (OrdinalType j=0; j<n; j++) {
        blas.SCAL(m, 0.0, B+j*ldb, OrdinalType(1));
        blas.AXPY(m, alpha, &gemm_AB[j*m], OrdinalType(1), B+j*ldb,
                  OrdinalType(1));
      }
  }
  else
    blas.TRMM(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);
}

template <typename OrdinalType, typename FadType>
template <typename alpha_type, typename A_type>
void
Sacado::Fad::BLAS<OrdinalType,FadType>::
Fad_TRSM(Teuchos::ESide side,
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
         const OrdinalType n_dot) const
{
#ifdef SACADO_DEBUG
  // Check sizes are consistent
  TEUCHOS_TEST_FOR_EXCEPTION((n_alpha_dot != n_dot && n_alpha_dot != 0) ||
                     (n_A_dot != n_dot && n_A_dot != 0) ||
                     (n_B_dot != n_dot && n_B_dot != 0),
                     std::logic_error,
                     "BLAS::Fad_TRSM(): All arguments must have " <<
                     "the same number of derivative components, or none");
#endif
  OrdinalType n_A_cols = m;
  if (side == Teuchos::RIGHT_SIDE) {
    n_A_cols = n;
  }

  // Compute [Bd_1 .. Bd_n] = alpha*[Bd_1 .. Bd_n]
  if (n_B_dot > 0) {
    if (ldb_dot == m)
      blas.SCAL(m*n*n_B_dot, alpha, B_dot, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_B_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.SCAL(m, alpha, B_dot+i*ldb_dot*n+j*ldb_dot, OrdinalType(1));
  }

  // Compute [Bd_1 .. Bd_n] += [alphad_1*B .. alphad_n*B]
  if (n_alpha_dot > 0) {
    if (ldb == m && ldb_dot == m)
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        blas.AXPY(m*n, alpha_dot[i], B, OrdinalType(1),
                  B_dot+i*ldb_dot*n, OrdinalType(1));
    else
      for (OrdinalType i=0; i<n_alpha_dot; i++)
        for (OrdinalType j=0; j<n; j++)
          blas.AXPY(m, alpha_dot[i], B+j*ldb, OrdinalType(1),
                    B_dot+i*ldb_dot*n+j*ldb_dot, OrdinalType(1));
  }

  // Solve A*X = alpha*B
  blas.TRSM(side, uplo, transa, diag, m, n, alpha, A, lda, B, ldb);

  // Compute [Bd_1 .. Bd_n] -= [Ad_1*X .. Ad_n*X]
  if (n_A_dot > 0) {
    if (gemm_AB.size() != std::size_t(m*n))
      gemm_AB.resize(m*n);
    for (OrdinalType i=0; i<n_A_dot; i++) {
      if (ldb == m)
        blas.COPY(m*n, B, OrdinalType(1), &gemm_AB[0], OrdinalType(1));
      else
        for (OrdinalType j=0; j<n; j++)
          blas.COPY(m, B+j*ldb, OrdinalType(1), &gemm_AB[j*m], OrdinalType(1));
      blas.TRMM(side, uplo, transa, Teuchos::NON_UNIT_DIAG, m, n, 1.0,
                A_dot+i*lda_dot*n_A_cols, lda_dot, &gemm_AB[0],
                OrdinalType(m));
      if (ldb_dot == m)
        blas.AXPY(m*n, -1.0, &gemm_AB[0], OrdinalType(1),
                  B_dot+i*ldb_dot*n, OrdinalType(1));
      else
        for (OrdinalType j=0; j<n; j++)
          blas.AXPY(m, -1.0, &gemm_AB[j*m], OrdinalType(1),
                    B_dot+i*ldb_dot*n+j*ldb_dot, OrdinalType(1));
    }
  }

  // Solve A*[Xd_1 .. Xd_n] = [Bd_1 .. Bd_n]
  if (side == Teuchos::LEFT_SIDE)
    blas.TRSM(side, uplo, transa, diag, m, n*n_dot, 1.0, A, lda, B_dot,
              ldb_dot);
  else
    for (OrdinalType i=0; i<n_dot; i++)
      blas.TRSM(side, uplo, transa, diag, m, n, 1.0, A, lda, B_dot+i*ldb_dot*n,
                ldb_dot);
}
