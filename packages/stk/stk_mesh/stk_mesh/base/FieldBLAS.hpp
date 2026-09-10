// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_MESH_BASE_FIELDBLAS_HPP
#define STK_MESH_BASE_FIELDBLAS_HPP

#include <stk_util/stk_config.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/baseImpl/FieldBLASImpl.hpp>
#include <stk_mesh/baseImpl/NgpFieldBLASImpl.hpp>

#include "stk_util/ngp/NgpSpaces.hpp"

namespace stk::mesh {

//==============================================================================
// axpy: y[i] = a*x[i] + y[i]
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_axpy(T alpha,
    const FieldBase& xField,
    const FieldBase& yField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  field_blas::impl::field_axpy_impl(alpha, xField, yField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField, const Selector& selector, const typename NgpSpace::exec_space execSpace) {
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  ngp_field_blas::impl::field_axpy_impl<NgpSpace>(alpha, xField, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField, const Selector& selector) {
  auto execSpace = typename NgpSpace::exec_space();
  field_axpy<NgpSpace>(alpha, xField, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField, const typename NgpSpace::exec_space execSpace) {
  const Selector selector = selectField(xField) & selectField(yField);
  field_axpy<NgpSpace>(alpha, xField, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField) {
  auto execSpace = typename NgpSpace::exec_space();
  field_axpy<NgpSpace>(alpha, xField, yField, execSpace);
}

template <typename T>
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  field_axpy<ngp::HostSpace>(alpha, xField, yField, selector);
}

template <typename T>
void field_axpy(T alpha, const FieldBase& xField, const FieldBase& yField)
{
  field_axpy<ngp::HostSpace>(alpha, xField, yField);
}
//==============================================================================
// axpby: y[i] = a*x[i] + b*y[i]
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_axpby(T alpha,
    const FieldBase& xField,
    T beta,
    const FieldBase& yField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  field_blas::impl::field_axpby_impl(alpha, xField, beta, yField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const Selector& selector, const typename NgpSpace::exec_space execSpace) {
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  ngp_field_blas::impl::field_axpby_impl<NgpSpace>(alpha, xField, beta, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const Selector& selector) {
  auto execSpace = typename NgpSpace::exec_space();
  field_axpby<NgpSpace>(alpha, xField, beta, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const typename NgpSpace::exec_space execSpace) {
  const Selector selector = selectField(xField) & selectField(yField);
  field_axpby<NgpSpace>(alpha, xField, beta, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField) {
  auto execSpace = typename NgpSpace::exec_space();
  field_axpby<NgpSpace>(alpha, xField, beta, yField, execSpace);
}

template <typename T>
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const Selector& selector)
{
  field_axpby<ngp::HostSpace>(alpha, xField, beta, yField, selector);
}

template <typename T>
void field_axpby(T alpha, const FieldBase& xField, T beta, const FieldBase& yField)
{
  field_axpby<ngp::HostSpace>(alpha, xField, beta, yField);
}

//==============================================================================
// axpbyz: z[i] = a*x[i] + b*y[i]
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_axpbyz(T alpha,
    const FieldBase& xField,
    T beta,
    const FieldBase& yField,
    const FieldBase& zField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, zField));

  field_blas::impl::field_axpbyz_impl(alpha, xField, beta, yField, zField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_axpbyz(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const FieldBase& zField, const Selector& selector, const typename NgpSpace::exec_space execSpace) {
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, zField));

  ngp_field_blas::impl::field_axpbyz_impl<NgpSpace>(alpha, xField, beta, yField, zField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpbyz(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const FieldBase& zField, const Selector& selector) {
  auto execSpace = typename NgpSpace::exec_space();
  field_axpbyz<NgpSpace>(alpha, xField, beta, yField, zField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpbyz(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const FieldBase& zField, const typename NgpSpace::exec_space execSpace) {
  const Selector selector = selectField(xField) & selectField(yField) & selectField(zField);
  field_axpbyz<NgpSpace>(alpha, xField, beta, yField, zField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_axpbyz(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const FieldBase& zField) {
  auto execSpace = typename NgpSpace::exec_space();
  field_axpbyz<NgpSpace>(alpha, xField, beta, yField, zField, execSpace);
}

template <typename T>
void field_axpbyz(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
  field_axpbyz<ngp::HostSpace>(alpha, xField, beta, yField, zField, selector);
}

template <typename T>
void field_axpbyz(T alpha, const FieldBase& xField, T beta, const FieldBase& yField, const FieldBase& zField)
{
  field_axpbyz<ngp::HostSpace>(alpha, xField, beta, yField, zField);
}

//==============================================================================
// product: z[i] = x[i] * y[i]
//
template <typename NgpSpace>
  requires ngp::is_host_space<NgpSpace>
void field_product(const FieldBase& xField,
    const FieldBase& yField,
    const FieldBase& zField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));
  STK_ThrowRequire(field_blas::impl::is_compatible(yField, zField));

  field_blas::impl::field_product_impl(xField, yField, zField, selector);
}

template <typename NgpSpace> requires ngp::is_device_space<NgpSpace>
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));
  STK_ThrowRequire(field_blas::impl::is_compatible(yField, zField));

  ngp_field_blas::impl::field_product_impl<NgpSpace>(xField, yField, zField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_product<NgpSpace>(xField, yField, zField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField) & selectField(yField) & selectField(zField);
  field_product<NgpSpace>(xField, yField, zField, selector, execSpace);
}


template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_product<NgpSpace>(xField, yField, zField, execSpace);
}

template <typename T>
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
  field_product<ngp::HostSpace>(xField, yField, zField, selector);
}

template <typename T>
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField)
{
  field_product<ngp::HostSpace>(xField, yField, zField);
}

inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField, const Selector& selector)
{
  field_product<ngp::HostSpace>(xField, yField, zField, selector);
}

inline
void field_product(const FieldBase& xField, const FieldBase& yField, const FieldBase& zField)
{
  field_product<ngp::HostSpace>(xField, yField, zField);
}

//==============================================================================
// copy: y[i] = x[i]
//
template <typename NgpSpace>
  requires ngp::is_host_space<NgpSpace>
void field_copy(const FieldBase& xField,
    const FieldBase& yField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  field_blas::impl::field_copy_impl(xField, yField, selector);
}

template <typename NgpSpace> requires ngp::is_device_space<NgpSpace>
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  ngp_field_blas::impl::field_copy_impl<NgpSpace>(xField, yField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_copy<NgpSpace>(xField, yField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_copy(const FieldBase& xField, const FieldBase& yField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_copy<NgpSpace>(xField, yField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_copy(const FieldBase& xField, const FieldBase& yField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_copy<NgpSpace>(xField, yField, execSpace);
}

template <typename T>
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  field_copy<ngp::HostSpace>(xField, yField, selector);
}

inline
void field_copy(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  field_copy<ngp::HostSpace>(xField, yField, selector);
}

inline
void field_copy(const FieldBase& xField, const FieldBase& yField)
{
  field_copy<ngp::HostSpace>(xField, yField);
}

//==============================================================================
// dot: global_sum( sum_i( x[i]*y[i] ) )
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_dot(T& result,
    const FieldBase& xField,
    const FieldBase& yField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  field_blas::impl::field_dot_impl(result, xField, yField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField, const Selector& selector, const typename NgpSpace::exec_space execSpace) {
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  ngp_field_blas::impl::field_dot_impl<NgpSpace>(result, xField, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField, const Selector& selector) {
  auto execSpace = typename NgpSpace::exec_space();
  field_dot<NgpSpace>(result, xField, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField, const typename NgpSpace::exec_space execSpace) {
  const Selector selector = selectField(xField) & selectField(yField);
  field_dot<NgpSpace>(result, xField, yField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField) {
  auto execSpace = typename NgpSpace::exec_space();
  field_dot<NgpSpace>(result, xField, yField, execSpace);
}

template <typename T, Layout xLayout, Layout yLayout>
T field_dot(const Field<T, xLayout>& xField, const Field<T, yLayout>& yField, const Selector& selector)
{
  T result{};
  field_dot<ngp::HostSpace>(result, xField, yField, selector);
  return result;
}

template <typename T, Layout xLayout, Layout yLayout>
T field_dot(const Field<T, xLayout>& xField, const Field<T, yLayout>& yField)
{
  T result{};
  field_dot<ngp::HostSpace>(result, xField, yField);
  return result;
}

template <typename T>
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  field_dot<ngp::HostSpace>(result, xField, yField, selector);
}

template <typename T>
void field_dot(T& result, const FieldBase& xField, const FieldBase& yField)
{
  field_dot<ngp::HostSpace>(result, xField, yField);
}

//==============================================================================
// nrm2: sqrt( global_sum( sum_i( x[i]*x[i] )))
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_nrm2(
    T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_nrm2_impl(result, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_nrm2(T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_nrm2_impl<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_nrm2(T& result, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_nrm2<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_nrm2(T& result, const FieldBase& xField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_nrm2<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_nrm2(T& result, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_nrm2<NgpSpace>(result, xField, execSpace);
}

template <typename T, Layout xLayout>
T field_nrm2(const Field<T, xLayout>& xField, const Selector& selector)
{
  T result;
  field_nrm2<ngp::HostSpace>(result, xField, selector);
  return result;
}

template <typename T, Layout xLayout>
T field_nrm2(const Field<T, xLayout>& xField)
{
  T result;
  field_nrm2<ngp::HostSpace>(result, xField);
  return result;
}

template <typename T>
void field_nrm2(T& result, const FieldBase& xField, const Selector& selector)
{
  field_nrm2<ngp::HostSpace>(result, xField, selector);
}

template <typename T>
void field_nrm2(T& result, const FieldBase& xField)
{
  field_nrm2<ngp::HostSpace>(result, xField);
}

//==============================================================================
// scale: x[i] = a*x[i]
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_scale(
    T alpha, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_scale_impl(alpha, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_scale(T alpha, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_scale_impl<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_scale(T alpha, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_scale<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_scale(T alpha, const FieldBase& xField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_scale<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_scale(T alpha, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_scale<NgpSpace>(alpha, xField, execSpace);
}

template <typename T>
void field_scale(T alpha, const FieldBase& xField, const Selector& selector)
{
  field_scale<ngp::HostSpace>(alpha, xField, selector);
}

template <typename T>
void field_scale(T alpha, const FieldBase& xField)
{
  field_scale<ngp::HostSpace>(alpha, xField);
}

//==============================================================================
// fill: x[i] = a
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill(
    T alpha, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_fill_impl(alpha, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_fill_impl<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_fill<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xField, execSpace);
}

template <typename T>
void field_fill(T alpha, const FieldBase& xField, const Selector& selector)
{
  field_fill<ngp::HostSpace>(alpha, xField, selector);
}

template <typename T>
void field_fill(T alpha, const FieldBase& xField)
{
  field_fill<ngp::HostSpace>(alpha, xField);
}

//==============================================================================
// fill: x[i, c] = a
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill(T alpha,
    const FieldBase& xField,
    int component,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_fill_impl(alpha, xField, component, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, int component, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_fill_impl<NgpSpace>(alpha, xField, component, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, int component, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xField, component, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, int component, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_fill<NgpSpace>(alpha, xField, component, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const FieldBase& xField, int component)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xField, component, execSpace);
}

template <typename T>
void field_fill(T alpha, const FieldBase& xField, int component, const Selector& selector)
{
  field_fill<ngp::HostSpace>(alpha, xField, component, selector);
}

template <typename T>
void field_fill(T alpha, const FieldBase& xField, int component)
{
  field_fill<ngp::HostSpace>(alpha, xField, component);
}

//==============================================================================
// fill: x_j[i] = a
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill(T alpha,
    const std::vector<const FieldBase*>& xFields,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
    field_blas::impl::field_fill_impl(alpha, *xField, selector);
  }
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
  }
  ngp_field_blas::impl::field_fill_impl<NgpSpace>(alpha, xFields, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xFields, selector, execSpace);
}

template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill(
    T alpha, const std::vector<const FieldBase*>& xFields, const typename NgpSpace::exec_space /*execSpace*/)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
    const Selector selector = selectField(*xField);
    field_blas::impl::field_fill_impl(alpha, *xField, selector);
  }
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, const typename NgpSpace::exec_space execSpace)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
  }
  Selector selector = Selector(*xFields.front());
  std::for_each(std::cbegin(xFields)+1, std::cend(xFields), [&selector](const auto* field) {
    selector &= Selector(*field);
  });
  ngp_field_blas::impl::field_fill_impl<NgpSpace>(alpha, xFields, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xFields, execSpace);
}

template <typename T>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, const Selector& selector)
{
  field_fill<ngp::HostSpace>(alpha, xFields, selector);
}

template <typename T>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields)
{
  field_fill<ngp::HostSpace>(alpha, xFields);
}

//==============================================================================
// fill: x_j[i, c] = a
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill(T alpha,
    const std::vector<const FieldBase*>& xFields,
    int component,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
    field_blas::impl::field_fill_impl(alpha, *xField, component, selector);
  }
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, int component, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
  }
  ngp_field_blas::impl::field_fill_impl<NgpSpace>(alpha, xFields, component, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, int component, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xFields, component, selector, execSpace);
}

template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill(T alpha,
    const std::vector<const FieldBase*>& xFields,
    int component,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
    const Selector selector = selectField(*xField);
    field_blas::impl::field_fill_impl(alpha, *xField, component, selector);
  }
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, int component, const typename NgpSpace::exec_space execSpace)
{
  for (const FieldBase* xField : xFields) {
    STK_ThrowRequire(field_blas::impl::is_compatible<T>(*xField));
  }
  Selector selector = Selector(*xFields.front());
  std::for_each(std::cbegin(xFields)+1, std::cend(xFields), [&selector](const auto* field) {
    selector &= Selector(*field);
  });
  ngp_field_blas::impl::field_fill_impl<NgpSpace>(alpha, xFields, component, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, int component)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill<NgpSpace>(alpha, xFields, component, execSpace);
}

template <typename T>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, int component, const Selector& selector)
{
  field_fill<ngp::HostSpace>(alpha, xFields, component, selector);
}

template <typename T>
void field_fill(T alpha, const std::vector<const FieldBase*>& xFields, int component)
{
  field_fill<ngp::HostSpace>(alpha, xFields, component);
}

//==============================================================================
// fill_component: x[i, comp] = a[comp]
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_fill_component(const T* alpha,
    const FieldBase& xField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_fill_component_impl(alpha, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill_component(const T* alpha, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill_component<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill_component(const T* alpha, const FieldBase& xField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_fill_component<NgpSpace>(alpha, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_fill_component(const T* alpha, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_fill_component<NgpSpace>(alpha, xField, execSpace);
}

template <typename T>
void field_fill_component(const T* alpha, const FieldBase& xField, const Selector& selector)
{
  field_fill_component<ngp::HostSpace>(alpha, xField, selector);
}

template <typename T>
void field_fill_component(const T* alpha, const FieldBase& xField)
{
  field_fill_component<ngp::HostSpace>(alpha, xField);
}

//==============================================================================
// swap: x[i] = y[i]
//       y[i] = x[i]
//
template <typename NgpSpace>
  requires ngp::is_host_space<NgpSpace>
void field_swap(const FieldBase& xField,
    const FieldBase& yField,
    const Selector& selector,
    const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  field_blas::impl::field_swap_impl(xField, yField, selector);
}

template <typename NgpSpace> requires ngp::is_device_space<NgpSpace>
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible(xField, yField));

  ngp_field_blas::impl::field_swap_impl<NgpSpace>(xField, yField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_swap<NgpSpace>(xField, yField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_swap(const FieldBase& xField, const FieldBase& yField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField) & selectField(yField);
  field_swap<NgpSpace>(xField, yField, selector, execSpace);
}

template <typename NgpSpace> requires ngp::is_ngp_space<NgpSpace>
void field_swap(const FieldBase& xField, const FieldBase& yField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_swap<NgpSpace>(xField, yField, execSpace);
}

template <typename T>
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  field_swap<ngp::HostSpace>(xField, yField, selector);
}

inline
void field_swap(const FieldBase& xField, const FieldBase& yField, const Selector& selector)
{
  field_swap<ngp::HostSpace>(xField, yField, selector);
}

inline
void field_swap(const FieldBase& xField, const FieldBase& yField)
{
  field_swap<ngp::HostSpace>(xField, yField);
}

//==============================================================================
// asum: global_sum( sum_i( abs(x[i]) ))
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_asum(
    T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_asum_impl(result, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_asum(T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_asum_impl<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_asum(T& result, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_asum<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_asum(T& result, const FieldBase& xField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_asum<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_asum(T& result, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_asum<NgpSpace>(result, xField, execSpace);
}

template <typename T, Layout xLayout>
T field_asum(const Field<T, xLayout>& xField, const Selector& selector)
{
  T result;
  field_asum<ngp::HostSpace>(result, xField, selector);
  return result;
}

template <typename T, Layout xLayout>
inline
T field_asum(const Field<T, xLayout>& xField)
{
  T result;
  field_asum<ngp::HostSpace>(result, xField);
  return result;
}

template <typename T>
void field_asum(T& result, const FieldBase& xField, const Selector& selector)
{
  field_asum<ngp::HostSpace>(result, xField, selector);
}

template <typename T>
void field_asum(T& result, const FieldBase& xFieldBase)
{
  field_asum<ngp::HostSpace>(result, xFieldBase);
}

//==============================================================================
// amax: global_max( max_i( abs(x[i]) ))
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_amax(
    T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_amax_impl(result, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_amax(T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_amax_impl<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_amax(T& result, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_amax<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_amax(T& result, const FieldBase& xField, const typename NgpSpace::exec_space execSpace)
{
  const Selector selector = selectField(xField);
  field_amax<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_amax(T& result, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_amax<NgpSpace>(result, xField, execSpace);
}

template <typename T, Layout xLayout>
T field_amax(const Field<T, xLayout>& xField, const Selector& selector)
{
  T result;
  field_amax<ngp::HostSpace>(result, xField, selector);
  return result;
}

template <typename T, Layout xLayout>
T field_amax(const Field<T, xLayout>& xField)
{
  T result;
  field_amax<ngp::HostSpace>(result, xField);
  return result;
}

template <typename T>
void field_amax(T& result, const FieldBase& xField, const Selector& selector)
{
  field_amax<ngp::HostSpace>(result, xField, selector);
}

template <typename T>
void field_amax(T& result, const FieldBase& xField)
{
  field_amax<ngp::HostSpace>(result, xField);
}

//==============================================================================
// amin: global_min( min_i( abs(x[i]) ))
//
template <typename NgpSpace, typename T>
  requires ngp::is_host_space<NgpSpace>
void field_amin(
    T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space /*execSpace*/)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  field_blas::impl::field_amin_impl(result, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_device_space<NgpSpace>
void field_amin(T& result, const FieldBase& xField, const Selector& selector, const typename NgpSpace::exec_space execSpace)
{
  STK_ThrowRequire(field_blas::impl::is_compatible<T>(xField));

  ngp_field_blas::impl::field_amin_impl<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_amin(T& result, const FieldBase& xField, const Selector& selector)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_amin<NgpSpace>(result, xField, selector, execSpace);
}

template <typename NgpSpace, typename T>
  requires ngp::is_ngp_space<NgpSpace>
void field_amin(T& result, const FieldBase& xField, const typename NgpSpace::exec_space /*execSpace*/)
{
  const Selector selector = selectField(xField);
  field_amin<NgpSpace>(result, xField, selector);
}

template <typename NgpSpace, typename T> requires ngp::is_ngp_space<NgpSpace>
void field_amin(T& result, const FieldBase& xField)
{
  auto execSpace = typename NgpSpace::exec_space();
  field_amin<NgpSpace>(result, xField, execSpace);
}

template <typename T, Layout xLayout>
T field_amin(const Field<T, xLayout>& xField, const Selector& selector)
{
  T result;
  field_amin<ngp::HostSpace>(result, xField, selector);
  return result;
}

template <typename T, Layout xLayout>
inline
T field_amin(const Field<T, xLayout>& xField)
{
  T result;
  field_amin<ngp::HostSpace>(result, xField);
  return result;
}

template <typename T>
void field_amin(T& result, const FieldBase& xField, const Selector& selector)
{
  field_amin<ngp::HostSpace>(result, xField, selector);
}

template <typename T>
void field_amin(T& result, const FieldBase& xField)
{
  field_amin<ngp::HostSpace>(result, xField);
}

} // stk::mesh

#endif // STK_MESH_BASE_FIELDBLAS_HPP

