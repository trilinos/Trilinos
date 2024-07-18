// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

template <typename Scalar>
int mainCuda(bool test_flat, bool test_orig, bool test_lin, bool test_block,
             bool symmetric, int device_id)
{
  return 0 ;
}

template int mainCuda<float>(bool, bool, bool, bool, bool, int);
template int mainCuda<double>(bool, bool, bool, bool, bool, int);
