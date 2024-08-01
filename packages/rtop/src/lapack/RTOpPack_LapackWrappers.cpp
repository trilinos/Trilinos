// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_LapackWrappers.hpp"

const Teuchos::Tuple<char,RTOpPack::NUM_ETRANS_ARGS>
RTOpPack::transpMap = Teuchos::tuple('N', 'T', 'C');
