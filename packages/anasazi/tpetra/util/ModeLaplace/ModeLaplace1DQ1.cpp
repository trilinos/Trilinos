// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//

#include "ModeLaplace1DQ1.hpp"
#include <Tpetra_Map_fwd.hpp>

using LO = typename Tpetra::Map<>::local_ordinal_type;
using GO = typename Tpetra::Map<>::global_ordinal_type;
using Node = typename Tpetra::Map<>::node_type;

template<>
const int ModeLaplace1DQ1<double,LO,GO,Node>::dofEle = 2;
template<>
const int ModeLaplace1DQ1<double,LO,GO,Node>::maxConnect = 3;
#ifndef M_PI
template<>
const double ModeLaplace1DQ1<double,LO,GO,Node>::M_PI = 3.14159265358979323846;
#endif

