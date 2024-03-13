// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_REPORTHANDLER_H_
#define KRINO_INCLUDE_AKRI_REPORTHANDLER_H_

#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/parallel/ParallelReduceBool.hpp>

#define ParallelThrowRequire(parallel,expr)            STK_ThrowRequire(stk::is_true_on_all_procs(parallel,expr))
#define ParallelThrowRequireMsg(parallel,expr,message) STK_ThrowRequireMsg(stk::is_true_on_all_procs(parallel,expr),message)
#define ParallelThrowAssert(parallel,expr)             STK_ThrowAssert(stk::is_true_on_all_procs(parallel,expr))
#define ParallelThrowAssertMsg(parallel,expr,message)  STK_ThrowAssertMsg(stk::is_true_on_all_procs(parallel,expr),message)


#endif /* KRINO_INCLUDE_AKRI_REPORTHANDLER_H_ */
