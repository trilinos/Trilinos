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

//BEGINVersion
#include "gtest/gtest.h"
#include <stk_util/stk_config.h>
#ifdef STK_HAVE_KOKKOS
//currently stk *ALWAYS* has Kokkos, but in the future perhaps some
//portions of stk won't depend on Kokkos...
#include <Kokkos_Core.hpp>
#endif
#if defined(STK_HAVE_STKIO) && defined(STK_HAS_SEACAS_IOSS)
#include <Ioss_Version.h>
#endif
#if defined(HAVE_STK_Trilinos) || defined(STK_BUILT_FOR_SIERRA)
#include <Trilinos_version.h>
#endif
#include <stk_util/Version.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <iostream>
#include <string>

TEST(stkHowTo, reportVersion)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) != 1) { GTEST_SKIP(); }

  const std::string stk_version = stk::version_string();
  std::cout << "This program is using STK Version: " << stk_version << std::endl;
  std::cout << "Value of STK_VERSION macro: " << STK_VERSION << std::endl;
  std::cout << "Kokkos Version: " << KOKKOS_VERSION << std::endl;
#if defined(STK_HAVE_STKIO) && defined(STK_HAS_SEACAS_IOSS)
  std::cout << "SEACAS/IOSS Version: " << Ioss::Version() << std::endl;
#endif
#ifdef TRILINOS_VERSION_STRING
  std::cout << "Trilinos Version: " << TRILINOS_VERSION_STRING << std::endl;
#endif
}
//ENDVersion

