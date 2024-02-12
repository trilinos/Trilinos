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

#include "gtest/gtest.h"
#include "stk_util/parallel/MPI.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace {

#if defined ( STK_HAS_MPI )
TEST(MPIDatatype, UnsignedBasicTypes)
{
  EXPECT_EQ(MPI_BYTE              , sierra::MPI::Datatype<unsigned char>::type());
  EXPECT_EQ(MPI_UNSIGNED          , sierra::MPI::Datatype<unsigned int>::type());
  EXPECT_EQ(MPI_UNSIGNED_SHORT    , sierra::MPI::Datatype<unsigned short>::type());
  EXPECT_EQ(MPI_UNSIGNED_LONG_LONG, sierra::MPI::Datatype<unsigned long long>::type());

  MPI_Datatype goldValue_unsigned_long{MPI_DATATYPE_NULL};
  goldValue_unsigned_long = std::is_same<uint64_t, unsigned long>::value ? MPI_UINT64_T : MPI_UNSIGNED_LONG;
  EXPECT_EQ(goldValue_unsigned_long, sierra::MPI::Datatype<unsigned long>::type());

  EXPECT_EQ(MPI_UINT64_T, sierra::MPI::Datatype<uint64_t>::type());

  MPI_Datatype goldValue_size_t{MPI_DATATYPE_NULL};
  if(std::is_same<uint64_t, size_t>::value){
    goldValue_size_t = MPI_UINT64_T;
  } else if(std::is_same<unsigned long, size_t>::value){
    goldValue_size_t = MPI_UNSIGNED_LONG;
  } else if(std::is_same<unsigned, size_t>::value){
    goldValue_size_t = MPI_UNSIGNED;
  }
  EXPECT_EQ(goldValue_size_t, sierra::MPI::Datatype<size_t>::type());
}
#endif

}
