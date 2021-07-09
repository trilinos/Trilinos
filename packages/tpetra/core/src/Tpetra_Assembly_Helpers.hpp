// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_ASSEMBLY_HELPERS_HPP
#define TPETRA_ASSEMBLY_HELPERS_HPP

namespace Tpetra {

namespace Impl {
// Helper function to to apply an operation to each member of  a
// c++11 parameter pack since parameter expansion only happens
// within functions, constructors, and initializer_lists
template <typename... Args>
inline void foreach_pack(Args &&... args) {}
} // namespace Impl


template <typename... Args>
void beginAssembly(Args &&... args)
{
  // use the comma operator to transform a potentially void function call
  // into a argument to allow proper parameter expansion for c++11
  Impl::foreach_pack( (args.beginAssembly(),1)... );

  // using c++17 the code would be
  // (args.beginAssembly()...);
}

template <typename... Args>
void endAssembly(Args &&... args)
{
  // use the comma operator to transform a potentially void function call
  // into a argument to allow proper parameter expansion for c++11
  Impl::foreach_pack( (args.endAssembly(),1)... );

  // using c++17 the code would be
  // (args.endAssembly()...);

}

template <typename... Args>
void beginModify(Args &&... args)
{
  // use the comma operator to transform a potentially void function call
  // into a argument to allow proper parameter expansion for c++11
  Impl::foreach_pack( (args.beginModify(),1)... );

  // using c++17 the code would be
  // (args.beginModify()...);
}

template <typename... Args>
void endModify(Args &&... args)
{
  // use the comma operator to transform a potentially void function call
  // into a argument to allow proper parameter expansion for c++11
  Impl::foreach_pack( (args.endModify(),1)... );

  // using c++17 the code would be
  // (args.endModify()...);

}

}// namespace Tpetra

#endif // TPETRA_ASSEMBLY_HELPERS_HPP
