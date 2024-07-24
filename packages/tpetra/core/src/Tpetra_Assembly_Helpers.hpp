// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
