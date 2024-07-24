// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_ACCESS_HPP
#define TPETRA_ACCESS_HPP

namespace Tpetra
{
namespace Access
{
  // Structs for Access tags, these should not be used by user code
  struct ReadOnlyStruct {};
  struct OverwriteAllStruct {};
  struct ReadWriteStruct {};


  //Tag indicating intent to read up-to-date data, but not modify.
  constexpr struct ReadOnlyStruct  ReadOnly  = ReadOnlyStruct();
  //Tag indicating intent to completely overwrite existing data.
  constexpr struct OverwriteAllStruct OverwriteAll = OverwriteAllStruct();
  //Tag indicating intent to both read up-to-date data and modify it.
  constexpr struct ReadWriteStruct ReadWrite = ReadWriteStruct();
}
}

#endif

