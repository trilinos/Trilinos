// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_ACCESS_HPP
#define XPETRA_ACCESS_HPP

namespace Xpetra {
namespace Access {
// Structs for Access tags, these should not be used by user code
struct ReadOnlyStruct {};
struct OverwriteAllStruct {};
struct ReadWriteStruct {};

// Tag indicating intent to read up-to-date data, but not modify.
constexpr struct ReadOnlyStruct ReadOnly = ReadOnlyStruct();
// Tag indicating intent to completely overwrite existing data.
constexpr struct OverwriteAllStruct OverwriteAll = OverwriteAllStruct();
// Tag indicating intent to both read up-to-date data and modify it.
constexpr struct ReadWriteStruct ReadWrite = ReadWriteStruct();
}  // namespace Access
}  // namespace Xpetra
#endif  // XPETRA_ACCESS_HPP
