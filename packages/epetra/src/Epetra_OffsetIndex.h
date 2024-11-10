/*
//@HEADER
// ************************************************************************
//
//               Epetra: Linear Algebra Services Package
//                 Copyright 2011 Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_OFFSETINDEX_H
#define EPETRA_OFFSETINDEX_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_Object.h"

class Epetra_Import;
class Epetra_Export;
class Epetra_CrsGraph;
class Epetra_Distributor;

//! Epetra_OffsetIndex: This class builds index for efficient mapping of data from one Epetra_CrsGraph based object to another.

/*! Epetra_OffsetIndex generates and index of offsets allowing direct access to data
 for Import/Export operations on Epetra_CrsGraph based objects such as Epetra_CrsMatrix.
*/

class EPETRA_LIB_DLL_EXPORT Epetra_OffsetIndex: public Epetra_Object {

 public:

  //! Constructs a Epetra_OffsetIndex object from the graphs and an importer.
  Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph,
                      const Epetra_CrsGraph & TargetGraph,
                      Epetra_Import & Importer );

  //! Constructs a Epetra_OffsetIndex object from the graphs and an exporter.
  Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph,
                      const Epetra_CrsGraph & TargetGraph,
                      Epetra_Export & Exporter );

  //! Epetra_OffsetIndex copy constructor.
  Epetra_OffsetIndex(const Epetra_OffsetIndex & Indexor);

  //! Epetra_OffsetIndex destructor.
  virtual ~Epetra_OffsetIndex(void);

  //! @name Print object to an output stream
  //@{
  virtual void Print(std::ostream & os) const;
  //@}

  //! Accessor
  int ** SameOffsets() const { return SameOffsets_; }

  //! Accessor
  int ** PermuteOffsets() const { return PermuteOffsets_; }

  //! Accessor
  int ** RemoteOffsets() const { return RemoteOffsets_; }

 private:

  template<typename int_type>
  void GenerateLocalOffsets_( const Epetra_CrsGraph & SourceGraph,
                              const Epetra_CrsGraph & TargetGraph,
                              const int * PermuteLIDs );

  template<typename int_type>
  void GenerateRemoteOffsets_( const Epetra_CrsGraph & SourceGraph,
                               const Epetra_CrsGraph & TargetGraph,
                               const int * ExportLIDs,
                               const int * RemoteLIDs,
                               Epetra_Distributor & Distor );

  //! Epetra_OffsetIndex copy constructor.
  Epetra_OffsetIndex & operator=(const Epetra_OffsetIndex & Indexor);
 public:

  int NumSame_;
  int ** SameOffsets_;
  int NumPermute_;
  int ** PermuteOffsets_;
  int NumExport_;
  int NumRemote_;
  int ** RemoteOffsets_;

  bool DataOwned_;
};

#endif /* EPETRA_OFFSETINDEX_H */
