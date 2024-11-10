/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_METISREORDERING_H
#define IFPACK_METISREORDERING_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Reordering.h"
#include "Epetra_RowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Comm;
class Ifpack_Graph;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Import;

//! Ifpack_METISReordering: A class to reorder a graph using METIS.

class Ifpack_METISReordering : public Ifpack_Reordering {

public:

  //! Constructor.
  Ifpack_METISReordering();

  //! Destructor.
  virtual ~Ifpack_METISReordering() {};

  //! Sets integer parameters `Name'.
  virtual int SetParameter(const std::string Name, const int Value)
  {
    if (Name == "partitioner: use symmetric graph")
      UseSymmetricGraph_ = (bool)Value;
    return(0);
  }

  //! Sets double parameters `Name'.
  virtual int SetParameter(const std::string /* Name */, const double /* Value */)
  {
    return(0);
  };

  //! Sets all the parameters for the partitioner (none at moment).
  virtual int SetParameters(Teuchos::ParameterList& List)
  {
    UseSymmetricGraph_ = List.get("partitioner: use symmetric graph",
                                  UseSymmetricGraph_);

    return(0);
  }

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Ifpack_Graph& Graph);

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Epetra_RowMatrix& Matrix);

  //! Returns \c true is the reordering object has been successfully initialized, false otherwise.
  virtual bool IsComputed() const
  {
    return(IsComputed_);
  }

  //! Returns the reordered index of row \c i.
  virtual int Reorder(const int i) const;

  //! Returns the inverse reordered index of row \c i.
  virtual int InvReorder(const int i) const;

  //! Applies reordering to multivector Xorig, whose local length equals the number of local rows, stores result in X.
  virtual int P(const Epetra_MultiVector& Xorig,
                Epetra_MultiVector& X) const;

  //! Applies inverse reordering to multivector Xorig, whose local length equals the number of local rows, stores result in X.
  virtual int Pinv(const Epetra_MultiVector& Xorig,
                   Epetra_MultiVector& X) const;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& Print(std::ostream& os) const;

private:
  //! If \c true, the graph has to be symmetrized before calling METIS.
  bool UseSymmetricGraph_;
  //! Number of local rows in the graph.
  int NumMyRows_;
  //! If \c true, the reordering has been successfully computed.
  bool IsComputed_;
  //! Contains the reordering.
  std::vector<int> Reorder_;
  //! Contains the inverse reordering.
  std::vector<int> InvReorder_;

}; // class Ifpack_METISReordering

#endif // IFPACK_METISREORDERING_H
