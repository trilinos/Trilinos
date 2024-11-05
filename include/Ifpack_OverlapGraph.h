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

#ifndef IFPACK_OVERLAPGRAPH_H
#define IFPACK_OVERLAPGRAPH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "Ifpack_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Import.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_RowMatrix;

namespace Teuchos {
  class ParameterList;
}

//! Ifpack_OverlapGraph: Constructs a graph for use with Ifpack preconditioners.

class Ifpack_OverlapGraph: public Epetra_Object {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Epetra_CrsGraph.
  /*! Creates an Ifpack_OverlapGraph object from the user graph.
    \param In
           UserMatrixGraph_in - Graph from user matrix.
  */
  Ifpack_OverlapGraph(const Teuchos::RefCountPtr<const Epetra_CrsGraph>& UserMatrixGraph_in, int OverlapLevel_in);

  //! Constructor using Epetra_RowMatrix.
  /*! Creates an Ifpack_OverlapGraph object from the user graph implicitly defined by the
         Epetra_RowMatrix interface.
    \param In
            RowMatrix - An object that has implemented the Epetra_RowMatrix interface.
  */
  Ifpack_OverlapGraph(const Teuchos::RefCountPtr<const Epetra_RowMatrix>& UserMatrix_in, int OverlapLevel_in);

  //! Copy constructor.
  Ifpack_OverlapGraph(const Ifpack_OverlapGraph & Source);

  //! Ifpack_CrsIlut Destructor
  virtual ~Ifpack_OverlapGraph() {};
  //@}

  //@{ \name Attribute access methods.

  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the configure argument
     '--enable-ifpack-teuchos' was used.
     This method recognizes the name: level_overlap, which is case insensitive.
     The ParameterEntry must have type int.
  */
  int SetParameters(const Teuchos::ParameterList& parameterlist,
                    bool cerr_warning_if_unused=false);

  //! Returns the overlap graph object.
  const Epetra_CrsGraph & OverlapGraph() const {return(*OverlapGraph_);}

  //! Returns the RowMap associated with the overlap graph.
  const Epetra_BlockMap & OverlapRowMap() const {return(*OverlapRowMap_);}

  //! Returns the overlap graph object.
  const Epetra_Import & OverlapImporter() const {return(*OverlapImporter_);}

  //! Returns the level of overlap used to create this graph.
  /*! The graph created by this class uses a recursive definition 0f overlap.
      Level one overlap is created by copying all off-processor rows that are
      reached to be at least one column of the rows that are on processor.
      Level two overlap is the same process used on the level one graph.
  */
  int OverlapLevel() const {return(OverlapLevel_);}
  //@}

  //@{ \name Epetra_Object print method (allows use of << operator with this class).

  void Print(std::ostream& os) const {
    using std::endl;

    os << endl;
    if (UserMatrix_!=Teuchos::null)
      os << "Overlap Graph created using the user's Epetra_RowMatrix object" << endl;
    else
      os << "Overlap Graph created using the user's Epetra_CrsGraph object" << endl;

    os << " Level of Overlap = " << OverlapLevel_ << endl;
    OverlapGraph_->Print(os);
    return;
  }
  //@}

 protected:

  int ConstructOverlapGraph(const Teuchos::RefCountPtr<const Epetra_CrsGraph>& UserMatrixGraph);
  Teuchos::RefCountPtr<Epetra_CrsGraph> OverlapGraph_;
  Teuchos::RefCountPtr<const Epetra_CrsGraph> UserMatrixGraph_;
  Teuchos::RefCountPtr<const Epetra_RowMatrix> UserMatrix_;
  Teuchos::RefCountPtr<Epetra_BlockMap> OverlapRowMap_;
  Teuchos::RefCountPtr<Epetra_Import> OverlapImporter_;
  int OverlapLevel_;
  bool IsOverlapped_;
};
#endif // IFPACK_OVERLAPGRAPH_H
