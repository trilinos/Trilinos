/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef IFPACK2_OVERLAPGRAPH_HPP
#define IFPACK2_OVERLAPGRAPH_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_RCP.hpp"
#include "Ifpack2_CreateOverlapGraph.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace Ifpack2 {

//! Ifpack2::OverlapGraph constructs an overlapped graph.

template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class OverlapGraph : public Teuchos::Describable {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Tpetra::CrsGraph.
  /*! Creates an Ifpack2::OverlapGraph object from the user graph. 
    \param In
           UserMatrixGraph_in - Graph from user matrix.
  */
  OverlapGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& UserMatrixGraph_in, int OverlapLevel_in);

  //! Constructor using Tpetra_RowMatrix.
  /*! Creates an Ifpack2_OverlapGraph object from the user graph implicitly defined by the
	 Tpetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Tpetra_RowMatrix interface.
  */
//  Ifpack2_OverlapGraph(const Teuchos::RCP<const Tpetra_RowMatrix>& UserMatrix_in, int OverlapLevel_in);
  
  //! Copy constructor.
  OverlapGraph(const OverlapGraph<LocalOrdinal,GlobalOrdinal,Node> & Source);

  //! Destructor
  virtual ~OverlapGraph() {};
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
  const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> & getOverlapGraph() const {return(*OverlapGraph_);}
    
  //! Returns the RowMap associated with the overlap graph.
  const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> & getOverlapRowMap() const {return(*OverlapRowMap_);}
    
  //! Returns the overlap graph object.
  const Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> & getOverlapImporter() const {return(*OverlapImporter_);}
    
  //! Returns the level of overlap used to create this graph.
  /*! The graph created by this class uses a recursive definition 0f overlap.
      Level one overlap is created by copying all off-processor rows that are
      reached to be at least one column of the rows that are on processor.
      Level two overlap is the same process used on the level one graph.
  */
  int OverlapLevel() const {return(OverlapLevel_);}
  //@}

  //@{ \name Tpetra_Object print method (allows use of << operator with this class).

//  void Print(ostream& os) const {
//    os << endl;
//    if (UserMatrix_!=Teuchos::null) 
//      os << "Overlap Graph created using the user's Tpetra_RowMatrix object" << endl;
//    else
//      os << "Overlap Graph created using the user's Tpetra_CrsGraph object" << endl;
//    
//    os << " Level of Overlap = " << OverlapLevel_ << endl;
//    OverlapGraph_->Print(os);
//    return;
//  }
  //@}

 protected:

  int ConstructOverlapGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& UserMatrixGraph);
  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > OverlapGraph_;
  Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > UserMatrixGraph_;
//  Teuchos::RCP<const Tpetra_RowMatrix> UserMatrix_;
  Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > OverlapRowMap_;
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > OverlapImporter_;
  int OverlapLevel_;
  bool IsOverlapped_;
};

template<class LocalOrdinal, class GlobalOrdinal, class Node>
OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>::OverlapGraph(const Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& UserMatrixGraph_in, int OverlapLevel_in)
 : UserMatrixGraph_(UserMatrixGraph_in),
   OverlapLevel_(OverlapLevel_in),
   IsOverlapped_(OverlapLevel_in>0 && UserMatrixGraph_in->getDomainMap()->isDistributed())
{
  OverlapGraph_ = CreateOverlapGraph(UserMatrixGraph_, OverlapLevel_);
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>::OverlapGraph(const OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>& Source)
 : UserMatrixGraph_(Source.UserMatrixGraph_),
   OverlapRowMap_(Source.OverlapRowMap_),
   OverlapLevel_(Source.OverlapLevel_),
   IsOverlapped_(Source.IsOverlapped_)
{
  if (IsOverlapped_) {
    if (OverlapGraph_!=Teuchos::null) OverlapGraph_ = Teuchos::rcp(new Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>(*OverlapGraph_));
    if (OverlapRowMap_!=Teuchos::null) OverlapRowMap_ = Teuchos::rcp(new Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>(*OverlapRowMap_));
  }
}

}//namespace Ifpack2

#endif // IFPACK2_OVERLAPGRAPH_HPP
