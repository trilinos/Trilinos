/*@HEADER
// ***********************************************************************
// 
//       Tifpack: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef TIFPACK_OVERLAPGRAPH_HPP
#define TIFPACK_OVERLAPGRAPH_HPP

#include "Tifpack_ConfigDefs.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_RCP.hpp"
#include "Tifpack_CreateOverlapGraph.hpp"

class Tpetra_Comm;
class Tpetra_BlockMap;
class Tpetra_RowMatrix;

namespace Teuchos {
  class ParameterList;
}

namespace Tifpack {

//! Tifpack_OverlapGraph: Constructs a graph for use with Tifpack preconditioners.

template<class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class OverlapGraph : public Teuchos::Describable {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Tpetra_CrsGraph.
  /*! Creates an Tifpack_OverlapGraph object from the user graph. 
    \param In
           UserMatrixGraph_in - Graph from user matrix.
  */
  OverlapGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& UserMatrixGraph_in, int OverlapLevel_in);

  //! Constructor using Tpetra_RowMatrix.
  /*! Creates an Tifpack_OverlapGraph object from the user graph implicitly defined by the
	 Tpetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Tpetra_RowMatrix interface.
  */
//  Tifpack_OverlapGraph(const Teuchos::RCP<const Tpetra_RowMatrix>& UserMatrix_in, int OverlapLevel_in);
  
  //! Copy constructor.
  OverlapGraph(const OverlapGraph<LocalOrdinal,GlobalOrdinal,Node> & Source);

  //! Tifpack_CrsIlut Destructor
  virtual ~OverlapGraph() {};
  //@}

  //@{ \name Atribute access methods.
    
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
  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > OverlapGraph_;
  Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > UserMatrixGraph_;
//  Teuchos::RCP<const Tpetra_RowMatrix> UserMatrix_;
  Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > OverlapRowMap_;
  Teuchos::RCP<Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > OverlapImporter_;
  int OverlapLevel_;
  bool IsOverlapped_;
};

template<class LocalOrdinal, class GlobalOrdinal, class Node>
OverlapGraph<LocalOrdinal,GlobalOrdinal,Node>::OverlapGraph(const Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> >& UserMatrixGraph_in, int OverlapLevel_in)
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

}//namespace Tifpack

#endif // TIFPACK_OVERLAPGRAPH_HPP
