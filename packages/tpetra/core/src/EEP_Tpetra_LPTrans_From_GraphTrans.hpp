//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER
                                                                                                    
#ifndef TPETRA_LINEARPROBLEM_GRAPHTRANS_H
#define TPETRA_LINEARPROBLEM_GRAPHTRANS_H

#include <EEP_Tpetra_Transform.hpp>

#include <Tpetra_Export_decl.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_LinearProblem_decl.hpp>
#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_Map_decl.hpp>

namespace Tpetra {

//! Tpetra::LinearProblem_GraphTrans: Adaptation of a Tpetra::CrsGraph Transform to a Tpetra::LinearProblem Transform
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class LinearProblem_GraphTrans : public SameTypeTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
{
  StructuralSameTypeTransform< Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > & graphTrans_;

  Tpetra::Import<LocalOrdinal,
                 GlobalOrdinal,
                 Node> * Importer_;
  Tpetra::Export<LocalOrdinal,
                 GlobalOrdinal,
                 Node> * MatExporter_;
  Tpetra::Export<LocalOrdinal,
                 GlobalOrdinal,
                 Node> * VecExporter_;

  Tpetra::LinearProblem<Scalar,
                        LocalOrdinal,
                        GlobalOrdinal,
                        Node> * OldProblem_;
  const Tpetra::CrsGraph<LocalOrdinal,
                         GlobalOrdinal,
                         Node> * OldGraph_;
  Tpetra::CrsMatrix<Scalar,
                    LocalOrdinal,
                    GlobalOrdinal,
                    Node> * OldMatrix_;
  Tpetra::MultiVector<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> * OldLHS_;
  Tpetra::MultiVector<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> * OldRHS_;
  const Tpetra::Map<LocalOrdinal,
                    GlobalOrdinal,
                    Node> * OldRowMap_;

  Tpetra::LinearProblem<Scalar,
                        LocalOrdinal,
                        GlobalOrdinal,
                        Node> * NewProblem_;
  Tpetra::CrsMatrix<Scalar,
                    LocalOrdinal,
                    GlobalOrdinal,
                    Node> * NewMatrix_;
  Tpetra::MultiVector<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> * NewLHS_;
  Tpetra::MultiVector<Scalar,
                      LocalOrdinal,
                      GlobalOrdinal,
                      Node> * NewRHS_;

 public:

  //! EpetraExt::LinearProblem_GraphTrans Destructor
  ~LinearProblem_GraphTrans();

  //! EpetraExt::LinearProblem_GraphTrans Constructor
  /*! Constructs a LinearProblem Transform based on the input CrsGraph Transform
      \param In
      graph_trans - Base Tpetra::CrsGraph Transform from which a consistent Tpetra::LinearProblem  Transform is generated
 */
  LinearProblem_GraphTrans( StructuralSameTypeTransform< Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > & graph_trans )
  : graphTrans_(graph_trans),
    Importer_(0),
    MatExporter_(0),
    VecExporter_(0),
    OldProblem_(0),
    OldGraph_(0),
    OldMatrix_(0),
    OldLHS_(0),
    OldRHS_(0),
    OldRowMap_(0),
    NewProblem_(0),
    NewMatrix_(0),
    NewLHS_(0),
    NewRHS_(0)
  {
    std::cout << "EEP Passing through LinearProblem_GraphTrans::constructor()..." << std::endl;
  }

  //! Constructs an Tpetra::LinearProblem from the original using the same row transformation given by the Tpetra::CrsGraph Transform
  /*! 
      \param In
      orig - Original Tpetra::LinearProblem to be transformed.
      \return Tpetra::LinearProblem generated by transformation operation
  */
  typename SameTypeTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewTypeRef operator()( typename SameTypeTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalTypeRef orig );

  //! Forward migration of data from original to transformed object
  bool fwd();

  //! Reverse migration of data from transformed to original object
  bool rvs();

};

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
LinearProblem_GraphTrans<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
~LinearProblem_GraphTrans()
{
  if( MatExporter_ ) delete MatExporter_;
  if( VecExporter_ ) delete VecExporter_;
  if( Importer_ ) delete Importer_;

  if( NewProblem_ ) delete NewProblem_;
  if( NewRHS_ ) delete NewRHS_;
  if( NewLHS_ ) delete NewLHS_;
  if( NewMatrix_ ) delete NewMatrix_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SameTypeTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::NewTypeRef
LinearProblem_GraphTrans<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()( typename SameTypeTransform< Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::OriginalTypeRef orig )
{
  std::cout << "EEP Entering tpetra/core/src/EEP_Tpetra_LPTrans_From_GraphTrans.hpp LinearProblem_GraphTrans::operator()..." << std::endl;
  
  OldProblem_ = &orig;
  OldMatrix_ = dynamic_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>*>( orig.getMatrix().get() );
  OldGraph_ = dynamic_cast<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>*>( OldMatrix_->getGraph().get() );
  OldRHS_ = orig.getRHS().get();
  OldLHS_ = orig.getLHS().get();
  OldRowMap_ = dynamic_cast<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>*>( OldMatrix_->getRowMap().get() );

  std::cout << "EEP In tpetra/core/src/EEP_Tpetra_LPTrans_From_GraphTrans.hpp LinearProblem_GraphTrans::operator(), pos 001" << std::endl;
  Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> & NewGraph = graphTrans_( *(const_cast<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>*>(OldGraph_)) );
  std::cout << "EEP In tpetra/core/src/EEP_Tpetra_LPTrans_From_GraphTrans.hpp LinearProblem_GraphTrans::operator(), pos 002"
            << ": NewGraph = " << NewGraph
	    << std::endl;
  NewMatrix_ = new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>( Teuchos::rcp<Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>(&NewGraph) );
#if 0 // AquiToDo // EEP____
  Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> & NewRowMap = const_cast<Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>&>(NewGraph.getRowMap());

  NewRHS_ = new Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>( NewRowMap, 1 );
  NewLHS_ = new Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>( NewRowMap, 1 );

  MatExporter_ = new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>( *OldRowMap_, NewRowMap );
  VecExporter_ = new Tpetra::Export<LocalOrdinal, GlobalOrdinal, Node>( *OldRowMap_, NewRowMap );
  Importer_ = new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>( *OldRowMap_, NewRowMap );

  std::cout << "EEP In tpetra/core/src/EEP_Tpetra_LPTrans_From_GraphTrans.hpp LinearProblem_GraphTrans::operator(), pos 003" << std::endl;
  NewProblem_ = new Tpetra::LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>( NewMatrix_, NewLHS_, NewRHS_ );
#endif
  std::cout << "EEP In tpetra/core/src/EEP_Tpetra_LPTrans_From_GraphTrans.hpp LinearProblem_GraphTrans::operator(), pos 004" << std::endl;

  std::cout << "EEP Leaving tpetra/core/src/EEP_Tpetra_LPTrans_From_GraphTrans.hpp LinearProblem_GraphTrans::operator()" << std::endl;
  return *NewProblem_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool
LinearProblem_GraphTrans<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
fwd()
{
  NewLHS_->doExport( *OldLHS_, *VecExporter_, INSERT );
  NewRHS_->doExport( *OldRHS_, *VecExporter_, INSERT );
  NewMatrix_->doExport( *OldMatrix_, *MatExporter_, INSERT );

  return true;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
bool
LinearProblem_GraphTrans<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
rvs()
{
  OldLHS_->doImport( *NewLHS_, *Importer_, INSERT );
//  OldRHS_->Import( *NewRHS_, *Importer_, Insert ); // As in the original epetraext/src/transform/ file
//  OldMatrix_->Import( *NewMatrix_, *Importer_, Insert ); // As in the original epetraext/src/transform/ file

  return true;
}

} // namespace Tpetra

#endif // TPETRA_LINEARPROBLEM_GRAPHTRANS_H
