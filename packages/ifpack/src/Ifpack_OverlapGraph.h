/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * October 20, 2002, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef IFPACK_OVERLAPGRAPH_H
#define IFPACK_OVERLAPGRAPH_H

#include "Epetra_Object.h"
#include "Epetra_CrsGraph.h"
class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_RowMatrix;
class Epetra_Import;

//! Ifpack_OverlapGraph: Constructs a graph for use with Ifpack preconditioners.

class Ifpack_OverlapGraph: public Epetra_Object {

 public:
  //@{ \name Constructors/Destructor
  //! Constructor using Epetra_CrsGraph.
  /*! Creates an Ifpack_OverlapGraph object from the user graph. 
    \param In
           UserMatrixGraph - Graph from user matrix.
  */
  Ifpack_OverlapGraph(const Epetra_CrsGraph * UserMatrixGraph, int OverlapLevel);

  //! Constructor using Epetra_RowMatrix.
  /*! Creates an Ifpack_OverlapGraph object from the user graph implicitly defined by the
	 Epetra_RowMatrix interface. 
    \param In
            RowMatrix - An object that has implemented the Epetra_RowMatrix interface.
  */
  Ifpack_OverlapGraph(const Epetra_RowMatrix * UserMatrix, int OverlapLevel);
  
  //! Copy constructor.
  Ifpack_OverlapGraph(const Ifpack_OverlapGraph & Source);

  //! Ifpack_CrsIlut Destructor
  virtual ~Ifpack_OverlapGraph();
  //@}

  //@{ \name Atribute access methods.
    
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

  void Print(ostream& os) const {

  os << endl;
  if (UserMatrix_!=0) 
    os << "Overlap Graph created using the user's Epetra_RowMatrix object" << endl;
  else
        os << "Overlap Graph created using the user's Epetra_CrsGraph object" << endl;

  os << " Level of Overlap = " << OverlapLevel_ << endl;
  OverlapGraph_->Print(os);
  return;
}
  //@}

 protected:

  int ConstructOverlapGraph(const Epetra_CrsGraph * UserMatrixGraph);
  Epetra_CrsGraph * OverlapGraph_;
  const Epetra_CrsGraph * UserMatrixGraph_;
  const Epetra_RowMatrix * UserMatrix_;
  Epetra_BlockMap * OverlapRowMap_;
  Epetra_Import * OverlapImporter_;
  int OverlapLevel_;
  bool IsOverlapped_;
};
#endif // IFPACK_OVERLAPGRAPH_H
