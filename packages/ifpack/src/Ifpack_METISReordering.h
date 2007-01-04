#ifndef IFPACK_METISREORDERING_H
#define IFPACK_METISREORDERING_H

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Reordering.h"
#include "Teuchos_ParameterList.hpp"
class Epetra_Comm;
class Epetra_RowMatrix;
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
  ~Ifpack_METISReordering() {};

  //! Sets integer parameters `Name'.
  virtual int SetParameter(const string Name, const int Value)
  {
    if (Name == "partitioner: use symmetric graph")
      UseSymmetricGraph_ = (bool)Value;
    return(0);
  }
 
  //! Sets double parameters `Name'.
  virtual int SetParameter(const string Name, const double Value)
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
  virtual ostream& Print(std::ostream& os) const;

private:
  //! If \c true, the graph has to be symmetrized before calling METIS.
  bool UseSymmetricGraph_;
  //! Number of local rows in the graph.
  int NumMyRows_;
  //! If \c true, the reordering has been successfully computed.
  bool IsComputed_;
  //! Contains the reordering.
  vector<int> Reorder_;
  //! Contains the inverse reordering.
  vector<int> InvReorder_;
  Ifpack_Graph* Graph_;

}; // class Ifpack_METISReordering

#endif // IFPACK_METISREORDERING_H
