#ifndef IFPACK_REORDERING_H
#define IFPACK_REORDERING_H

#include "Ifpack_ConfigDefs.h"

namespace Teuchos {
  class ParameterList;
}
class Epetra_MultiVector;
class Ifpack_Graph;
class Epetra_RowMatrix;

//! Ifpack_Reordering: basic class for reordering for a Ifpack_Graph.
/*!
Class Ifpack_Reordering is a pure virtual class that defines the 
structure of all Ifpack reordering.

An Ifpack_Reordering object is a tool, used in class Ifpack_Preconditioner,
to reorder the localized matrix (with or without overlap). As its
basic usage is for localized matrices, this class takes care of
reordering the \e local rows only.

If IFPACK is configure with Teuchos support, method SetParameters()
should be adopted. Otherwise, users can set parameters (one at-a-time),
using methods SetParameter(), for integers and doubles.

Ifpack_Preconditioner objects overload the << operator. Derived
classes should specify a Print() method, that will be used in
operator <<.

\author Marzio Sala, SNL 9214

\date Oct-04
*/

class Ifpack_Reordering {

public:

  //! Destructor.
  virtual ~Ifpack_Reordering() {};
  
  //! Sets integer parameters `Name'.
  virtual int SetParameter(const string Name, const int Value) = 0;
 
  //! Sets double parameters `Name'.
  virtual int SetParameter(const string Name, const double Value) = 0;

#ifdef HAVE_IFPACK_TEUCHOS  
  //! Sets all parameters.
  virtual int SetParameters(Teuchos::ParameterList& List) = 0;
#endif

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Ifpack_Graph& Graph) = 0;

  //! Computes all it is necessary to initialize the reordering object.
  virtual int Compute(const Epetra_RowMatrix& Matrix) = 0;

  //! Returns \c true is the reordering object has been successfully initialized, false otherwise.
  virtual bool IsComputed() const = 0;

  //! Returns the reordered index of row \c i.
  virtual int Reorder(const int i) const = 0;

  //! Returns the inverse reordered index of row \c i.
  virtual int InvReorder(const int i) const = 0;

  //! Applies reordering to multivector X, whose local length equals the number of local rows.
  virtual int P(const Epetra_MultiVector& Xorig,
		Epetra_MultiVector& X) const = 0;

  //! Applies inverse reordering to multivector X, whose local length equals the number of local rows.
  virtual int Pinv(const Epetra_MultiVector& Xorig,
		   Epetra_MultiVector& X) const = 0;

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual ostream& Print(std::ostream& os) const = 0;

}; 

inline ostream& operator<<(ostream& os, const Ifpack_Reordering& obj)
{
  return(obj.Print(os));
}

#endif
