#ifndef GENSQP_YUEPETRAVECTOR_H
#define GENSQP_YUEPETRAVECTOR_H

#include "GenSQP_Vector.hpp"
#include "Epetra_MultiVector.h"


/** \class GenSQP::YUEpetraVector
    \brief The GenSQP::Vector / (y,u) Epetra_MultiVector adapter class.
      
    Holds a pointer to two Epetra_MultiVectors y_epetra_vec and u_epetra_vec and implements the member functions 
    of the GenSQP::Vector class.
    Common use: optimal control. The vector y_epetra_vec represents the state variables, the vector
    u_epetra_vec represents the control variables.
*/


namespace GenSQP {

class YUEpetraVector : public Vector {

private:
  
  Teuchos::RefCountPtr<Epetra_MultiVector>  y_epetra_vec_;
  Teuchos::RefCountPtr<Epetra_MultiVector>  u_epetra_vec_;

public:

  YUEpetraVector( const Teuchos::RefCountPtr<Epetra_MultiVector> &y_epetra_vec,
                  const Teuchos::RefCountPtr<Epetra_MultiVector> &u_epetra_vec );

  /** \name Overridden from Vector */
  //@{

  double innerProd( const Vector &x ) const;

  void linComb( const double &alpha, const Vector &x, const double &beta );

  void Scale( const double &alpha );

  void Set( const double &alpha );

  void Set( const double &alpha, const Vector &x );

  Teuchos::RefCountPtr<Vector> createVector() const;

  //@}
  
  /** Returns a reference counted pointer to the private y_epetra_vec data container ("state variables").
  */
  Teuchos::RefCountPtr<const Epetra_MultiVector> getYVector() const;

  /** Returns a reference counted pointer to the private u_epetra_vec data container ("control variables").
  */
  Teuchos::RefCountPtr<const Epetra_MultiVector> getUVector() const;


}; // class YUEpetraVector

} // namespace GenSQP

#endif
