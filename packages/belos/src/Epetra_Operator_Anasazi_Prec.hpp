// File Epetra_Operator_Anasazi_Prec.hpp: interface for the AnasaziPetra class.
//
#ifndef EPETRA_OPERATOR_ANASAZI_PREC_H
#define EPETRA_OPERATOR_ANASAZI_PREC_H

#include "Epetra_Operator.h"

#include <cassert>

#include "AnasaziMatrix.hpp"
#include "AnasaziPrecondition.hpp"

///////////////////////////////////////////////////////////////
//--------template class Epetra_Operator_Anasazi_Prec--------------------

template <class TYPE> 
class Epetra_Operator_Anasazi_Prec : public AnasaziPrecondition<TYPE> {

public:

	Epetra_Operator_Anasazi_Prec(Epetra_Operator *prec, bool precflag){
		prec_ = prec;
		precflag_ = precflag;
	}

    void ApplyPrecondition ( const AnasaziMultiVec<TYPE>& x, AnasaziMultiVec<TYPE>& y ) const;

private:
   Epetra_Operator * prec_;
   bool precflag_;

};
//
// implementation of the Epetra_Operator_Anasazi_Prec class.
//
/////////////////////////////////////////////////////////////
//
// AnasaziPrecond matrix multiply
//
template <class TYPE>
void Epetra_Operator_Anasazi_Prec<TYPE>::ApplyPrecondition ( const AnasaziMultiVec<TYPE>& x, 
						    AnasaziMultiVec<TYPE>& y) const {
  AnasaziMultiVec<TYPE> & temp_x = const_cast<AnasaziMultiVec<TYPE> &>(x);
  Epetra_MultiVector* vec_x = dynamic_cast<Epetra_MultiVector* >(&temp_x);
  Epetra_MultiVector* vec_y = dynamic_cast<Epetra_MultiVector* >(&y);
  
  assert( vec_x || vec_y );

  if (vec_x && vec_y){
	  if (precflag_){
		  //cout << "Applying preconditoner" << endl;
		  assert(prec_->ApplyInverse( *vec_x, *vec_y)==0);
	  }
	  else{
		  //cout << " No preconditioning used" << endl;
		  //cout << " Doing a vector copy" << endl;
		  assert(vec_y->Scale( 1.0, *vec_x)==0);
	  }
  }
  //
}

#endif // EPETRA_OPERATOR_ANASAZI_PREC_H
