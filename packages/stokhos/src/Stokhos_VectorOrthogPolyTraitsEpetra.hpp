// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_VECTORORTHOGPOLYTRAITSEPETRA_HPP
#define STOKHOS_VECTORORTHOGPOLYTRAITSEPETRA_HPP

#include "Stokhos_VectorOrthogPoly.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

namespace Stokhos {

  //! Cloner for Epetra_Vector coefficients
  class EpetraVectorCloner {
  public:
    EpetraVectorCloner(const Epetra_BlockMap& map_) : map(&map_), vec(NULL) {}
    EpetraVectorCloner(const Epetra_Vector& vec_) : map(NULL), vec(&vec_) {}
    Teuchos::RCP<Epetra_Vector> clone(int i) const {
      if (map) 
	return Teuchos::rcp(new Epetra_Vector(*map));
      else 
	return Teuchos::rcp(new Epetra_Vector(*vec));
    }
  protected:
    const Epetra_BlockMap* map;
    const Epetra_Vector* vec;
  };

  //! Cloner for Epetra_Vector coefficients
  /*!
   * Epetra_Operator's cannot be cloned, thus the definition is empty and
   * will cause a compiler error if a clone is attempted.
   */
  class EpetraOperatorCloner {};

  //! Cloner for Epetra_CrsMatrix coefficients
  class EpetraCrsMatrixCloner {
  public:
    EpetraCrsMatrixCloner(const Epetra_CrsMatrix& mat_) : mat(mat_) {}
    Teuchos::RCP<Epetra_CrsMatrix> clone(int i) const {
      return Teuchos::rcp(new Epetra_CrsMatrix(mat));
    }
  protected:
    const Epetra_CrsMatrix& mat;
  };
  
  //! Specialization of VectorOrthogPolyTraits to Epetra_Vector coefficients
  template <> 
  class VectorOrthogPolyTraits<Epetra_Vector> {
  public:
    
    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Typename of cloner
    typedef EpetraVectorCloner cloner_type;

    //! Initialize vector
    static void init(Epetra_Vector& vec, double val) { vec.PutScalar(val); }

    //! Update vector
    static void update(Epetra_Vector& vec, double a, const Epetra_Vector& x) {
      vec.Update(a,x,1.0); 
    }

  };

  //! Specialization of VectorOrthogPolyTraits to Epetra_CrsMatrix coefficients
  template <> 
  class VectorOrthogPolyTraits<Epetra_CrsMatrix> {
  public:
    
    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Typename of cloner
    typedef EpetraCrsMatrixCloner cloner_type;

    //! Initialize matrix
    static void init(Epetra_CrsMatrix& mat, double val) { mat.PutScalar(val); }

    //! Update matrix
    static void update(Epetra_CrsMatrix& mat, double a, 
		       const Epetra_CrsMatrix& x) {
      int num_col;
      for (int i=0; i<mat.NumMyRows(); i++) {
	mat.NumMyRowEntries(i, num_col);
	for (int j=0; j<num_col; j++)
	  mat[i][j] += a*x[i][j];
      }
    }

  };

  //! Specialization of VectorOrthogPolyTraits to Epetra_Operator coefficients
  template <> 
  class VectorOrthogPolyTraits<Epetra_Operator> {
  public:
    
    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Typename of cloner
    typedef EpetraOperatorCloner cloner_type;

    //! Initialize operator
    static void init(Epetra_Operator& op, double val) { 
      Epetra_CrsMatrix& mat = dynamic_cast<Epetra_CrsMatrix&>(op);
      VectorOrthogPolyTraits<Epetra_CrsMatrix>::init(mat, val);
    }

    //! Update operator
    static void update(Epetra_Operator& op, double a, 
		       const Epetra_Operator& x_op) {
      Epetra_CrsMatrix& mat = dynamic_cast<Epetra_CrsMatrix&>(op);
      const Epetra_CrsMatrix& x = dynamic_cast<const Epetra_CrsMatrix&>(x_op);
      VectorOrthogPolyTraits<Epetra_CrsMatrix>::update(mat, a, x);
    }

  };
}

#endif // STOKHOS_VECTORORTHOGPOLYTRAITSEPETRA_HPP
