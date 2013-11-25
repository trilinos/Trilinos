// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_VECTORORTHOGPOLYTRAITSEPETRA_HPP
#define STOKHOS_VECTORORTHOGPOLYTRAITSEPETRA_HPP

#include <iostream>
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "EpetraExt_BlockVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"

namespace Stokhos {

  //! Cloner for Epetra_Vector coefficients
  class EpetraVectorCloner {
  public:
    EpetraVectorCloner(const Epetra_BlockMap& map_) : 
      map(&map_), vec(NULL), block_vec(NULL) {}
    EpetraVectorCloner(const Epetra_Vector& vec_) : 
      map(NULL), vec(&vec_), block_vec(NULL) {}
    EpetraVectorCloner(EpetraExt::BlockVector& block_vec_) : 
      map(NULL), vec(NULL), block_vec(&block_vec_) {}
    Teuchos::RCP<Epetra_Vector> clone(int i) const {
      if (map) 
	return Teuchos::rcp(new Epetra_Vector(*map));
      else if (vec)
	return Teuchos::rcp(new Epetra_Vector(*vec));
      else
	return block_vec->GetBlock(i);
    }
  protected:
    const Epetra_BlockMap* map;
    const Epetra_Vector* vec;
    EpetraExt::BlockVector *block_vec;
    bool view;
  };

  //! Cloner for Epetra_MultiVector coefficients
  class EpetraMultiVectorCloner {
  public:
    EpetraMultiVectorCloner(const Epetra_BlockMap& map_,
			    int num_vectors) : 
      map(&map_), vec(NULL), num_vecs(num_vectors)  {}
    EpetraMultiVectorCloner(const Epetra_MultiVector& vec_) : 
      map(NULL), vec(&vec_), num_vecs(vec_.NumVectors()) {}
    Teuchos::RCP<Epetra_MultiVector> clone(int i) const {
      if (map) 
	return Teuchos::rcp(new Epetra_MultiVector(*map, num_vecs));
      else 
	return Teuchos::rcp(new Epetra_MultiVector(*vec));
    }
  protected:
    const Epetra_BlockMap* map;
    const Epetra_MultiVector* vec;
    int num_vecs;
  };

  //! Cloner for Epetra_Operator coefficients
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
  
  //! Specialization of ProductContainerTraits to Epetra_Vector coefficients
  template <> 
  class ProductContainerTraits<Epetra_Vector> {
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

    //! Print vector
    static std::ostream& print(std::ostream& os, const Epetra_Vector& vec) {
      vec.Print(os);
      return os;
    }

  };

  //! Specialization of ProductContainerTraits to Epetra_MultiVector coefficients
  template <> 
  class ProductContainerTraits<Epetra_MultiVector> {
  public:
    
    //! Typename of values
    typedef double value_type;

    //! Typename of ordinals
    typedef int ordinal_type;

    //! Typename of cloner
    typedef EpetraMultiVectorCloner cloner_type;

    //! Initialize vector
    static void init(Epetra_MultiVector& vec, double val) { 
      vec.PutScalar(val); }

    //! Update vector
    static void update(Epetra_MultiVector& vec, double a, 
		       const Epetra_MultiVector& x) {
      vec.Update(a,x,1.0); 
    }

    //! Print vector
    static std::ostream& print(std::ostream& os, 
			       const Epetra_MultiVector& vec) {
      vec.Print(os);
      return os;
    }

  };

  //! Specialization of ProductContainerTraits to Epetra_CrsMatrix coefficients
  template <> 
  class ProductContainerTraits<Epetra_CrsMatrix> {
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

    //! Print matrix
    static std::ostream& print(std::ostream& os, const Epetra_CrsMatrix& mat) {
      mat.Print(os);
      return os;
    }

  };

  //! Specialization of ProductContainerTraits to Epetra_Operator coefficients
  template <> 
  class ProductContainerTraits<Epetra_Operator> {
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
      ProductContainerTraits<Epetra_CrsMatrix>::init(mat, val);
    }

    //! Update operator
    static void update(Epetra_Operator& op, double a, 
		       const Epetra_Operator& x_op) {
      Epetra_CrsMatrix& mat = dynamic_cast<Epetra_CrsMatrix&>(op);
      const Epetra_CrsMatrix& x = dynamic_cast<const Epetra_CrsMatrix&>(x_op);
      ProductContainerTraits<Epetra_CrsMatrix>::update(mat, a, x);
    }

    //! Print operator
    static std::ostream& print(std::ostream& os, const Epetra_Operator& op) {
      os << "Epetra_Operator" << std::endl;
      const Epetra_CrsMatrix& mat = dynamic_cast<const Epetra_CrsMatrix&>(op);
      ProductContainerTraits<Epetra_CrsMatrix>::print(os, mat);
      return os;
    }

  };

}

#endif // STOKHOS_VECTORORTHOGPOLYTRAITSEPETRA_HPP
