//@HEADER
/*
************************************************************************

              Komplex: Complex Linear Solver Package 
                Copyright (2002) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

//! Komplex_Ordering: The Komplex Ordering Class.
/*! The Komplex_Ordering class defines the ordering for a Komplex_RowMatrix object.

The different K forms (K1, K2, K3, K4, K14, and K23) of an equivalent real formulation can 
easily convert back and forth by going from one K form to the canonical form to another K 
form.  The Komplex_Ordering that each Komplex_RowMatrix object has is what determines the 
conversions.  Let Kanon stand for the canonical form of a complex matrix in equivalent real 
formulation.  Then any K form is equivalent to:
           P_l * Kanon * P_r * D_r,
where P_l, P_r are specific permutation matrices and D_r is a specific right diagonal scaling matrix.
This is helpful because certain forms are advantageous in certain conditions.  To be able to convert
back and forth during preconditioning and solving should allow for faster, more accurate solutions.
*/

class Komplex_Ordering {

  public:
  //@{ \name Constructors/Destructor.
  //! Komplex_Ordering Default Constructor.
  /*! Creates an empty Komplex_Ordering instance.
  */
  Komplex_Ordering();

  //! Komplex_Ordering Constructor to create the ordering for a Komplex_RowMatrix object.
  /*! Creates a Komplex_Ordering instance for a given Komplex_RowMatrix object.
  */
  Komplex_Ordering(Komplex_RowMatrix * A);

  //! Komplex_Ordering Copy Constructor.
  /*! Makes a copy of an existing Komplex_Ordering instance.
  */
  Komplex_Ordering(const Komplex_Ordering & Ordering);

  //! Komplex_Ordering Destructor.
  /*! Completely deletes a Komplex_Ordering object.
  */
  virtual ~Komplex_Ordering();
  //@}

  //@{ \name Matrix converter methods

  //! Convert a Komplex_RowMatrix to canonical form.
  /*! Converts a given Komplex_RowMatrix into canonical form, returning the answer
      in a user-provided RowMatrix object.
      \param In 
             KForm - RowMatrix to be converted to canonical form.
      \param Out
             Canonical - RowMatrix in canonical form.
   	\return Integer error code, set to 0 if successful.
  */
  int ToCanonical(Komplex_RowMatrix * KForm, Komplex_RowMatrix & Canonical);

  //! Convert a Komplex_RowMatrix in canonical form to K form.
  /*! Converts a given Komplex_RowMatrix in canonical form to K form, returning the answer
      in a user-provided RowMatrix object.
 	\param In
		 Canonical - RowMatrix to be converted to K form.
 	\param Out
		 KForm - RowMatrix in K form.
	\return Integer error code, set to 0 if successful.
  */
  int ToKForm(Komplex_RowMatrix * Canonical, Komplex_RowMatrix & KForm);

  private:
  Epetra_MultiVector * DiagRight_;
  Epetra_MultiVector * PermLeft_;
  Epetra_MultiVector * PermRight_;
};