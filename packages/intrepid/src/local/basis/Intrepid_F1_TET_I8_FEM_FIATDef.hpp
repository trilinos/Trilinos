// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov) or
//                    Robert Kirby (robert.c.kirby@ttu.edu)
//                    FIAT (fiat-dev@fenics.org) (must join mailing list first, see www.fenics.org)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_F1_TET_I8_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 7 for 1-forms on TET cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F1_TET_I8_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 1 , 0 , 0 , 8 , 1 , 0 , 1 , 8 , 1 , 0 , 2 , 8 , 1 , 0 , 3 , 8 , 1 , 0 , 4 , 8 , 1 , 0 , 5 , 8 , 1 , 0 , 6 , 8 , 1 , 0 , 7 , 8 , 1 , 1 , 0 , 8 , 1 , 1 , 1 , 8 , 1 , 1 , 2 , 8 , 1 , 1 , 3 , 8 , 1 , 1 , 4 , 8 , 1 , 1 , 5 , 8 , 1 , 1 , 6 , 8 , 1 , 1 , 7 , 8 , 1 , 2 , 0 , 8 , 1 , 2 , 1 , 8 , 1 , 2 , 2 , 8 , 1 , 2 , 3 , 8 , 1 , 2 , 4 , 8 , 1 , 2 , 5 , 8 , 1 , 2 , 6 , 8 , 1 , 2 , 7 , 8 , 1 , 3 , 0 , 8 , 1 , 3 , 1 , 8 , 1 , 3 , 2 , 8 , 1 , 3 , 3 , 8 , 1 , 3 , 4 , 8 , 1 , 3 , 5 , 8 , 1 , 3 , 6 , 8 , 1 , 3 , 7 , 8 , 1 , 4 , 0 , 8 , 1 , 4 , 1 , 8 , 1 , 4 , 2 , 8 , 1 , 4 , 3 , 8 , 1 , 4 , 4 , 8 , 1 , 4 , 5 , 8 , 1 , 4 , 6 , 8 , 1 , 4 , 7 , 8 , 1 , 5 , 0 , 8 , 1 , 5 , 1 , 8 , 1 , 5 , 2 , 8 , 1 , 5 , 3 , 8 , 1 , 5 , 4 , 8 , 1 , 5 , 5 , 8 , 1 , 5 , 6 , 8 , 1 , 5 , 7 , 8 , 2 , 0 , 0 , 56 , 2 , 0 , 1 , 56 , 2 , 0 , 2 , 56 , 2 , 0 , 3 , 56 , 2 , 0 , 4 , 56 , 2 , 0 , 5 , 56 , 2 , 0 , 6 , 56 , 2 , 0 , 7 , 56 , 2 , 0 , 8 , 56 , 2 , 0 , 9 , 56 , 2 , 0 , 10 , 56 , 2 , 0 , 11 , 56 , 2 , 0 , 12 , 56 , 2 , 0 , 13 , 56 , 2 , 0 , 14 , 56 , 2 , 0 , 15 , 56 , 2 , 0 , 16 , 56 , 2 , 0 , 17 , 56 , 2 , 0 , 18 , 56 , 2 , 0 , 19 , 56 , 2 , 0 , 20 , 56 , 2 , 0 , 21 , 56 , 2 , 0 , 22 , 56 , 2 , 0 , 23 , 56 , 2 , 0 , 24 , 56 , 2 , 0 , 25 , 56 , 2 , 0 , 26 , 56 , 2 , 0 , 27 , 56 , 2 , 0 , 28 , 56 , 2 , 0 , 29 , 56 , 2 , 0 , 30 , 56 , 2 , 0 , 31 , 56 , 2 , 0 , 32 , 56 , 2 , 0 , 33 , 56 , 2 , 0 , 34 , 56 , 2 , 0 , 35 , 56 , 2 , 0 , 36 , 56 , 2 , 0 , 37 , 56 , 2 , 0 , 38 , 56 , 2 , 0 , 39 , 56 , 2 , 0 , 40 , 56 , 2 , 0 , 41 , 56 , 2 , 0 , 42 , 56 , 2 , 0 , 43 , 56 , 2 , 0 , 44 , 56 , 2 , 0 , 45 , 56 , 2 , 0 , 46 , 56 , 2 , 0 , 47 , 56 , 2 , 0 , 48 , 56 , 2 , 0 , 49 , 56 , 2 , 0 , 50 , 56 , 2 , 0 , 51 , 56 , 2 , 0 , 52 , 56 , 2 , 0 , 53 , 56 , 2 , 0 , 54 , 56 , 2 , 0 , 55 , 56 , 2 , 1 , 0 , 56 , 2 , 1 , 1 , 56 , 2 , 1 , 2 , 56 , 2 , 1 , 3 , 56 , 2 , 1 , 4 , 56 , 2 , 1 , 5 , 56 , 2 , 1 , 6 , 56 , 2 , 1 , 7 , 56 , 2 , 1 , 8 , 56 , 2 , 1 , 9 , 56 , 2 , 1 , 10 , 56 , 2 , 1 , 11 , 56 , 2 , 1 , 12 , 56 , 2 , 1 , 13 , 56 , 2 , 1 , 14 , 56 , 2 , 1 , 15 , 56 , 2 , 1 , 16 , 56 , 2 , 1 , 17 , 56 , 2 , 1 , 18 , 56 , 2 , 1 , 19 , 56 , 2 , 1 , 20 , 56 , 2 , 1 , 21 , 56 , 2 , 1 , 22 , 56 , 2 , 1 , 23 , 56 , 2 , 1 , 24 , 56 , 2 , 1 , 25 , 56 , 2 , 1 , 26 , 56 , 2 , 1 , 27 , 56 , 2 , 1 , 28 , 56 , 2 , 1 , 29 , 56 , 2 , 1 , 30 , 56 , 2 , 1 , 31 , 56 , 2 , 1 , 32 , 56 , 2 , 1 , 33 , 56 , 2 , 1 , 34 , 56 , 2 , 1 , 35 , 56 , 2 , 1 , 36 , 56 , 2 , 1 , 37 , 56 , 2 , 1 , 38 , 56 , 2 , 1 , 39 , 56 , 2 , 1 , 40 , 56 , 2 , 1 , 41 , 56 , 2 , 1 , 42 , 56 , 2 , 1 , 43 , 56 , 2 , 1 , 44 , 56 , 2 , 1 , 45 , 56 , 2 , 1 , 46 , 56 , 2 , 1 , 47 , 56 , 2 , 1 , 48 , 56 , 2 , 1 , 49 , 56 , 2 , 1 , 50 , 56 , 2 , 1 , 51 , 56 , 2 , 1 , 52 , 56 , 2 , 1 , 53 , 56 , 2 , 1 , 54 , 56 , 2 , 1 , 55 , 56 , 2 , 2 , 0 , 56 , 2 , 2 , 1 , 56 , 2 , 2 , 2 , 56 , 2 , 2 , 3 , 56 , 2 , 2 , 4 , 56 , 2 , 2 , 5 , 56 , 2 , 2 , 6 , 56 , 2 , 2 , 7 , 56 , 2 , 2 , 8 , 56 , 2 , 2 , 9 , 56 , 2 , 2 , 10 , 56 , 2 , 2 , 11 , 56 , 2 , 2 , 12 , 56 , 2 , 2 , 13 , 56 , 2 , 2 , 14 , 56 , 2 , 2 , 15 , 56 , 2 , 2 , 16 , 56 , 2 , 2 , 17 , 56 , 2 , 2 , 18 , 56 , 2 , 2 , 19 , 56 , 2 , 2 , 20 , 56 , 2 , 2 , 21 , 56 , 2 , 2 , 22 , 56 , 2 , 2 , 23 , 56 , 2 , 2 , 24 , 56 , 2 , 2 , 25 , 56 , 2 , 2 , 26 , 56 , 2 , 2 , 27 , 56 , 2 , 2 , 28 , 56 , 2 , 2 , 29 , 56 , 2 , 2 , 30 , 56 , 2 , 2 , 31 , 56 , 2 , 2 , 32 , 56 , 2 , 2 , 33 , 56 , 2 , 2 , 34 , 56 , 2 , 2 , 35 , 56 , 2 , 2 , 36 , 56 , 2 , 2 , 37 , 56 , 2 , 2 , 38 , 56 , 2 , 2 , 39 , 56 , 2 , 2 , 40 , 56 , 2 , 2 , 41 , 56 , 2 , 2 , 42 , 56 , 2 , 2 , 43 , 56 , 2 , 2 , 44 , 56 , 2 , 2 , 45 , 56 , 2 , 2 , 46 , 56 , 2 , 2 , 47 , 56 , 2 , 2 , 48 , 56 , 2 , 2 , 49 , 56 , 2 , 2 , 50 , 56 , 2 , 2 , 51 , 56 , 2 , 2 , 52 , 56 , 2 , 2 , 53 , 56 , 2 , 2 , 54 , 56 , 2 , 2 , 55 , 56 , 2 , 3 , 0 , 56 , 2 , 3 , 1 , 56 , 2 , 3 , 2 , 56 , 2 , 3 , 3 , 56 , 2 , 3 , 4 , 56 , 2 , 3 , 5 , 56 , 2 , 3 , 6 , 56 , 2 , 3 , 7 , 56 , 2 , 3 , 8 , 56 , 2 , 3 , 9 , 56 , 2 , 3 , 10 , 56 , 2 , 3 , 11 , 56 , 2 , 3 , 12 , 56 , 2 , 3 , 13 , 56 , 2 , 3 , 14 , 56 , 2 , 3 , 15 , 56 , 2 , 3 , 16 , 56 , 2 , 3 , 17 , 56 , 2 , 3 , 18 , 56 , 2 , 3 , 19 , 56 , 2 , 3 , 20 , 56 , 2 , 3 , 21 , 56 , 2 , 3 , 22 , 56 , 2 , 3 , 23 , 56 , 2 , 3 , 24 , 56 , 2 , 3 , 25 , 56 , 2 , 3 , 26 , 56 , 2 , 3 , 27 , 56 , 2 , 3 , 28 , 56 , 2 , 3 , 29 , 56 , 2 , 3 , 30 , 56 , 2 , 3 , 31 , 56 , 2 , 3 , 32 , 56 , 2 , 3 , 33 , 56 , 2 , 3 , 34 , 56 , 2 , 3 , 35 , 56 , 2 , 3 , 36 , 56 , 2 , 3 , 37 , 56 , 2 , 3 , 38 , 56 , 2 , 3 , 39 , 56 , 2 , 3 , 40 , 56 , 2 , 3 , 41 , 56 , 2 , 3 , 42 , 56 , 2 , 3 , 43 , 56 , 2 , 3 , 44 , 56 , 2 , 3 , 45 , 56 , 2 , 3 , 46 , 56 , 2 , 3 , 47 , 56 , 2 , 3 , 48 , 56 , 2 , 3 , 49 , 56 , 2 , 3 , 50 , 56 , 2 , 3 , 51 , 56 , 2 , 3 , 52 , 56 , 2 , 3 , 53 , 56 , 2 , 3 , 54 , 56 , 2 , 3 , 55 , 56 , 3 , 0 , 0 , 168 , 3 , 0 , 1 , 168 , 3 , 0 , 2 , 168 , 3 , 0 , 3 , 168 , 3 , 0 , 4 , 168 , 3 , 0 , 5 , 168 , 3 , 0 , 6 , 168 , 3 , 0 , 7 , 168 , 3 , 0 , 8 , 168 , 3 , 0 , 9 , 168 , 3 , 0 , 10 , 168 , 3 , 0 , 11 , 168 , 3 , 0 , 12 , 168 , 3 , 0 , 13 , 168 , 3 , 0 , 14 , 168 , 3 , 0 , 15 , 168 , 3 , 0 , 16 , 168 , 3 , 0 , 17 , 168 , 3 , 0 , 18 , 168 , 3 , 0 , 19 , 168 , 3 , 0 , 20 , 168 , 3 , 0 , 21 , 168 , 3 , 0 , 22 , 168 , 3 , 0 , 23 , 168 , 3 , 0 , 24 , 168 , 3 , 0 , 25 , 168 , 3 , 0 , 26 , 168 , 3 , 0 , 27 , 168 , 3 , 0 , 28 , 168 , 3 , 0 , 29 , 168 , 3 , 0 , 30 , 168 , 3 , 0 , 31 , 168 , 3 , 0 , 32 , 168 , 3 , 0 , 33 , 168 , 3 , 0 , 34 , 168 , 3 , 0 , 35 , 168 , 3 , 0 , 36 , 168 , 3 , 0 , 37 , 168 , 3 , 0 , 38 , 168 , 3 , 0 , 39 , 168 , 3 , 0 , 40 , 168 , 3 , 0 , 41 , 168 , 3 , 0 , 42 , 168 , 3 , 0 , 43 , 168 , 3 , 0 , 44 , 168 , 3 , 0 , 45 , 168 , 3 , 0 , 46 , 168 , 3 , 0 , 47 , 168 , 3 , 0 , 48 , 168 , 3 , 0 , 49 , 168 , 3 , 0 , 50 , 168 , 3 , 0 , 51 , 168 , 3 , 0 , 52 , 168 , 3 , 0 , 53 , 168 , 3 , 0 , 54 , 168 , 3 , 0 , 55 , 168 , 3 , 0 , 56 , 168 , 3 , 0 , 57 , 168 , 3 , 0 , 58 , 168 , 3 , 0 , 59 , 168 , 3 , 0 , 60 , 168 , 3 , 0 , 61 , 168 , 3 , 0 , 62 , 168 , 3 , 0 , 63 , 168 , 3 , 0 , 64 , 168 , 3 , 0 , 65 , 168 , 3 , 0 , 66 , 168 , 3 , 0 , 67 , 168 , 3 , 0 , 68 , 168 , 3 , 0 , 69 , 168 , 3 , 0 , 70 , 168 , 3 , 0 , 71 , 168 , 3 , 0 , 72 , 168 , 3 , 0 , 73 , 168 , 3 , 0 , 74 , 168 , 3 , 0 , 75 , 168 , 3 , 0 , 76 , 168 , 3 , 0 , 77 , 168 , 3 , 0 , 78 , 168 , 3 , 0 , 79 , 168 , 3 , 0 , 80 , 168 , 3 , 0 , 81 , 168 , 3 , 0 , 82 , 168 , 3 , 0 , 83 , 168 , 3 , 0 , 84 , 168 , 3 , 0 , 85 , 168 , 3 , 0 , 86 , 168 , 3 , 0 , 87 , 168 , 3 , 0 , 88 , 168 , 3 , 0 , 89 , 168 , 3 , 0 , 90 , 168 , 3 , 0 , 91 , 168 , 3 , 0 , 92 , 168 , 3 , 0 , 93 , 168 , 3 , 0 , 94 , 168 , 3 , 0 , 95 , 168 , 3 , 0 , 96 , 168 , 3 , 0 , 97 , 168 , 3 , 0 , 98 , 168 , 3 , 0 , 99 , 168 , 3 , 0 , 100 , 168 , 3 , 0 , 101 , 168 , 3 , 0 , 102 , 168 , 3 , 0 , 103 , 168 , 3 , 0 , 104 , 168 , 3 , 0 , 105 , 168 , 3 , 0 , 106 , 168 , 3 , 0 , 107 , 168 , 3 , 0 , 108 , 168 , 3 , 0 , 109 , 168 , 3 , 0 , 110 , 168 , 3 , 0 , 111 , 168 , 3 , 0 , 112 , 168 , 3 , 0 , 113 , 168 , 3 , 0 , 114 , 168 , 3 , 0 , 115 , 168 , 3 , 0 , 116 , 168 , 3 , 0 , 117 , 168 , 3 , 0 , 118 , 168 , 3 , 0 , 119 , 168 , 3 , 0 , 120 , 168 , 3 , 0 , 121 , 168 , 3 , 0 , 122 , 168 , 3 , 0 , 123 , 168 , 3 , 0 , 124 , 168 , 3 , 0 , 125 , 168 , 3 , 0 , 126 , 168 , 3 , 0 , 127 , 168 , 3 , 0 , 128 , 168 , 3 , 0 , 129 , 168 , 3 , 0 , 130 , 168 , 3 , 0 , 131 , 168 , 3 , 0 , 132 , 168 , 3 , 0 , 133 , 168 , 3 , 0 , 134 , 168 , 3 , 0 , 135 , 168 , 3 , 0 , 136 , 168 , 3 , 0 , 137 , 168 , 3 , 0 , 138 , 168 , 3 , 0 , 139 , 168 , 3 , 0 , 140 , 168 , 3 , 0 , 141 , 168 , 3 , 0 , 142 , 168 , 3 , 0 , 143 , 168 , 3 , 0 , 144 , 168 , 3 , 0 , 145 , 168 , 3 , 0 , 146 , 168 , 3 , 0 , 147 , 168 , 3 , 0 , 148 , 168 , 3 , 0 , 149 , 168 , 3 , 0 , 150 , 168 , 3 , 0 , 151 , 168 , 3 , 0 , 152 , 168 , 3 , 0 , 153 , 168 , 3 , 0 , 154 , 168 , 3 , 0 , 155 , 168 , 3 , 0 , 156 , 168 , 3 , 0 , 157 , 168 , 3 , 0 , 158 , 168 , 3 , 0 , 159 , 168 , 3 , 0 , 160 , 168 , 3 , 0 , 161 , 168 , 3 , 0 , 162 , 168 , 3 , 0 , 163 , 168 , 3 , 0 , 164 , 168 , 3 , 0 , 165 , 168 , 3 , 0 , 166 , 168 , 3 , 0 , 167 , 168 };
  
  // Basis-independent function sets tag and enum data in the static arrays:
  Intrepid::setEnumTagData(tagToEnum_,
                           enumToTag_,
                           tags,
                           numDof_,
                           tagSize,
                           posScDim,
                           posScId,
                           posBfId);
}

template<class Scalar> 
void Basis_F1_TET_I8_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 7 (I7) has 440 basis functions on a tetrahedron that are 1-forms in 3D
  int    numFields = 440;
  EField fieldType = FIELD_FORM_1;
  int    spaceDim  = 3;

  // temporaries
  int countPt  = 0;               // point counter
  Teuchos::Array<int> indexV(3);  // multi-index for values
  Teuchos::Array<int> indexD(3);  // multi-index for curl

  // Shape the FieldContainer for the output values using these values:
  outputValues.resize(numPoints,
                     numFields,
                     fieldType,
                     operatorType,
                     spaceDim);

#ifdef HAVE_INTREPID_DEBUG
  for (countPt=0; countPt<numPoints; countPt++) {
    // Verify that all points are inside the TET reference cell
    TEST_FOR_EXCEPTION( !MultiCell<Scalar>::inReferenceCell(CELL_TET, inputPoints[countPt]),
                        std::invalid_argument,
                        ">>> ERROR (Basis_F1_TET_I8_FEM_FIAT): Evaluation point is outside the TET reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(165,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(440,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,8,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<440;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<440;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm2_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<440;i++) {
          outputValues(countPt,i,2) = result(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TET_I8_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(165,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(165,165);
      Teuchos::SerialDenseMatrix<int,Scalar> curlcomp(440,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,8,inputPoints,expansions);
      curlcomp.putScalar(0.0);
      // x component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm2_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm1_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<440;i++) {
          outputValues(countPt,i,0) = curlcomp(i,countPt);
        }
      }
      // y component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm2_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<440;i++) {
          outputValues(countPt,i,1) = curlcomp(i,countPt);
        }
      }
      // z component of curl
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,0.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      curlcomp.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,-1.0,*vdm0_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<440;i++) {
          outputValues(countPt,i,2) = curlcomp(i,countPt);
        }
      }
    }
    break;

    case OPERATOR_DIV:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F1_TET_I8_FEM_FIAT: operator not implemented" );
    }
    break; 
    default:
      TEST_FOR_EXCEPTION( ( (operatorType != OPERATOR_VALUE) &&
                            (operatorType != OPERATOR_GRAD)  &&
                            (operatorType != OPERATOR_CURL)  &&
                            (operatorType != OPERATOR_DIV)   &&
                            (operatorType != OPERATOR_D1)    &&
                            (operatorType != OPERATOR_D2)    &&
                            (operatorType != OPERATOR_D3)    &&
                            (operatorType != OPERATOR_D4)    &&
                            (operatorType != OPERATOR_D5)    &&
                            (operatorType != OPERATOR_D6)    &&
                            (operatorType != OPERATOR_D7)    &&
                            (operatorType != OPERATOR_D8)    &&
                            (operatorType != OPERATOR_D9)    &&
                            (operatorType != OPERATOR_D10) ),
                          std::invalid_argument,
                          ">>> ERROR (Basis_F1_TET_I8_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F1_TET_I8_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F1_TET_I8_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F1_TET_I8_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F1_TET_I8_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F1_TET_I8_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F1_TET_I8_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F1_TET_I8_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TET;
}



template<class Scalar>
EBasis Basis_F1_TET_I8_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F1_TET_I8_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F1_TET_I8_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

