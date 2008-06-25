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

/** \file   Intrepid_F2_TET_I10_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 9 for 2-forms on TET cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F2_TET_I10_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 2 , 0 , 0 , 55 , 2 , 0 , 1 , 55 , 2 , 0 , 2 , 55 , 2 , 0 , 3 , 55 , 2 , 0 , 4 , 55 , 2 , 0 , 5 , 55 , 2 , 0 , 6 , 55 , 2 , 0 , 7 , 55 , 2 , 0 , 8 , 55 , 2 , 0 , 9 , 55 , 2 , 0 , 10 , 55 , 2 , 0 , 11 , 55 , 2 , 0 , 12 , 55 , 2 , 0 , 13 , 55 , 2 , 0 , 14 , 55 , 2 , 0 , 15 , 55 , 2 , 0 , 16 , 55 , 2 , 0 , 17 , 55 , 2 , 0 , 18 , 55 , 2 , 0 , 19 , 55 , 2 , 0 , 20 , 55 , 2 , 0 , 21 , 55 , 2 , 0 , 22 , 55 , 2 , 0 , 23 , 55 , 2 , 0 , 24 , 55 , 2 , 0 , 25 , 55 , 2 , 0 , 26 , 55 , 2 , 0 , 27 , 55 , 2 , 0 , 28 , 55 , 2 , 0 , 29 , 55 , 2 , 0 , 30 , 55 , 2 , 0 , 31 , 55 , 2 , 0 , 32 , 55 , 2 , 0 , 33 , 55 , 2 , 0 , 34 , 55 , 2 , 0 , 35 , 55 , 2 , 0 , 36 , 55 , 2 , 0 , 37 , 55 , 2 , 0 , 38 , 55 , 2 , 0 , 39 , 55 , 2 , 0 , 40 , 55 , 2 , 0 , 41 , 55 , 2 , 0 , 42 , 55 , 2 , 0 , 43 , 55 , 2 , 0 , 44 , 55 , 2 , 0 , 45 , 55 , 2 , 0 , 46 , 55 , 2 , 0 , 47 , 55 , 2 , 0 , 48 , 55 , 2 , 0 , 49 , 55 , 2 , 0 , 50 , 55 , 2 , 0 , 51 , 55 , 2 , 0 , 52 , 55 , 2 , 0 , 53 , 55 , 2 , 0 , 54 , 55 , 2 , 1 , 0 , 55 , 2 , 1 , 1 , 55 , 2 , 1 , 2 , 55 , 2 , 1 , 3 , 55 , 2 , 1 , 4 , 55 , 2 , 1 , 5 , 55 , 2 , 1 , 6 , 55 , 2 , 1 , 7 , 55 , 2 , 1 , 8 , 55 , 2 , 1 , 9 , 55 , 2 , 1 , 10 , 55 , 2 , 1 , 11 , 55 , 2 , 1 , 12 , 55 , 2 , 1 , 13 , 55 , 2 , 1 , 14 , 55 , 2 , 1 , 15 , 55 , 2 , 1 , 16 , 55 , 2 , 1 , 17 , 55 , 2 , 1 , 18 , 55 , 2 , 1 , 19 , 55 , 2 , 1 , 20 , 55 , 2 , 1 , 21 , 55 , 2 , 1 , 22 , 55 , 2 , 1 , 23 , 55 , 2 , 1 , 24 , 55 , 2 , 1 , 25 , 55 , 2 , 1 , 26 , 55 , 2 , 1 , 27 , 55 , 2 , 1 , 28 , 55 , 2 , 1 , 29 , 55 , 2 , 1 , 30 , 55 , 2 , 1 , 31 , 55 , 2 , 1 , 32 , 55 , 2 , 1 , 33 , 55 , 2 , 1 , 34 , 55 , 2 , 1 , 35 , 55 , 2 , 1 , 36 , 55 , 2 , 1 , 37 , 55 , 2 , 1 , 38 , 55 , 2 , 1 , 39 , 55 , 2 , 1 , 40 , 55 , 2 , 1 , 41 , 55 , 2 , 1 , 42 , 55 , 2 , 1 , 43 , 55 , 2 , 1 , 44 , 55 , 2 , 1 , 45 , 55 , 2 , 1 , 46 , 55 , 2 , 1 , 47 , 55 , 2 , 1 , 48 , 55 , 2 , 1 , 49 , 55 , 2 , 1 , 50 , 55 , 2 , 1 , 51 , 55 , 2 , 1 , 52 , 55 , 2 , 1 , 53 , 55 , 2 , 1 , 54 , 55 , 2 , 2 , 0 , 55 , 2 , 2 , 1 , 55 , 2 , 2 , 2 , 55 , 2 , 2 , 3 , 55 , 2 , 2 , 4 , 55 , 2 , 2 , 5 , 55 , 2 , 2 , 6 , 55 , 2 , 2 , 7 , 55 , 2 , 2 , 8 , 55 , 2 , 2 , 9 , 55 , 2 , 2 , 10 , 55 , 2 , 2 , 11 , 55 , 2 , 2 , 12 , 55 , 2 , 2 , 13 , 55 , 2 , 2 , 14 , 55 , 2 , 2 , 15 , 55 , 2 , 2 , 16 , 55 , 2 , 2 , 17 , 55 , 2 , 2 , 18 , 55 , 2 , 2 , 19 , 55 , 2 , 2 , 20 , 55 , 2 , 2 , 21 , 55 , 2 , 2 , 22 , 55 , 2 , 2 , 23 , 55 , 2 , 2 , 24 , 55 , 2 , 2 , 25 , 55 , 2 , 2 , 26 , 55 , 2 , 2 , 27 , 55 , 2 , 2 , 28 , 55 , 2 , 2 , 29 , 55 , 2 , 2 , 30 , 55 , 2 , 2 , 31 , 55 , 2 , 2 , 32 , 55 , 2 , 2 , 33 , 55 , 2 , 2 , 34 , 55 , 2 , 2 , 35 , 55 , 2 , 2 , 36 , 55 , 2 , 2 , 37 , 55 , 2 , 2 , 38 , 55 , 2 , 2 , 39 , 55 , 2 , 2 , 40 , 55 , 2 , 2 , 41 , 55 , 2 , 2 , 42 , 55 , 2 , 2 , 43 , 55 , 2 , 2 , 44 , 55 , 2 , 2 , 45 , 55 , 2 , 2 , 46 , 55 , 2 , 2 , 47 , 55 , 2 , 2 , 48 , 55 , 2 , 2 , 49 , 55 , 2 , 2 , 50 , 55 , 2 , 2 , 51 , 55 , 2 , 2 , 52 , 55 , 2 , 2 , 53 , 55 , 2 , 2 , 54 , 55 , 2 , 3 , 0 , 55 , 2 , 3 , 1 , 55 , 2 , 3 , 2 , 55 , 2 , 3 , 3 , 55 , 2 , 3 , 4 , 55 , 2 , 3 , 5 , 55 , 2 , 3 , 6 , 55 , 2 , 3 , 7 , 55 , 2 , 3 , 8 , 55 , 2 , 3 , 9 , 55 , 2 , 3 , 10 , 55 , 2 , 3 , 11 , 55 , 2 , 3 , 12 , 55 , 2 , 3 , 13 , 55 , 2 , 3 , 14 , 55 , 2 , 3 , 15 , 55 , 2 , 3 , 16 , 55 , 2 , 3 , 17 , 55 , 2 , 3 , 18 , 55 , 2 , 3 , 19 , 55 , 2 , 3 , 20 , 55 , 2 , 3 , 21 , 55 , 2 , 3 , 22 , 55 , 2 , 3 , 23 , 55 , 2 , 3 , 24 , 55 , 2 , 3 , 25 , 55 , 2 , 3 , 26 , 55 , 2 , 3 , 27 , 55 , 2 , 3 , 28 , 55 , 2 , 3 , 29 , 55 , 2 , 3 , 30 , 55 , 2 , 3 , 31 , 55 , 2 , 3 , 32 , 55 , 2 , 3 , 33 , 55 , 2 , 3 , 34 , 55 , 2 , 3 , 35 , 55 , 2 , 3 , 36 , 55 , 2 , 3 , 37 , 55 , 2 , 3 , 38 , 55 , 2 , 3 , 39 , 55 , 2 , 3 , 40 , 55 , 2 , 3 , 41 , 55 , 2 , 3 , 42 , 55 , 2 , 3 , 43 , 55 , 2 , 3 , 44 , 55 , 2 , 3 , 45 , 55 , 2 , 3 , 46 , 55 , 2 , 3 , 47 , 55 , 2 , 3 , 48 , 55 , 2 , 3 , 49 , 55 , 2 , 3 , 50 , 55 , 2 , 3 , 51 , 55 , 2 , 3 , 52 , 55 , 2 , 3 , 53 , 55 , 2 , 3 , 54 , 55 , 3 , 0 , 0 , 495 , 3 , 0 , 1 , 495 , 3 , 0 , 2 , 495 , 3 , 0 , 3 , 495 , 3 , 0 , 4 , 495 , 3 , 0 , 5 , 495 , 3 , 0 , 6 , 495 , 3 , 0 , 7 , 495 , 3 , 0 , 8 , 495 , 3 , 0 , 9 , 495 , 3 , 0 , 10 , 495 , 3 , 0 , 11 , 495 , 3 , 0 , 12 , 495 , 3 , 0 , 13 , 495 , 3 , 0 , 14 , 495 , 3 , 0 , 15 , 495 , 3 , 0 , 16 , 495 , 3 , 0 , 17 , 495 , 3 , 0 , 18 , 495 , 3 , 0 , 19 , 495 , 3 , 0 , 20 , 495 , 3 , 0 , 21 , 495 , 3 , 0 , 22 , 495 , 3 , 0 , 23 , 495 , 3 , 0 , 24 , 495 , 3 , 0 , 25 , 495 , 3 , 0 , 26 , 495 , 3 , 0 , 27 , 495 , 3 , 0 , 28 , 495 , 3 , 0 , 29 , 495 , 3 , 0 , 30 , 495 , 3 , 0 , 31 , 495 , 3 , 0 , 32 , 495 , 3 , 0 , 33 , 495 , 3 , 0 , 34 , 495 , 3 , 0 , 35 , 495 , 3 , 0 , 36 , 495 , 3 , 0 , 37 , 495 , 3 , 0 , 38 , 495 , 3 , 0 , 39 , 495 , 3 , 0 , 40 , 495 , 3 , 0 , 41 , 495 , 3 , 0 , 42 , 495 , 3 , 0 , 43 , 495 , 3 , 0 , 44 , 495 , 3 , 0 , 45 , 495 , 3 , 0 , 46 , 495 , 3 , 0 , 47 , 495 , 3 , 0 , 48 , 495 , 3 , 0 , 49 , 495 , 3 , 0 , 50 , 495 , 3 , 0 , 51 , 495 , 3 , 0 , 52 , 495 , 3 , 0 , 53 , 495 , 3 , 0 , 54 , 495 , 3 , 0 , 55 , 495 , 3 , 0 , 56 , 495 , 3 , 0 , 57 , 495 , 3 , 0 , 58 , 495 , 3 , 0 , 59 , 495 , 3 , 0 , 60 , 495 , 3 , 0 , 61 , 495 , 3 , 0 , 62 , 495 , 3 , 0 , 63 , 495 , 3 , 0 , 64 , 495 , 3 , 0 , 65 , 495 , 3 , 0 , 66 , 495 , 3 , 0 , 67 , 495 , 3 , 0 , 68 , 495 , 3 , 0 , 69 , 495 , 3 , 0 , 70 , 495 , 3 , 0 , 71 , 495 , 3 , 0 , 72 , 495 , 3 , 0 , 73 , 495 , 3 , 0 , 74 , 495 , 3 , 0 , 75 , 495 , 3 , 0 , 76 , 495 , 3 , 0 , 77 , 495 , 3 , 0 , 78 , 495 , 3 , 0 , 79 , 495 , 3 , 0 , 80 , 495 , 3 , 0 , 81 , 495 , 3 , 0 , 82 , 495 , 3 , 0 , 83 , 495 , 3 , 0 , 84 , 495 , 3 , 0 , 85 , 495 , 3 , 0 , 86 , 495 , 3 , 0 , 87 , 495 , 3 , 0 , 88 , 495 , 3 , 0 , 89 , 495 , 3 , 0 , 90 , 495 , 3 , 0 , 91 , 495 , 3 , 0 , 92 , 495 , 3 , 0 , 93 , 495 , 3 , 0 , 94 , 495 , 3 , 0 , 95 , 495 , 3 , 0 , 96 , 495 , 3 , 0 , 97 , 495 , 3 , 0 , 98 , 495 , 3 , 0 , 99 , 495 , 3 , 0 , 100 , 495 , 3 , 0 , 101 , 495 , 3 , 0 , 102 , 495 , 3 , 0 , 103 , 495 , 3 , 0 , 104 , 495 , 3 , 0 , 105 , 495 , 3 , 0 , 106 , 495 , 3 , 0 , 107 , 495 , 3 , 0 , 108 , 495 , 3 , 0 , 109 , 495 , 3 , 0 , 110 , 495 , 3 , 0 , 111 , 495 , 3 , 0 , 112 , 495 , 3 , 0 , 113 , 495 , 3 , 0 , 114 , 495 , 3 , 0 , 115 , 495 , 3 , 0 , 116 , 495 , 3 , 0 , 117 , 495 , 3 , 0 , 118 , 495 , 3 , 0 , 119 , 495 , 3 , 0 , 120 , 495 , 3 , 0 , 121 , 495 , 3 , 0 , 122 , 495 , 3 , 0 , 123 , 495 , 3 , 0 , 124 , 495 , 3 , 0 , 125 , 495 , 3 , 0 , 126 , 495 , 3 , 0 , 127 , 495 , 3 , 0 , 128 , 495 , 3 , 0 , 129 , 495 , 3 , 0 , 130 , 495 , 3 , 0 , 131 , 495 , 3 , 0 , 132 , 495 , 3 , 0 , 133 , 495 , 3 , 0 , 134 , 495 , 3 , 0 , 135 , 495 , 3 , 0 , 136 , 495 , 3 , 0 , 137 , 495 , 3 , 0 , 138 , 495 , 3 , 0 , 139 , 495 , 3 , 0 , 140 , 495 , 3 , 0 , 141 , 495 , 3 , 0 , 142 , 495 , 3 , 0 , 143 , 495 , 3 , 0 , 144 , 495 , 3 , 0 , 145 , 495 , 3 , 0 , 146 , 495 , 3 , 0 , 147 , 495 , 3 , 0 , 148 , 495 , 3 , 0 , 149 , 495 , 3 , 0 , 150 , 495 , 3 , 0 , 151 , 495 , 3 , 0 , 152 , 495 , 3 , 0 , 153 , 495 , 3 , 0 , 154 , 495 , 3 , 0 , 155 , 495 , 3 , 0 , 156 , 495 , 3 , 0 , 157 , 495 , 3 , 0 , 158 , 495 , 3 , 0 , 159 , 495 , 3 , 0 , 160 , 495 , 3 , 0 , 161 , 495 , 3 , 0 , 162 , 495 , 3 , 0 , 163 , 495 , 3 , 0 , 164 , 495 , 3 , 0 , 165 , 495 , 3 , 0 , 166 , 495 , 3 , 0 , 167 , 495 , 3 , 0 , 168 , 495 , 3 , 0 , 169 , 495 , 3 , 0 , 170 , 495 , 3 , 0 , 171 , 495 , 3 , 0 , 172 , 495 , 3 , 0 , 173 , 495 , 3 , 0 , 174 , 495 , 3 , 0 , 175 , 495 , 3 , 0 , 176 , 495 , 3 , 0 , 177 , 495 , 3 , 0 , 178 , 495 , 3 , 0 , 179 , 495 , 3 , 0 , 180 , 495 , 3 , 0 , 181 , 495 , 3 , 0 , 182 , 495 , 3 , 0 , 183 , 495 , 3 , 0 , 184 , 495 , 3 , 0 , 185 , 495 , 3 , 0 , 186 , 495 , 3 , 0 , 187 , 495 , 3 , 0 , 188 , 495 , 3 , 0 , 189 , 495 , 3 , 0 , 190 , 495 , 3 , 0 , 191 , 495 , 3 , 0 , 192 , 495 , 3 , 0 , 193 , 495 , 3 , 0 , 194 , 495 , 3 , 0 , 195 , 495 , 3 , 0 , 196 , 495 , 3 , 0 , 197 , 495 , 3 , 0 , 198 , 495 , 3 , 0 , 199 , 495 , 3 , 0 , 200 , 495 , 3 , 0 , 201 , 495 , 3 , 0 , 202 , 495 , 3 , 0 , 203 , 495 , 3 , 0 , 204 , 495 , 3 , 0 , 205 , 495 , 3 , 0 , 206 , 495 , 3 , 0 , 207 , 495 , 3 , 0 , 208 , 495 , 3 , 0 , 209 , 495 , 3 , 0 , 210 , 495 , 3 , 0 , 211 , 495 , 3 , 0 , 212 , 495 , 3 , 0 , 213 , 495 , 3 , 0 , 214 , 495 , 3 , 0 , 215 , 495 , 3 , 0 , 216 , 495 , 3 , 0 , 217 , 495 , 3 , 0 , 218 , 495 , 3 , 0 , 219 , 495 , 3 , 0 , 220 , 495 , 3 , 0 , 221 , 495 , 3 , 0 , 222 , 495 , 3 , 0 , 223 , 495 , 3 , 0 , 224 , 495 , 3 , 0 , 225 , 495 , 3 , 0 , 226 , 495 , 3 , 0 , 227 , 495 , 3 , 0 , 228 , 495 , 3 , 0 , 229 , 495 , 3 , 0 , 230 , 495 , 3 , 0 , 231 , 495 , 3 , 0 , 232 , 495 , 3 , 0 , 233 , 495 , 3 , 0 , 234 , 495 , 3 , 0 , 235 , 495 , 3 , 0 , 236 , 495 , 3 , 0 , 237 , 495 , 3 , 0 , 238 , 495 , 3 , 0 , 239 , 495 , 3 , 0 , 240 , 495 , 3 , 0 , 241 , 495 , 3 , 0 , 242 , 495 , 3 , 0 , 243 , 495 , 3 , 0 , 244 , 495 , 3 , 0 , 245 , 495 , 3 , 0 , 246 , 495 , 3 , 0 , 247 , 495 , 3 , 0 , 248 , 495 , 3 , 0 , 249 , 495 , 3 , 0 , 250 , 495 , 3 , 0 , 251 , 495 , 3 , 0 , 252 , 495 , 3 , 0 , 253 , 495 , 3 , 0 , 254 , 495 , 3 , 0 , 255 , 495 , 3 , 0 , 256 , 495 , 3 , 0 , 257 , 495 , 3 , 0 , 258 , 495 , 3 , 0 , 259 , 495 , 3 , 0 , 260 , 495 , 3 , 0 , 261 , 495 , 3 , 0 , 262 , 495 , 3 , 0 , 263 , 495 , 3 , 0 , 264 , 495 , 3 , 0 , 265 , 495 , 3 , 0 , 266 , 495 , 3 , 0 , 267 , 495 , 3 , 0 , 268 , 495 , 3 , 0 , 269 , 495 , 3 , 0 , 270 , 495 , 3 , 0 , 271 , 495 , 3 , 0 , 272 , 495 , 3 , 0 , 273 , 495 , 3 , 0 , 274 , 495 , 3 , 0 , 275 , 495 , 3 , 0 , 276 , 495 , 3 , 0 , 277 , 495 , 3 , 0 , 278 , 495 , 3 , 0 , 279 , 495 , 3 , 0 , 280 , 495 , 3 , 0 , 281 , 495 , 3 , 0 , 282 , 495 , 3 , 0 , 283 , 495 , 3 , 0 , 284 , 495 , 3 , 0 , 285 , 495 , 3 , 0 , 286 , 495 , 3 , 0 , 287 , 495 , 3 , 0 , 288 , 495 , 3 , 0 , 289 , 495 , 3 , 0 , 290 , 495 , 3 , 0 , 291 , 495 , 3 , 0 , 292 , 495 , 3 , 0 , 293 , 495 , 3 , 0 , 294 , 495 , 3 , 0 , 295 , 495 , 3 , 0 , 296 , 495 , 3 , 0 , 297 , 495 , 3 , 0 , 298 , 495 , 3 , 0 , 299 , 495 , 3 , 0 , 300 , 495 , 3 , 0 , 301 , 495 , 3 , 0 , 302 , 495 , 3 , 0 , 303 , 495 , 3 , 0 , 304 , 495 , 3 , 0 , 305 , 495 , 3 , 0 , 306 , 495 , 3 , 0 , 307 , 495 , 3 , 0 , 308 , 495 , 3 , 0 , 309 , 495 , 3 , 0 , 310 , 495 , 3 , 0 , 311 , 495 , 3 , 0 , 312 , 495 , 3 , 0 , 313 , 495 , 3 , 0 , 314 , 495 , 3 , 0 , 315 , 495 , 3 , 0 , 316 , 495 , 3 , 0 , 317 , 495 , 3 , 0 , 318 , 495 , 3 , 0 , 319 , 495 , 3 , 0 , 320 , 495 , 3 , 0 , 321 , 495 , 3 , 0 , 322 , 495 , 3 , 0 , 323 , 495 , 3 , 0 , 324 , 495 , 3 , 0 , 325 , 495 , 3 , 0 , 326 , 495 , 3 , 0 , 327 , 495 , 3 , 0 , 328 , 495 , 3 , 0 , 329 , 495 , 3 , 0 , 330 , 495 , 3 , 0 , 331 , 495 , 3 , 0 , 332 , 495 , 3 , 0 , 333 , 495 , 3 , 0 , 334 , 495 , 3 , 0 , 335 , 495 , 3 , 0 , 336 , 495 , 3 , 0 , 337 , 495 , 3 , 0 , 338 , 495 , 3 , 0 , 339 , 495 , 3 , 0 , 340 , 495 , 3 , 0 , 341 , 495 , 3 , 0 , 342 , 495 , 3 , 0 , 343 , 495 , 3 , 0 , 344 , 495 , 3 , 0 , 345 , 495 , 3 , 0 , 346 , 495 , 3 , 0 , 347 , 495 , 3 , 0 , 348 , 495 , 3 , 0 , 349 , 495 , 3 , 0 , 350 , 495 , 3 , 0 , 351 , 495 , 3 , 0 , 352 , 495 , 3 , 0 , 353 , 495 , 3 , 0 , 354 , 495 , 3 , 0 , 355 , 495 , 3 , 0 , 356 , 495 , 3 , 0 , 357 , 495 , 3 , 0 , 358 , 495 , 3 , 0 , 359 , 495 , 3 , 0 , 360 , 495 , 3 , 0 , 361 , 495 , 3 , 0 , 362 , 495 , 3 , 0 , 363 , 495 , 3 , 0 , 364 , 495 , 3 , 0 , 365 , 495 , 3 , 0 , 366 , 495 , 3 , 0 , 367 , 495 , 3 , 0 , 368 , 495 , 3 , 0 , 369 , 495 , 3 , 0 , 370 , 495 , 3 , 0 , 371 , 495 , 3 , 0 , 372 , 495 , 3 , 0 , 373 , 495 , 3 , 0 , 374 , 495 , 3 , 0 , 375 , 495 , 3 , 0 , 376 , 495 , 3 , 0 , 377 , 495 , 3 , 0 , 378 , 495 , 3 , 0 , 379 , 495 , 3 , 0 , 380 , 495 , 3 , 0 , 381 , 495 , 3 , 0 , 382 , 495 , 3 , 0 , 383 , 495 , 3 , 0 , 384 , 495 , 3 , 0 , 385 , 495 , 3 , 0 , 386 , 495 , 3 , 0 , 387 , 495 , 3 , 0 , 388 , 495 , 3 , 0 , 389 , 495 , 3 , 0 , 390 , 495 , 3 , 0 , 391 , 495 , 3 , 0 , 392 , 495 , 3 , 0 , 393 , 495 , 3 , 0 , 394 , 495 , 3 , 0 , 395 , 495 , 3 , 0 , 396 , 495 , 3 , 0 , 397 , 495 , 3 , 0 , 398 , 495 , 3 , 0 , 399 , 495 , 3 , 0 , 400 , 495 , 3 , 0 , 401 , 495 , 3 , 0 , 402 , 495 , 3 , 0 , 403 , 495 , 3 , 0 , 404 , 495 , 3 , 0 , 405 , 495 , 3 , 0 , 406 , 495 , 3 , 0 , 407 , 495 , 3 , 0 , 408 , 495 , 3 , 0 , 409 , 495 , 3 , 0 , 410 , 495 , 3 , 0 , 411 , 495 , 3 , 0 , 412 , 495 , 3 , 0 , 413 , 495 , 3 , 0 , 414 , 495 , 3 , 0 , 415 , 495 , 3 , 0 , 416 , 495 , 3 , 0 , 417 , 495 , 3 , 0 , 418 , 495 , 3 , 0 , 419 , 495 , 3 , 0 , 420 , 495 , 3 , 0 , 421 , 495 , 3 , 0 , 422 , 495 , 3 , 0 , 423 , 495 , 3 , 0 , 424 , 495 , 3 , 0 , 425 , 495 , 3 , 0 , 426 , 495 , 3 , 0 , 427 , 495 , 3 , 0 , 428 , 495 , 3 , 0 , 429 , 495 , 3 , 0 , 430 , 495 , 3 , 0 , 431 , 495 , 3 , 0 , 432 , 495 , 3 , 0 , 433 , 495 , 3 , 0 , 434 , 495 , 3 , 0 , 435 , 495 , 3 , 0 , 436 , 495 , 3 , 0 , 437 , 495 , 3 , 0 , 438 , 495 , 3 , 0 , 439 , 495 , 3 , 0 , 440 , 495 , 3 , 0 , 441 , 495 , 3 , 0 , 442 , 495 , 3 , 0 , 443 , 495 , 3 , 0 , 444 , 495 , 3 , 0 , 445 , 495 , 3 , 0 , 446 , 495 , 3 , 0 , 447 , 495 , 3 , 0 , 448 , 495 , 3 , 0 , 449 , 495 , 3 , 0 , 450 , 495 , 3 , 0 , 451 , 495 , 3 , 0 , 452 , 495 , 3 , 0 , 453 , 495 , 3 , 0 , 454 , 495 , 3 , 0 , 455 , 495 , 3 , 0 , 456 , 495 , 3 , 0 , 457 , 495 , 3 , 0 , 458 , 495 , 3 , 0 , 459 , 495 , 3 , 0 , 460 , 495 , 3 , 0 , 461 , 495 , 3 , 0 , 462 , 495 , 3 , 0 , 463 , 495 , 3 , 0 , 464 , 495 , 3 , 0 , 465 , 495 , 3 , 0 , 466 , 495 , 3 , 0 , 467 , 495 , 3 , 0 , 468 , 495 , 3 , 0 , 469 , 495 , 3 , 0 , 470 , 495 , 3 , 0 , 471 , 495 , 3 , 0 , 472 , 495 , 3 , 0 , 473 , 495 , 3 , 0 , 474 , 495 , 3 , 0 , 475 , 495 , 3 , 0 , 476 , 495 , 3 , 0 , 477 , 495 , 3 , 0 , 478 , 495 , 3 , 0 , 479 , 495 , 3 , 0 , 480 , 495 , 3 , 0 , 481 , 495 , 3 , 0 , 482 , 495 , 3 , 0 , 483 , 495 , 3 , 0 , 484 , 495 , 3 , 0 , 485 , 495 , 3 , 0 , 486 , 495 , 3 , 0 , 487 , 495 , 3 , 0 , 488 , 495 , 3 , 0 , 489 , 495 , 3 , 0 , 490 , 495 , 3 , 0 , 491 , 495 , 3 , 0 , 492 , 495 , 3 , 0 , 493 , 495 , 3 , 0 , 494 , 495 };
  
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
void Basis_F2_TET_I10_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 9 (I9) has 715 basis functions on a tetrahedron that are 2-forms in 3D
  int    numFields = 715;
  EField fieldType = FIELD_FORM_2;
  int    spaceDim  = 3;

  // temporaries
  int countPt  = 0;               // point counter
  Teuchos::Array<int> indexV(3);  // multi-index for values
  Teuchos::Array<int> indexD(2);  // multi-index for divergences

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
                        ">>> ERROR (Basis_F2_TET_I10_FEM_FIAT): Evaluation point is outside the TET reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(286,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(715,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,10,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<715;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<715;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm2_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<715;i++) {
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
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TET_I10_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TET_I10_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_DIV:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(286,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(286,286);
      Teuchos::SerialDenseMatrix<int,Scalar> div(715,numPoints);
      div.putScalar(0.0);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,10,inputPoints,expansions);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm2_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<715;i++) {
          outputValues(countPt,i) = div(i,countPt);
        }
      }
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
                          ">>> ERROR (Basis_F2_TET_I10_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F2_TET_I10_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F2_TET_I10_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F2_TET_I10_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F2_TET_I10_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F2_TET_I10_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F2_TET_I10_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F2_TET_I10_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TET;
}



template<class Scalar>
EBasis Basis_F2_TET_I10_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F2_TET_I10_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F2_TET_I10_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

