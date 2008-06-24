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

/** \file   Intrepid_F2_TET_I9_FEM_FIATDef.hpp
    \brief  Definition file for FEM basis functions of degree 8 for 2-forms on TET cells.
    \author Created by R. Kirby via the FIAT project
*/

namespace Intrepid {

template<class Scalar>
void Basis_F2_TET_I9_FEM_FIAT<Scalar>::initialize() {

  // Basis-dependent initializations
  int tagSize  = 4;         // size of DoF tag
  int posScDim = 0;         // position in the tag, counting from 0, of the subcell dim 
  int posScId  = 1;         // position in the tag, counting from 0, of the subcell id
  int posBfId  = 2;         // position in the tag, counting from 0, of DoF Id relative to the subcell

  // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration 
  int tags[] = { 2 , 0 , 0 , 45 , 2 , 0 , 1 , 45 , 2 , 0 , 2 , 45 , 2 , 0 , 3 , 45 , 2 , 0 , 4 , 45 , 2 , 0 , 5 , 45 , 2 , 0 , 6 , 45 , 2 , 0 , 7 , 45 , 2 , 0 , 8 , 45 , 2 , 0 , 9 , 45 , 2 , 0 , 10 , 45 , 2 , 0 , 11 , 45 , 2 , 0 , 12 , 45 , 2 , 0 , 13 , 45 , 2 , 0 , 14 , 45 , 2 , 0 , 15 , 45 , 2 , 0 , 16 , 45 , 2 , 0 , 17 , 45 , 2 , 0 , 18 , 45 , 2 , 0 , 19 , 45 , 2 , 0 , 20 , 45 , 2 , 0 , 21 , 45 , 2 , 0 , 22 , 45 , 2 , 0 , 23 , 45 , 2 , 0 , 24 , 45 , 2 , 0 , 25 , 45 , 2 , 0 , 26 , 45 , 2 , 0 , 27 , 45 , 2 , 0 , 28 , 45 , 2 , 0 , 29 , 45 , 2 , 0 , 30 , 45 , 2 , 0 , 31 , 45 , 2 , 0 , 32 , 45 , 2 , 0 , 33 , 45 , 2 , 0 , 34 , 45 , 2 , 0 , 35 , 45 , 2 , 0 , 36 , 45 , 2 , 0 , 37 , 45 , 2 , 0 , 38 , 45 , 2 , 0 , 39 , 45 , 2 , 0 , 40 , 45 , 2 , 0 , 41 , 45 , 2 , 0 , 42 , 45 , 2 , 0 , 43 , 45 , 2 , 0 , 44 , 45 , 2 , 1 , 0 , 45 , 2 , 1 , 1 , 45 , 2 , 1 , 2 , 45 , 2 , 1 , 3 , 45 , 2 , 1 , 4 , 45 , 2 , 1 , 5 , 45 , 2 , 1 , 6 , 45 , 2 , 1 , 7 , 45 , 2 , 1 , 8 , 45 , 2 , 1 , 9 , 45 , 2 , 1 , 10 , 45 , 2 , 1 , 11 , 45 , 2 , 1 , 12 , 45 , 2 , 1 , 13 , 45 , 2 , 1 , 14 , 45 , 2 , 1 , 15 , 45 , 2 , 1 , 16 , 45 , 2 , 1 , 17 , 45 , 2 , 1 , 18 , 45 , 2 , 1 , 19 , 45 , 2 , 1 , 20 , 45 , 2 , 1 , 21 , 45 , 2 , 1 , 22 , 45 , 2 , 1 , 23 , 45 , 2 , 1 , 24 , 45 , 2 , 1 , 25 , 45 , 2 , 1 , 26 , 45 , 2 , 1 , 27 , 45 , 2 , 1 , 28 , 45 , 2 , 1 , 29 , 45 , 2 , 1 , 30 , 45 , 2 , 1 , 31 , 45 , 2 , 1 , 32 , 45 , 2 , 1 , 33 , 45 , 2 , 1 , 34 , 45 , 2 , 1 , 35 , 45 , 2 , 1 , 36 , 45 , 2 , 1 , 37 , 45 , 2 , 1 , 38 , 45 , 2 , 1 , 39 , 45 , 2 , 1 , 40 , 45 , 2 , 1 , 41 , 45 , 2 , 1 , 42 , 45 , 2 , 1 , 43 , 45 , 2 , 1 , 44 , 45 , 2 , 2 , 0 , 45 , 2 , 2 , 1 , 45 , 2 , 2 , 2 , 45 , 2 , 2 , 3 , 45 , 2 , 2 , 4 , 45 , 2 , 2 , 5 , 45 , 2 , 2 , 6 , 45 , 2 , 2 , 7 , 45 , 2 , 2 , 8 , 45 , 2 , 2 , 9 , 45 , 2 , 2 , 10 , 45 , 2 , 2 , 11 , 45 , 2 , 2 , 12 , 45 , 2 , 2 , 13 , 45 , 2 , 2 , 14 , 45 , 2 , 2 , 15 , 45 , 2 , 2 , 16 , 45 , 2 , 2 , 17 , 45 , 2 , 2 , 18 , 45 , 2 , 2 , 19 , 45 , 2 , 2 , 20 , 45 , 2 , 2 , 21 , 45 , 2 , 2 , 22 , 45 , 2 , 2 , 23 , 45 , 2 , 2 , 24 , 45 , 2 , 2 , 25 , 45 , 2 , 2 , 26 , 45 , 2 , 2 , 27 , 45 , 2 , 2 , 28 , 45 , 2 , 2 , 29 , 45 , 2 , 2 , 30 , 45 , 2 , 2 , 31 , 45 , 2 , 2 , 32 , 45 , 2 , 2 , 33 , 45 , 2 , 2 , 34 , 45 , 2 , 2 , 35 , 45 , 2 , 2 , 36 , 45 , 2 , 2 , 37 , 45 , 2 , 2 , 38 , 45 , 2 , 2 , 39 , 45 , 2 , 2 , 40 , 45 , 2 , 2 , 41 , 45 , 2 , 2 , 42 , 45 , 2 , 2 , 43 , 45 , 2 , 2 , 44 , 45 , 2 , 3 , 0 , 45 , 2 , 3 , 1 , 45 , 2 , 3 , 2 , 45 , 2 , 3 , 3 , 45 , 2 , 3 , 4 , 45 , 2 , 3 , 5 , 45 , 2 , 3 , 6 , 45 , 2 , 3 , 7 , 45 , 2 , 3 , 8 , 45 , 2 , 3 , 9 , 45 , 2 , 3 , 10 , 45 , 2 , 3 , 11 , 45 , 2 , 3 , 12 , 45 , 2 , 3 , 13 , 45 , 2 , 3 , 14 , 45 , 2 , 3 , 15 , 45 , 2 , 3 , 16 , 45 , 2 , 3 , 17 , 45 , 2 , 3 , 18 , 45 , 2 , 3 , 19 , 45 , 2 , 3 , 20 , 45 , 2 , 3 , 21 , 45 , 2 , 3 , 22 , 45 , 2 , 3 , 23 , 45 , 2 , 3 , 24 , 45 , 2 , 3 , 25 , 45 , 2 , 3 , 26 , 45 , 2 , 3 , 27 , 45 , 2 , 3 , 28 , 45 , 2 , 3 , 29 , 45 , 2 , 3 , 30 , 45 , 2 , 3 , 31 , 45 , 2 , 3 , 32 , 45 , 2 , 3 , 33 , 45 , 2 , 3 , 34 , 45 , 2 , 3 , 35 , 45 , 2 , 3 , 36 , 45 , 2 , 3 , 37 , 45 , 2 , 3 , 38 , 45 , 2 , 3 , 39 , 45 , 2 , 3 , 40 , 45 , 2 , 3 , 41 , 45 , 2 , 3 , 42 , 45 , 2 , 3 , 43 , 45 , 2 , 3 , 44 , 45 , 3 , 0 , 0 , 360 , 3 , 0 , 1 , 360 , 3 , 0 , 2 , 360 , 3 , 0 , 3 , 360 , 3 , 0 , 4 , 360 , 3 , 0 , 5 , 360 , 3 , 0 , 6 , 360 , 3 , 0 , 7 , 360 , 3 , 0 , 8 , 360 , 3 , 0 , 9 , 360 , 3 , 0 , 10 , 360 , 3 , 0 , 11 , 360 , 3 , 0 , 12 , 360 , 3 , 0 , 13 , 360 , 3 , 0 , 14 , 360 , 3 , 0 , 15 , 360 , 3 , 0 , 16 , 360 , 3 , 0 , 17 , 360 , 3 , 0 , 18 , 360 , 3 , 0 , 19 , 360 , 3 , 0 , 20 , 360 , 3 , 0 , 21 , 360 , 3 , 0 , 22 , 360 , 3 , 0 , 23 , 360 , 3 , 0 , 24 , 360 , 3 , 0 , 25 , 360 , 3 , 0 , 26 , 360 , 3 , 0 , 27 , 360 , 3 , 0 , 28 , 360 , 3 , 0 , 29 , 360 , 3 , 0 , 30 , 360 , 3 , 0 , 31 , 360 , 3 , 0 , 32 , 360 , 3 , 0 , 33 , 360 , 3 , 0 , 34 , 360 , 3 , 0 , 35 , 360 , 3 , 0 , 36 , 360 , 3 , 0 , 37 , 360 , 3 , 0 , 38 , 360 , 3 , 0 , 39 , 360 , 3 , 0 , 40 , 360 , 3 , 0 , 41 , 360 , 3 , 0 , 42 , 360 , 3 , 0 , 43 , 360 , 3 , 0 , 44 , 360 , 3 , 0 , 45 , 360 , 3 , 0 , 46 , 360 , 3 , 0 , 47 , 360 , 3 , 0 , 48 , 360 , 3 , 0 , 49 , 360 , 3 , 0 , 50 , 360 , 3 , 0 , 51 , 360 , 3 , 0 , 52 , 360 , 3 , 0 , 53 , 360 , 3 , 0 , 54 , 360 , 3 , 0 , 55 , 360 , 3 , 0 , 56 , 360 , 3 , 0 , 57 , 360 , 3 , 0 , 58 , 360 , 3 , 0 , 59 , 360 , 3 , 0 , 60 , 360 , 3 , 0 , 61 , 360 , 3 , 0 , 62 , 360 , 3 , 0 , 63 , 360 , 3 , 0 , 64 , 360 , 3 , 0 , 65 , 360 , 3 , 0 , 66 , 360 , 3 , 0 , 67 , 360 , 3 , 0 , 68 , 360 , 3 , 0 , 69 , 360 , 3 , 0 , 70 , 360 , 3 , 0 , 71 , 360 , 3 , 0 , 72 , 360 , 3 , 0 , 73 , 360 , 3 , 0 , 74 , 360 , 3 , 0 , 75 , 360 , 3 , 0 , 76 , 360 , 3 , 0 , 77 , 360 , 3 , 0 , 78 , 360 , 3 , 0 , 79 , 360 , 3 , 0 , 80 , 360 , 3 , 0 , 81 , 360 , 3 , 0 , 82 , 360 , 3 , 0 , 83 , 360 , 3 , 0 , 84 , 360 , 3 , 0 , 85 , 360 , 3 , 0 , 86 , 360 , 3 , 0 , 87 , 360 , 3 , 0 , 88 , 360 , 3 , 0 , 89 , 360 , 3 , 0 , 90 , 360 , 3 , 0 , 91 , 360 , 3 , 0 , 92 , 360 , 3 , 0 , 93 , 360 , 3 , 0 , 94 , 360 , 3 , 0 , 95 , 360 , 3 , 0 , 96 , 360 , 3 , 0 , 97 , 360 , 3 , 0 , 98 , 360 , 3 , 0 , 99 , 360 , 3 , 0 , 100 , 360 , 3 , 0 , 101 , 360 , 3 , 0 , 102 , 360 , 3 , 0 , 103 , 360 , 3 , 0 , 104 , 360 , 3 , 0 , 105 , 360 , 3 , 0 , 106 , 360 , 3 , 0 , 107 , 360 , 3 , 0 , 108 , 360 , 3 , 0 , 109 , 360 , 3 , 0 , 110 , 360 , 3 , 0 , 111 , 360 , 3 , 0 , 112 , 360 , 3 , 0 , 113 , 360 , 3 , 0 , 114 , 360 , 3 , 0 , 115 , 360 , 3 , 0 , 116 , 360 , 3 , 0 , 117 , 360 , 3 , 0 , 118 , 360 , 3 , 0 , 119 , 360 , 3 , 0 , 120 , 360 , 3 , 0 , 121 , 360 , 3 , 0 , 122 , 360 , 3 , 0 , 123 , 360 , 3 , 0 , 124 , 360 , 3 , 0 , 125 , 360 , 3 , 0 , 126 , 360 , 3 , 0 , 127 , 360 , 3 , 0 , 128 , 360 , 3 , 0 , 129 , 360 , 3 , 0 , 130 , 360 , 3 , 0 , 131 , 360 , 3 , 0 , 132 , 360 , 3 , 0 , 133 , 360 , 3 , 0 , 134 , 360 , 3 , 0 , 135 , 360 , 3 , 0 , 136 , 360 , 3 , 0 , 137 , 360 , 3 , 0 , 138 , 360 , 3 , 0 , 139 , 360 , 3 , 0 , 140 , 360 , 3 , 0 , 141 , 360 , 3 , 0 , 142 , 360 , 3 , 0 , 143 , 360 , 3 , 0 , 144 , 360 , 3 , 0 , 145 , 360 , 3 , 0 , 146 , 360 , 3 , 0 , 147 , 360 , 3 , 0 , 148 , 360 , 3 , 0 , 149 , 360 , 3 , 0 , 150 , 360 , 3 , 0 , 151 , 360 , 3 , 0 , 152 , 360 , 3 , 0 , 153 , 360 , 3 , 0 , 154 , 360 , 3 , 0 , 155 , 360 , 3 , 0 , 156 , 360 , 3 , 0 , 157 , 360 , 3 , 0 , 158 , 360 , 3 , 0 , 159 , 360 , 3 , 0 , 160 , 360 , 3 , 0 , 161 , 360 , 3 , 0 , 162 , 360 , 3 , 0 , 163 , 360 , 3 , 0 , 164 , 360 , 3 , 0 , 165 , 360 , 3 , 0 , 166 , 360 , 3 , 0 , 167 , 360 , 3 , 0 , 168 , 360 , 3 , 0 , 169 , 360 , 3 , 0 , 170 , 360 , 3 , 0 , 171 , 360 , 3 , 0 , 172 , 360 , 3 , 0 , 173 , 360 , 3 , 0 , 174 , 360 , 3 , 0 , 175 , 360 , 3 , 0 , 176 , 360 , 3 , 0 , 177 , 360 , 3 , 0 , 178 , 360 , 3 , 0 , 179 , 360 , 3 , 0 , 180 , 360 , 3 , 0 , 181 , 360 , 3 , 0 , 182 , 360 , 3 , 0 , 183 , 360 , 3 , 0 , 184 , 360 , 3 , 0 , 185 , 360 , 3 , 0 , 186 , 360 , 3 , 0 , 187 , 360 , 3 , 0 , 188 , 360 , 3 , 0 , 189 , 360 , 3 , 0 , 190 , 360 , 3 , 0 , 191 , 360 , 3 , 0 , 192 , 360 , 3 , 0 , 193 , 360 , 3 , 0 , 194 , 360 , 3 , 0 , 195 , 360 , 3 , 0 , 196 , 360 , 3 , 0 , 197 , 360 , 3 , 0 , 198 , 360 , 3 , 0 , 199 , 360 , 3 , 0 , 200 , 360 , 3 , 0 , 201 , 360 , 3 , 0 , 202 , 360 , 3 , 0 , 203 , 360 , 3 , 0 , 204 , 360 , 3 , 0 , 205 , 360 , 3 , 0 , 206 , 360 , 3 , 0 , 207 , 360 , 3 , 0 , 208 , 360 , 3 , 0 , 209 , 360 , 3 , 0 , 210 , 360 , 3 , 0 , 211 , 360 , 3 , 0 , 212 , 360 , 3 , 0 , 213 , 360 , 3 , 0 , 214 , 360 , 3 , 0 , 215 , 360 , 3 , 0 , 216 , 360 , 3 , 0 , 217 , 360 , 3 , 0 , 218 , 360 , 3 , 0 , 219 , 360 , 3 , 0 , 220 , 360 , 3 , 0 , 221 , 360 , 3 , 0 , 222 , 360 , 3 , 0 , 223 , 360 , 3 , 0 , 224 , 360 , 3 , 0 , 225 , 360 , 3 , 0 , 226 , 360 , 3 , 0 , 227 , 360 , 3 , 0 , 228 , 360 , 3 , 0 , 229 , 360 , 3 , 0 , 230 , 360 , 3 , 0 , 231 , 360 , 3 , 0 , 232 , 360 , 3 , 0 , 233 , 360 , 3 , 0 , 234 , 360 , 3 , 0 , 235 , 360 , 3 , 0 , 236 , 360 , 3 , 0 , 237 , 360 , 3 , 0 , 238 , 360 , 3 , 0 , 239 , 360 , 3 , 0 , 240 , 360 , 3 , 0 , 241 , 360 , 3 , 0 , 242 , 360 , 3 , 0 , 243 , 360 , 3 , 0 , 244 , 360 , 3 , 0 , 245 , 360 , 3 , 0 , 246 , 360 , 3 , 0 , 247 , 360 , 3 , 0 , 248 , 360 , 3 , 0 , 249 , 360 , 3 , 0 , 250 , 360 , 3 , 0 , 251 , 360 , 3 , 0 , 252 , 360 , 3 , 0 , 253 , 360 , 3 , 0 , 254 , 360 , 3 , 0 , 255 , 360 , 3 , 0 , 256 , 360 , 3 , 0 , 257 , 360 , 3 , 0 , 258 , 360 , 3 , 0 , 259 , 360 , 3 , 0 , 260 , 360 , 3 , 0 , 261 , 360 , 3 , 0 , 262 , 360 , 3 , 0 , 263 , 360 , 3 , 0 , 264 , 360 , 3 , 0 , 265 , 360 , 3 , 0 , 266 , 360 , 3 , 0 , 267 , 360 , 3 , 0 , 268 , 360 , 3 , 0 , 269 , 360 , 3 , 0 , 270 , 360 , 3 , 0 , 271 , 360 , 3 , 0 , 272 , 360 , 3 , 0 , 273 , 360 , 3 , 0 , 274 , 360 , 3 , 0 , 275 , 360 , 3 , 0 , 276 , 360 , 3 , 0 , 277 , 360 , 3 , 0 , 278 , 360 , 3 , 0 , 279 , 360 , 3 , 0 , 280 , 360 , 3 , 0 , 281 , 360 , 3 , 0 , 282 , 360 , 3 , 0 , 283 , 360 , 3 , 0 , 284 , 360 , 3 , 0 , 285 , 360 , 3 , 0 , 286 , 360 , 3 , 0 , 287 , 360 , 3 , 0 , 288 , 360 , 3 , 0 , 289 , 360 , 3 , 0 , 290 , 360 , 3 , 0 , 291 , 360 , 3 , 0 , 292 , 360 , 3 , 0 , 293 , 360 , 3 , 0 , 294 , 360 , 3 , 0 , 295 , 360 , 3 , 0 , 296 , 360 , 3 , 0 , 297 , 360 , 3 , 0 , 298 , 360 , 3 , 0 , 299 , 360 , 3 , 0 , 300 , 360 , 3 , 0 , 301 , 360 , 3 , 0 , 302 , 360 , 3 , 0 , 303 , 360 , 3 , 0 , 304 , 360 , 3 , 0 , 305 , 360 , 3 , 0 , 306 , 360 , 3 , 0 , 307 , 360 , 3 , 0 , 308 , 360 , 3 , 0 , 309 , 360 , 3 , 0 , 310 , 360 , 3 , 0 , 311 , 360 , 3 , 0 , 312 , 360 , 3 , 0 , 313 , 360 , 3 , 0 , 314 , 360 , 3 , 0 , 315 , 360 , 3 , 0 , 316 , 360 , 3 , 0 , 317 , 360 , 3 , 0 , 318 , 360 , 3 , 0 , 319 , 360 , 3 , 0 , 320 , 360 , 3 , 0 , 321 , 360 , 3 , 0 , 322 , 360 , 3 , 0 , 323 , 360 , 3 , 0 , 324 , 360 , 3 , 0 , 325 , 360 , 3 , 0 , 326 , 360 , 3 , 0 , 327 , 360 , 3 , 0 , 328 , 360 , 3 , 0 , 329 , 360 , 3 , 0 , 330 , 360 , 3 , 0 , 331 , 360 , 3 , 0 , 332 , 360 , 3 , 0 , 333 , 360 , 3 , 0 , 334 , 360 , 3 , 0 , 335 , 360 , 3 , 0 , 336 , 360 , 3 , 0 , 337 , 360 , 3 , 0 , 338 , 360 , 3 , 0 , 339 , 360 , 3 , 0 , 340 , 360 , 3 , 0 , 341 , 360 , 3 , 0 , 342 , 360 , 3 , 0 , 343 , 360 , 3 , 0 , 344 , 360 , 3 , 0 , 345 , 360 , 3 , 0 , 346 , 360 , 3 , 0 , 347 , 360 , 3 , 0 , 348 , 360 , 3 , 0 , 349 , 360 , 3 , 0 , 350 , 360 , 3 , 0 , 351 , 360 , 3 , 0 , 352 , 360 , 3 , 0 , 353 , 360 , 3 , 0 , 354 , 360 , 3 , 0 , 355 , 360 , 3 , 0 , 356 , 360 , 3 , 0 , 357 , 360 , 3 , 0 , 358 , 360 , 3 , 0 , 359 , 360 };
  
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
void Basis_F2_TET_I9_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
          const Teuchos::Array< Point<Scalar> >& inputPoints,
          const EOperator                        operatorType) const {

  // Determine parameters to shape outputValues: number of points = size of input array
  int numPoints = inputPoints.size();
	
  // Incomplete polynomial basis of degree 8 (I8) has 540 basis functions on a tetrahedron that are 2-forms in 3D
  int    numFields = 540;
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
                        ">>> ERROR (Basis_F2_TET_I9_FEM_FIAT): Evaluation point is outside the TET reference cell");
  }
#endif
  switch(operatorType) {
    case OPERATOR_VALUE:   
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(220,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> result(540,numPoints);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,9,inputPoints,expansions);
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm0_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<540;i++) {
          outputValues(countPt,i,0) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm1_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<540;i++) {
          outputValues(countPt,i,1) = result(i,countPt);
        }
      }
      result.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1,*vdm2_,expansions,0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<540;i++) {
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
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TET_I9_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_CURL:
    {
      TEST_FOR_EXCEPTION( true , std::invalid_argument, ">>>ERROR F2_TET_I9_FEM_FIAT: operator not implemented" );
    }
    break; 
    case OPERATOR_DIV:
    {
      Teuchos::SerialDenseMatrix<int,Scalar> expansions(220,numPoints);
      Teuchos::SerialDenseMatrix<int,Scalar> tmp(220,220);
      Teuchos::SerialDenseMatrix<int,Scalar> div(540,numPoints);
      div.putScalar(0.0);
      OrthogonalExpansions<Scalar>::tabulate(CELL_TET,9,inputPoints,expansions);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats0_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm0_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats1_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm1_,tmp,1.0);
      tmp.multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*dmats2_,expansions,0.0);
      div.multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*vdm2_,tmp,1.0);
      for (countPt=0;countPt<numPoints;countPt++) {
        for (int i=0;i<540;i++) {
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
                          ">>> ERROR (Basis_F2_TET_I9_FEM_DEFAULT): Invalid operator type");
  }


}


template<class Scalar>
void Basis_F2_TET_I9_FEM_FIAT<Scalar>::getValues(FieldContainer<Scalar>&                  outputValues,
                                                    const Teuchos::Array< Point<Scalar> >& inputPoints,
                                                    const Cell<Scalar>&                    cell) const {
  TEST_FOR_EXCEPTION( (true),
                      std::logic_error,
                      ">>> ERROR (Basis_F2_TET_I9_FEM_FIAT): FEM Basis calling an FV/D member function");
}



template<class Scalar>
int Basis_F2_TET_I9_FEM_FIAT<Scalar>::getNumLocalDof() const {
    return numDof_;   
}



template<class Scalar>
int Basis_F2_TET_I9_FEM_FIAT<Scalar>::getLocalDofEnumeration(const LocalDofTag dofTag) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return tagToEnum_[dofTag.tag_[0]][dofTag.tag_[1]][dofTag.tag_[2]];
}



template<class Scalar>
LocalDofTag Basis_F2_TET_I9_FEM_FIAT<Scalar>::getLocalDofTag(int id) {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_[id];
}



template<class Scalar>
const Teuchos::Array<LocalDofTag> & Basis_F2_TET_I9_FEM_FIAT<Scalar>::getAllLocalDofTags() {
  if (!isSet_) {
    initialize();
    isSet_ = true;
  }
  return enumToTag_;
}



template<class Scalar>
ECell Basis_F2_TET_I9_FEM_FIAT<Scalar>::getCellType() const {
  return CELL_TET;
}



template<class Scalar>
EBasis Basis_F2_TET_I9_FEM_FIAT<Scalar>::getBasisType() const {
  return BASIS_FEM_FIAT;
}



template<class Scalar>
ECoordinates Basis_F2_TET_I9_FEM_FIAT<Scalar>::getCoordinateSystem() const {
  return COORDINATES_CARTESIAN;
}



template<class Scalar>
int Basis_F2_TET_I9_FEM_FIAT<Scalar>::getDegree() const {
  return 1;
}


}// namespace Intrepid

