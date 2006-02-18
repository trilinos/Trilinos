/*@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#ifndef THYRA_IFPACK_TYPES_HPP
#define THYRA_IFPACK_TYPES_HPP

#include "Ifpack_ConfigDefs.h"

namespace Thyra {

namespace Ifpack {

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
enum EPrecType {
  POINT_RELAXATION
  ,POINT_RELAXATION_STAND_ALONE
  ,BLOCK_RELAXATION
  ,BLOCK_RELAXATION_STAND_ALONE
#ifdef HAVE_IFPACK_AMESOS
  ,BLOCK_RELAXATION_STAND_ALONE_ILU
  ,BLOCK_RELAXATION_STAND_ALONE_AMESOS
  ,BLOCK_RELAXATION_AMESOS
  ,AMESOS
  ,AMESOS_STAND_ALONE
#endif // HAVE_IFPACK_AMESOS
  ,IC
  ,IC_STAND_ALONE
  ,ICT
  ,ICT_STAND_ALONE
  ,ILU
  ,ILU_STAND_ALONE
  ,ILUT
  ,ILUT_STAND_ALONE
#ifdef HAVE_IFPACK_SPARSKIT
  ,SPARSKIT
#endif // HAVE_IFPACK_SPARSKIT
}; 

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
const int numPrecTypes =
+4
#ifdef HAVE_IFPACK_AMESOS
+5
#endif
#ifdef HAVE_IFPACK_SPARSKIT
+1
#endif
+8
;

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
extern const EPrecType precTypeValues[numPrecTypes];

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
extern const char* precTypeNames[numPrecTypes];

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
extern const bool supportsUnsymmetric[numPrecTypes];

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
inline const char* toString(const EPrecType precType)
{ return precTypeNames[precType]; }

/** \brief .
\ingroup Ifpack_Thyra_adapters_grp
*/
EPrecType precTypeNameToEnum(const std::string& precTypeName, const std::string& paramName);

} // namespace Ifpack

} // namespace Thyra

#endif // THYRA_IFPACK_TYPES_HPP
