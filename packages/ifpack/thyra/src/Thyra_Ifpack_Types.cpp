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

#include "Thyra_Ifpack_Types.hpp"
#include "Teuchos_StringToIntMap.hpp"

namespace {

const Teuchos::StringToIntMap
precTypeNameToIntMap(
  "parameter \"Prec Type\"", Thyra::Ifpack::numPrecTypes, Thyra::Ifpack::precTypeNames
  );

} // namespace

namespace Thyra {

const Ifpack::EPrecType Ifpack::precTypeValues[Ifpack::numPrecTypes] =
{
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

const char* Ifpack::precTypeNames[Ifpack::numPrecTypes] =
{
  "point relaxation"
  ,"point relaxation stand-alone"
  ,"block relaxation"
  ,"block relaxation stand-alone"
#ifdef HAVE_IFPACK_AMESOS
  ,"block relaxation stand-alone (ILU)"
  ,"block relaxation stand-alone (Amesos)"
  ,"block relaxation (Amesos)"
  ,"Amesos"
  ,"LU"
#endif
  ,"IC"
  ,"IC stand-alone"
  ,"ICT"
  ,"ICT stand-alone"
  ,"ILU"
  ,"ILU stand-alone"
  ,"ILUT"
  ,"ILUT stand-alone"
#ifdef HAVE_IFPACK_SPARSKIT
  ,"SPARSKIT"
#endif
};

const bool Ifpack::supportsUnsymmetric[Ifpack::numPrecTypes] =
{
  true // point relaxation
  ,true // point relaxation stand-alone
  ,true // block relaxation
  ,true // block relaxation stand-alone
#ifdef HAVE_IFPACK_AMESOS
  ,true // block relaxation stand-alone (ILU)
  ,true // block relaxation stand-alone (Amesos)
  ,true // block relaxation (Amesos)
  ,true // Amesos
  ,true // LU
#endif
  ,false // IC
  ,false // IC stand-alone
  ,false // ICT
  ,false // ICT stand-alone
  ,true // ILU
  ,true // ILU stand-alone
  ,true // ILUT
  ,true // ILUT stand-alone
#ifdef HAVE_IFPACK_SPARSKIT
  ,true // SPARSKIT
#endif
};

Ifpack::EPrecType Ifpack::precTypeNameToEnum(
  const std::string& precTypeName
  ,const std::string& paramName
  )
{
  return static_cast<EPrecType>(precTypeNameToIntMap(precTypeName));
}

} // namespace Thyra
