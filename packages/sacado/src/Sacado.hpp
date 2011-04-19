// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_HPP
#define SACADO_HPP

// Version string
#include "Sacado_Version.hpp"

// Declarations of all overloaded math functions
#include "Sacado_MathFunctions.hpp"

// Traits for all of the Sacado classes -- Include these first so they are all
// defined before any nesting of AD classes
#include "Sacado_Fad_ExpressionTraits.hpp"
#include "Sacado_Fad_DFadTraits.hpp"
#include "Sacado_Fad_SFadTraits.hpp"
#include "Sacado_Fad_SLFadTraits.hpp"
#include "Sacado_Fad_DMFadTraits.hpp"
#include "Sacado_Fad_DVFadTraits.hpp"
#include "Sacado_ELRFad_ExpressionTraits.hpp"
#include "Sacado_ELRFad_DFadTraits.hpp"
#include "Sacado_ELRFad_SFadTraits.hpp"
#include "Sacado_ELRFad_SLFadTraits.hpp"
// #include "Sacado_CacheFad_ExpressionTraits.hpp"
// #include "Sacado_CacheFad_DFadTraits.hpp"
// #include "Sacado_CacheFad_SFadTraits.hpp"
// #include "Sacado_ELRCacheFad_SLFadTraits.hpp"
#include "Sacado_ELRCacheFad_ExpressionTraits.hpp"
#include "Sacado_ELRCacheFad_DFadTraits.hpp"
#include "Sacado_ELRCacheFad_SFadTraits.hpp"
#include "Sacado_ELRCacheFad_SLFadTraits.hpp"
#include "Sacado_LFad_LogicalSparseTraits.hpp"
#include "Sacado_ScalarFlopCounterTraits.hpp"
#include "Sacado_Tay_TaylorTraits.hpp"
#include "Sacado_trad_Traits.hpp"
#include "Sacado_trad2_Traits.hpp"
#include "Sacado_tradvec_Traits.hpp"

// Standard forward AD classes
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_Fad_SFad.hpp"
#include "Sacado_Fad_SLFad.hpp"
#include "Sacado_Fad_MemPoolManager.hpp"
#include "Sacado_Fad_DMFad.hpp"
#include "Sacado_LFad_LogicalSparse.hpp"
#include "Sacado_Fad_DVFad.hpp"
#include "Sacado_Fad_Vector.hpp"

// Expression-level-reverse forward AD classes
#include "Sacado_ELRFad_DFad.hpp"
#include "Sacado_ELRFad_SFad.hpp"
#include "Sacado_ELRFad_SLFad.hpp"

// Caching forward AD classes
// Not including CacheFad by default since AIX has issues with it.
// This class is not production anyway.
//#include "Sacado_CacheFad_DFad.hpp"
//#include "Sacado_CacheFad_SFad.hpp"
//#include "Sacado_CacheFad_SLFad.hpp"

// Caching expression-level reverse mode forward AD classes
#include "Sacado_ELRCacheFad_DFad.hpp"
#include "Sacado_ELRCacheFad_SFad.hpp"
#include "Sacado_ELRCacheFad_SLFad.hpp"

// Reverse AD classes
#include "Sacado_trad.hpp"
#include "Sacado_trad2.hpp"
#include "Sacado_tradvec.hpp"

// Taylor polynomial AD classes
#include "Sacado_Tay_Taylor.hpp"

// Flop-counting classes
#include "Sacado_ScalarFlopCounter.hpp"

#endif // SACADO_HPP 
