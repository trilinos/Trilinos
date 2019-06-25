// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/base/CoordinateSystems.hpp>
#include <sstream>                      // for operator<<, basic_ostream, etc
#include <stk_util/util/ReportHandler.hpp>  // for ThrowErrorMsgIf
#include <stk_util/util/string_case_compare.hpp>  // for not_equal_case
#include "Shards_Array.hpp"             // for ArrayDimTag::size_type, etc


namespace stk {
namespace mesh {

//----------------------------------------------------------------------

SHARDS_ARRAY_DIM_TAG_SIMPLE_IMPLEMENTATION( SimpleArrayTag )

//----------------------------------------------------------------------

const Cartesian2d & Cartesian2d::tag()
{ static const Cartesian2d self ; return self ; }

const char * Cartesian2d::name() const
{ static const char n[] = "Cartesian2d" ; return n ; }

//----------------------------------------------------------------------

const Cartesian3d & Cartesian3d::tag()
{ static const Cartesian3d self ; return self ; }

const char * Cartesian3d::name() const
{ static const char n[] = "Cartesian3d" ; return n ; }

//----------------------------------------------------------------------

const Cylindrical & Cylindrical::tag()
{ static const Cylindrical self ; return self ; }

const char * Cylindrical::name() const
{ static const char n[] = "Cylindrical" ; return n ; }

//----------------------------------------------------------------------

const FullTensor & FullTensor::tag()
{ static const FullTensor self ; return self ; }

const char * FullTensor::name() const
{ static const char n[] = "FullTensor" ; return n ; }

//----------------------------------------------------------------------

const FullTensor22 & FullTensor22::tag()
{ static const FullTensor22 self ; return self ; }

const char * FullTensor22::name() const
{ static const char n[] = "FullTensor22" ; return n ; }

//----------------------------------------------------------------------

const SymmetricTensor33 & SymmetricTensor33::tag()
{ static const SymmetricTensor33 self ; return self ; }

const char * SymmetricTensor::name() const
{ static const char n[] = "SymmetricTensor" ; return n ; }

//----------------------------------------------------------------------

const SymmetricTensor31 & SymmetricTensor31::tag()
{ static const SymmetricTensor31 self ; return self ; }

const char * SymmetricTensor31::name() const
{ static const char n[] = "SymmetricTensor31" ; return n ; }

//----------------------------------------------------------------------

const SymmetricTensor21 & SymmetricTensor21::tag()
{ static const SymmetricTensor21 self ; return self ; }

const char * SymmetricTensor21::name() const
{ static const char n[] = "SymmetricTensor21" ; return n ; }

//----------------------------------------------------------------------

const AsymmetricTensor03 & AsymmetricTensor03::tag()
{ static const AsymmetricTensor03 self ; return self ; }

const char * AsymmetricTensor03::name() const
{ static const char n[] = "AsymmetricTensor03" ; return n ; }

//----------------------------------------------------------------------

const Matrix22 & Matrix22::tag()
{ static const Matrix22 self ; return self ; }

const char * Matrix22::name() const
{ static const char n[] = "Matrix22" ; return n ; }

//----------------------------------------------------------------------

const Matrix33 & Matrix33::tag()
{ static const Matrix33 self ; return self ; }

const char * Matrix33::name() const
{ static const char n[] = "Matrix33" ; return n ; }

//----------------------------------------------------------------------

} // namespace mesh
} // namespace stk

