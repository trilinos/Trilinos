// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Domi_getValidParameters.hpp"
#include "Domi_Utils.hpp"

#include "Teuchos_Array.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

using std::string;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Array;
using Teuchos::tuple;
using Teuchos::ParameterList;
using Teuchos::ParameterEntryValidator;
using Teuchos::EnhancedNumberValidator;
using Teuchos::ArrayNumberValidator;
//using Teuchos::ScalarOrArrayNumberValidator;
using Teuchos::StringToIntegralParameterEntryValidator;

namespace Domi
{

/** \brief Return a ParameterList filled with all valid parameters,
 *         documentation strings, default values, and validators
 */
RCP< const ParameterList > getValidParameters() 
{
  static RCP< const ParameterList > result;

  if ( result.is_null() )
  {

    // Allocate a new non-const ParameterList
    ParameterList* plist = new ParameterList;

    ////////////////////////////////////////////////////////////////
    // "comm dimensions" parameter applies to MDComm, MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    RCP< EnhancedNumberValidator< int > > axisCommNumber =
      rcp(new EnhancedNumberValidator< int >());
    axisCommNumber->setMin(-1);

    RCP< const EnhancedNumberValidator< int > > constAxisCommNumber =
      rcp_const_cast< EnhancedNumberValidator< int > >(axisCommNumber);

    RCP< const ParameterEntryValidator > axisCommValidator =
      rcp< const ParameterEntryValidator >
      (new ArrayNumberValidator< int >(constAxisCommNumber));

    Array< int > commDims(1);
    commDims[0] = -1;
    plist->set("comm dimensions",
               commDims, 
               "An array of ints that specifies the size of the "
               "MDComm along each axis. If the 'dimensions' parameter is "
               "present, then the length of that parameter determines the "
               "number of dimensions.  If 'dimensions' is not present, then "
               "the length of this parameter determines the number of "
               "dimensions.  If the length of this parameter is shorter than "
               "the number of dimensions, then this parameter is extended with "
               "values of -1.  A negative value tells Domi to fill in a "
               "logical value based on the total number of processors.",
               axisCommValidator);

    ////////////////////////////////////////////////////////////////
    // "periodic" parameter applies to MDComm, MDMap, and MDVector
    ////////////////////////////////////////////////////////////////
    RCP< EnhancedNumberValidator< int > > periodicNumber =
      rcp(new EnhancedNumberValidator< int >());
    periodicNumber->setMin(0);
    periodicNumber->setMax(1);

    RCP< const EnhancedNumberValidator< int > > constPeriodicNumber =
      rcp_const_cast< EnhancedNumberValidator< int > >(periodicNumber);

    RCP< const ParameterEntryValidator > periodicValidator =
      rcp< const ParameterEntryValidator >
      (new ArrayNumberValidator< int >(constPeriodicNumber));

    //int periodic = 0;
    Array< int > periodic(1);
    periodic[0] = 0;
    plist->set("periodic",
               periodic,
               "A scalar or an array of int flags specifying whether axes are "
               "periodic. If a scalar is given, then all axes share that "
               "periodicity flag.  If an array is given and it is shorter than "
               "the length of commDims array, then the unspecified "
               "entries are given a default value of zero (not "
               "periodic).",
               periodicValidator);

    ////////////////////////////////////////////////////////////////
    // "replicated boundary" parameter applies to MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    RCP< EnhancedNumberValidator< int > > repBndryNumber =
      rcp(new EnhancedNumberValidator< int >());
    periodicNumber->setMin(0);
    periodicNumber->setMax(1);

    RCP< const EnhancedNumberValidator< int > > constRepBndryNumber =
      rcp_const_cast< EnhancedNumberValidator< int > >(repBndryNumber);

    RCP< const ParameterEntryValidator > repBndryValidator =
      rcp< const ParameterEntryValidator >
      (new ArrayNumberValidator< int >(constRepBndryNumber));

    Array< int > repBndry(1);
    periodic[0] = 0;
    plist->set("replicated boundary",
               periodic,
               "A scalar or an array of int flags specifying whether periodic "
               "boundaries have a replicated boundary (true) or unique grid "
               "points (false). If a scalar is given, then all axes share that "
               "replicated boundary flag.  If an array is given and it is "
               "shorter than the length of commDims array, then the "
               "unspecified entries are given a default value of zero (no "
               "replicating boundaries).",
               repBndryValidator);

    ////////////////////////////////////////////////////////////////
    // "dimensions" parameter applies to MDComm, MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    RCP< EnhancedNumberValidator< dim_type > > dimensionNumber =
      rcp(new EnhancedNumberValidator< dim_type >());
     dimensionNumber->setMin(0);

    RCP< const EnhancedNumberValidator< dim_type > > constDimensionNumber =
    	rcp_const_cast< EnhancedNumberValidator< dim_type > >(dimensionNumber);

     RCP< const ParameterEntryValidator > dimensionValidator =
     	rcp< const ParameterEntryValidator >
     	(new ArrayNumberValidator< dim_type >(constDimensionNumber));

    Array< dim_type > dimensions(1);
    dimensions[0] = 0;
    plist->set("dimensions",
               dimensions,
               "An array of ordinals specifying the global dimensions of "
               "the MDMap. If present for the MDComm constructor, the length "
               "of this parameter will set the number of dimensions.  If not "
               "present for the MDComm constructor, the number of dimensions "
               "will be set by the length of the 'comm dimensions' parameter.",
               dimensionValidator);

    // Both boundary pad and communication pad use the same number and
    // array validators, so just construct one EnhancedNumberValidator
    // and one ArrayNumberValidator.
    int          pad = 0;
    Array< int > pads;

    RCP< EnhancedNumberValidator< int > > padNumberValidator =
      rcp(new EnhancedNumberValidator< int >());
    padNumberValidator->setMin(0);

    RCP< const EnhancedNumberValidator< int > > constPadNumberValidator =
      rcp_const_cast< EnhancedNumberValidator< int > >(padNumberValidator);

    RCP< const ParameterEntryValidator > padArrayValidator =
      rcp< const ParameterEntryValidator >
      (new ArrayNumberValidator< int >(constPadNumberValidator));

    ////////////////////////////////////////////////////////////////
    // "boundary pad size" parameter applies to MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    plist->set("boundary pad size",
               pad,
               "An int that specifies the boundary padding size for all axes.",
               padNumberValidator);

    ////////////////////////////////////////////////////////////////
    // "boundary pad sizes" parameter applies to MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    plist->set("boundary pad sizes",
               pads,
               "An array of ints specifying the size of the boundary padding "
               "along each axis. All unspecified entries take the value of "
               "the 'boundary pad size' parameter, which defaults to zero.",
               padArrayValidator);

    ////////////////////////////////////////////////////////////////
    // "communication pad size" parameter applies to MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    plist->set("communication pad size",
               pad,
               "An int that specifies the communication padding size for all "
               "axes.",
               padNumberValidator);

    ////////////////////////////////////////////////////////////////
    // "communication pad sizes" parameter applies to MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    plist->set("communication pad sizes",
               pads,
               "An array of ints specifying the size of the communication "
               "padding along each axis. All unspecified entries take the "
               "value of the 'communication pad size' parameter, which "
               "defaults to zero.",
               padArrayValidator);

    ////////////////////////////////////////////////////////////////
    // "layout" parameter applies to MDMap and MDVector
    ////////////////////////////////////////////////////////////////
    string layout = "Default";

    Array< string >
      layoutOpts(tuple(string("C Order"),
                       string("Fortran Order"),
                       string("Row Major"),
                       string("Column Major"),
                       string("Last Index Fastest"),
                       string("First Index Fastest"),
                       string("Default")));

    Array< string >
      layoutDocs(tuple(string("C storage order (last index varies "
                              "fastest)"),
                       string("Fortran storage order (first index "
                              "varies fastest)"),
                       string("Row major storage order (last index "
                              "varies fastest)"),
                       string("Column major storage order (first "
                              "index varies fastest)"),
                       string("Last index varies fastest"),
                       string("First index varies fastest"),
                       string("Fortran storage order")));

    Array< int > layoutVals(tuple(0, 1, 0, 1, 0, 1, 1));

    RCP< const ParameterEntryValidator > layoutValidator =
      rcp(new StringToIntegralParameterEntryValidator< int >
                   (layoutOpts(),
                    layoutDocs(),
                    layoutVals(),
                    string("Default"),
                    false));

    plist->set("layout",
               layout,
               "A string indicating how the data is laid out in memory. "
               "Default is currently set to Fortran order.",
               layoutValidator);

    ////////////////////////////////////////////////////////////////
    // "leading dimension" parameter applies to MDVector
    ////////////////////////////////////////////////////////////////
    plist->set("leading dimension",
               pad,
               "Use the 'leading dimension' parameter to specify multiple "
               "degrees of freedom at each MDMap index.  This increases the "
               "dimension of the MDVector by one, and the new degrees of "
               "freedom are accessed with the first index.",
               padNumberValidator);

    ////////////////////////////////////////////////////////////////
    // "trailing dimension" parameter applies to MDVector
    ////////////////////////////////////////////////////////////////
    plist->set("trailing dimension",
               pad,
               "Use the 'trailing dimension' parameter to specify multiple "
               "degrees of freedom at each MDMap index.  This increases the "
               "dimension of the MDVector by one, and the new degrees of "
               "freedom are accessed with the last index.",
               padNumberValidator);

    // ParameterList construction is done, so wrap it with an RCP<
    // const ParameterList >
    result.reset(plist);
  }

  return result; 
}

}
