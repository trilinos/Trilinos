// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
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

    // "axis comm sizes" parameter applies to MDComm
    RCP< EnhancedNumberValidator< int > > axisCommNumber =
      rcp(new EnhancedNumberValidator< int >());
    axisCommNumber->setMin(-1);

    RCP< const EnhancedNumberValidator< int > > constAxisCommNumber =
      rcp_const_cast< EnhancedNumberValidator< int > >(axisCommNumber);

    RCP< const ParameterEntryValidator > axisCommValidator =
      rcp< const ParameterEntryValidator >
      (new ArrayNumberValidator< int >(constAxisCommNumber));

    Array< int > axisCommSizes(1);
    axisCommSizes[0] = -1;
    plist->set("axis comm sizes",
               axisCommSizes, 
               "An array of ints whose length is the number of dimensions "
               "of the problem and whose entries specify the size of the "
               "MDComm along each axis. A negative value tells Domi to "
               "fill in a logical value based on the total number of "
               "processors",
               axisCommValidator);

    // "periodic" parameter applies to MDComm
    RCP< EnhancedNumberValidator< int > > periodicNumber =
      rcp(new EnhancedNumberValidator< int >());
    periodicNumber->setMin(0);
    periodicNumber->setMax(1);

    RCP< const EnhancedNumberValidator< int > > constPeriodicNumber =
      rcp_const_cast< EnhancedNumberValidator< int > >(periodicNumber);

    RCP< const ParameterEntryValidator > periodicValidator =
      rcp< const ParameterEntryValidator >
      //(new ScalarOrArrayNumberValidator< int >(constPeriodicNumber));
      (new ArrayNumberValidator< int >(constPeriodicNumber));

    //int periodic = 0;
    Array< int > periodic(1);
    periodic[0] = 0;
    plist->set("periodic",
               periodic,
               "A scalar or an array of int flags specifying whether axes are "
               "periodic. If a scalar is given, then all axes share that "
               "periodicity flag.  If an array is given and it is shorter than "
               "the length of axisCommSizes array, then the unspecified "
               "entries are given a default value of zero (not "
               "periodic).",
               periodicValidator);

    // "dimensions" parameter applies to MDMap
    RCP< EnhancedNumberValidator< long int > > dimensionNumber =
      rcp(new EnhancedNumberValidator< long int >());
     dimensionNumber->setMin(0);

    RCP< const EnhancedNumberValidator< long int > > constDimensionNumber =
    	rcp_const_cast< EnhancedNumberValidator< long int > >(dimensionNumber);

     RCP< const ParameterEntryValidator > dimensionValidator =
     	rcp< const ParameterEntryValidator >
     	(new ArrayNumberValidator< long int >(constDimensionNumber));

    Array< long int > dimensions(1);
    dimensions[0] = 0;
    plist->set("dimensions",
               dimensions,
               "An array of ordinals specifying the global dimensions of "
               "the MDMap. The length of this array should be the same as "
               "the length of the axisCommSizes array. Note that this is "
               "an array of long int and it will need to be copied to an "
               "array of type GlobalOrd.",
               dimensionValidator);

    // Both boundary pad and communication pad use the same validator,
    // so just construct one
    Array< int > pad;

    RCP< EnhancedNumberValidator< int > > padNumber =
      rcp(new EnhancedNumberValidator< int >());
    padNumber->setMin(0);

    RCP< const EnhancedNumberValidator< int > > constPadNumber =
      rcp_const_cast< EnhancedNumberValidator< int > >(padNumber);

    RCP< const ParameterEntryValidator > padValidator =
      rcp< const ParameterEntryValidator >
      //(new ScalarOrArrayNumberValidator< int >(constPadNumber));
      (new ArrayNumberValidator< int >(constPadNumber));

    // "boundary pad" parameter applies to MDMap
    plist->set("boundary pad",
               pad,
               "A scalar or an array of ints specifying the size of the "
               "boundary padding. If a scalar is given, then all axes share "
               "that padding value.  An array can be used to specify padding "
               "along each axis.  All unspecified entries are assumed to be "
               "zero.",
               padValidator);

    // "communication pad" parameter applies to MDMap
    plist->set("communication pad",
               pad,
               "A scalar or an array of ints specifying the size of the "
               "communication padding. If a scalar is given, then all axes "
               "share that padding value.  An array can be used to specify "
               "padding along each axis.  All unspecified entries are assumed "
               "to be zero.",
               padValidator);

    // "layout" parameter applies to MDMap
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

    // ParameterList construction is done, so wrap it with an RCP<
    // const ParameterList >
    result.reset(plist);
  }

  return result; 
}

}
