// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// @HEADER


#include "Teuchos_VerbosityLevelCommandLineProcessorHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace {


using Teuchos::Array;
using Teuchos::tuple;


const Array<Teuchos::EVerbosityLevel>
verbosityLevelValues = tuple<Teuchos::EVerbosityLevel>(
  Teuchos::VERB_DEFAULT,
  Teuchos::VERB_NONE,
  Teuchos::VERB_LOW,
  Teuchos::VERB_MEDIUM,
  Teuchos::VERB_HIGH,
  Teuchos::VERB_EXTREME
  );


const Array<std::string>
verbosityLevelNamesStorage = tuple<std::string>(
  toString(Teuchos::VERB_DEFAULT),
  toString(Teuchos::VERB_NONE),
  toString(Teuchos::VERB_LOW),
  toString(Teuchos::VERB_MEDIUM),
  toString(Teuchos::VERB_HIGH),
  toString(Teuchos::VERB_EXTREME)
  );


Array<const char*> verbosityLevelNames;
// The above variable is Intialized below in order so that if exceptions are
// thrown then they will be caught in main()!


} // namespace


void Teuchos::setVerbosityLevelOption(
  const std::string &optionName,
  EVerbosityLevel *verbLevel,
  const std::string &docString,
  CommandLineProcessor *clp,
  const bool required
  )
{
  const int numVerbLevels = verbosityLevelValues.size();

  if ( !verbosityLevelNames.size() ) {
    verbosityLevelNames = tuple<const char*>(
      verbosityLevelNamesStorage[0].c_str(),
      verbosityLevelNamesStorage[1].c_str(),
      verbosityLevelNamesStorage[2].c_str(),
      verbosityLevelNamesStorage[3].c_str(),
      verbosityLevelNamesStorage[4].c_str(),
      verbosityLevelNamesStorage[5].c_str()
      );
  }
  
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT( optionName.length() );
  TEUCHOS_ASSERT( verbLevel );
  TEUCHOS_ASSERT( clp );
  TEUCHOS_ASSERT( implicit_cast<int>(verbosityLevelNamesStorage.size()) == numVerbLevels );
  TEUCHOS_ASSERT( implicit_cast<int>(verbosityLevelNames.size()) == numVerbLevels );
#endif
  clp->setOption(
    optionName.c_str(), verbLevel,
    numVerbLevels, &verbosityLevelValues[0], &verbosityLevelNames[0],
    docString.c_str(), required
    );
}
