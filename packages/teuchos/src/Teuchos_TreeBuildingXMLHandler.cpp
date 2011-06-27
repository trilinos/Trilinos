// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_TreeBuildingXMLHandler.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_TestForException.hpp"

using namespace Teuchos;

TreeBuildingXMLHandler::TreeBuildingXMLHandler()
	: root_(), current_(), path_()
{ 
	current_ = root_;
}

void TreeBuildingXMLHandler::characters(const std::string& chars)
{
  TEST_FOR_EXCEPTION(current_.isEmpty(), std::logic_error,
                     "TreeBuildingXMLHandler::trying to add content to an empty node");
                        
  current_.addContent(StrUtils::fixUnprintableCharacters(chars));
}

void TreeBuildingXMLHandler::startElement(const std::string& tag, 
																					const Map& attributes)
{
  XMLObject parent;
  
  if (current_.isEmpty()) 
    {
      root_ = XMLObject("root");
      current_ = root_;
    }
  parent = current_;
  path_.push(current_);
  current_ = XMLObject(tag);
  parent.addChild(current_);

  for (Map::const_iterator i=attributes.begin(); i != attributes.end(); ++i)
    {
      const std::string& key = (*i).first;
      const std::string& val = (*i).second;
      current_.addAttribute(key, val);
    }
}

int TreeBuildingXMLHandler::endElement(const std::string& tag)
{
  int error = 0;
  if (path_.size() > 0)
    {
      if (current_.getTag() != tag) {
        error = 1; // error: tags must be balanced
      }
      current_ = path_.top();
      path_.pop();
    }
  else 
    {
      error = 1; // error: cannot end element that wasn't started
    }
  return error;
}

