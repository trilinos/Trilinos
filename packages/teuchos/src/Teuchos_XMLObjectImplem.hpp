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

#ifndef TEUCHOS_XMLOBJECTIMPLEM_H
#define TEUCHOS_XMLOBJECTIMPLEM_H

/*! \file Teuchos_XMLObjectImplem.hpp
    \brief Low level implementation of XMLObject
*/

#include "Teuchos_map.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos
{

	class XMLObject;

	/** 
	 * \brief The XMLObjectImplem class takes care of the low-level implementation 
		details of XMLObject
	 */
	class XMLObjectImplem
		{
      typedef Teuchos::map<string, string> Map;

		public:
			//! Construct with a tag string 
			XMLObjectImplem(const string& string);

			//! Deep copy
			XMLObjectImplem* deepCopy() const ;

			//! Add a [name, value] attribute
      void addAttribute(const string& name, const string& value);

			//! Add a child XMLObject
			void addChild(const XMLObject& child);

			//! Add a content line
			void addContent(const string& contentLine);

			//! Write as a string
			string toString() const ;

			//! Return the tag string
			const string& getTag() const {return tag_;}

			//! Write the header
			string header() const ;

			//! Write the header terminated as <Header/>
			string terminatedHeader() const ;

			//! Write the footer
			string footer() const {return "</" + getTag() + ">";}

			//! Determine whether an attribute exists
			bool hasAttribute(const string& name) const 
				{return attributes_.find(name) != attributes_.end();}

			//! Look up an attribute by name
			const string& getAttribute(const string& name) const 
      {return (*(attributes_.find(name))).second;}

			//! Return the number of children
			int numChildren() const ;

			//! Look up a child by its index
			const XMLObject& getChild(int i) const ;

			//! Get the number of content lines
			int numContentLines() const {return content_.length();}

			//! Look up a content line by index
			const string& getContentLine(int i) const {return content_[i];}

      //!  Print to stream with the given indentation level
      void print(ostream& os, int indent) const ;

		private:

      //! Print content lines using the given indentation level
      void printContent(ostream& os, int indent) const ;

			string tag_;
			Map attributes_;
			Array<XMLObject> children_;
			Array<string> content_;
		};

}

#endif

