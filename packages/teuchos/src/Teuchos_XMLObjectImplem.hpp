#ifndef TEUCHOS_XMLOBJECTIMPLEM_H
#define TEUCHOS_XMLOBJECTIMPLEM_H

#include "Teuchos_Hashtable.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace Teuchos
{

	class XMLObject;

	/** 
	 * XMLObjectImplem is the low-level implementation of XMLObject
	 */
	class XMLObjectImplem
		{
		public:
			/** construct with a tag string */
			XMLObjectImplem(const string& string);

			/** deep copy */
			XMLObjectImplem* deepCopy() const ;

			/** add a [name, value] attribute */
			void addAttribute(const string& name, const string& value);

			/** add a child XMLObject */
			void addChild(const XMLObject& child);

			/** add a content line */
			void addContent(const string& contentLine);

			/** write as a string */
			string toString() const ;

			/** return the tag string */
			const string& getTag() const {return tag_;}

			/** Write the header */
			string header() const ;

			/** Write the footer */
			string footer() const {return "</" + getTag() + ">";}

			/** determine whether an attribute exists */
			bool hasAttribute(const string& name) const 
				{return attributes_.containsKey(name);}

			/** look up an attribute by name */
			const string& getAttribute(const string& name) const 
				{return attributes_.get(name);}

			/** return the number of children */
			int numChildren() const ;

			/** look up a child by its index */
			const XMLObject& getChild(int i) const ;

			/** get the number of content lines */
			int numContentLines() const {return content_.length();}

			/** look up a content line by index */
			const string& getContentLine(int i) const {return content_[i];}

		private:
			string tag_;
			Hashtable<string, string> attributes_;
			Array<XMLObject> children_;
			Array<string> content_;
		};

}

#endif

