#ifndef Teuchos_XMLOBJECT_H
#define Teuchos_XMLOBJECT_H

#include "Teuchos_XMLObjectImplem.hpp"
#include "Teuchos_Utils.hpp"

namespace Teuchos
{

	/** \ingroup XML 
	 * Representation of an XML data tree. XMLObject is a ref-counted
	 * handle to a XMLObjectImplem object, allowing storage by reference.
	 */

	class XMLObject
		{
		public:
			/** Constructors */
			//@{
			/** Empty ctor */
			XMLObject() : ptr_(null) {;}
			/** Create a node labeled by tag */
			XMLObject(const string& tag);

			/**
			 * construct with a pointer to the low-level representation. 
			 * This is used to allow construction of an XMLObject from the
			 * XMLObjectImplem* return value of ExceptionBase::toXML().
			 */
			XMLObject(XMLObjectImplem* ptr);
			//@}	

			/** \name Copying */
			//@{
			/** make a deep copy of this object */
			XMLObject deepCopy() const ;
			//@}

			/** Access to data stored in a node */
			//@{
			/** Return the tag of the current node */
			const string& getTag() const {return ptr_->getTag();}

			/** Find out if the current node has an attribute of the specified name */
			bool hasAttribute(const string& name) const 
				{return ptr_->hasAttribute(name);}

			/** Return the value of the attribute with the specified name */
			const string& getAttribute(const string& name) const 
				{return ptr_->getAttribute(name);}

			/** get an attribute, throwing an exception if it is not found */
			const string& getRequired(const string& name) const ;

			/** get a required attribute, returning it as a double */
			double getRequiredDouble(const string& name) const 
				{return atof(getRequired(name).c_str());}

			/** get a required attribute, returning it as a int */
			int getRequiredInt(const string& name) const 
				{return atoi(getRequired(name).c_str());}

			/** get a required attribute, returning it as a int */
			bool getRequiredBool(const string& name) const ;


			/** get an attribute, assigning a default value if the requested
			 * attribute does not exist */
			string getWithDefault(const string& name, 
														const string& defaultValue) const ;

			/** Return the number of child nodes owned by this node */
			int numChildren() const {return ptr_->numChildren();}

			/** Return the i-th child node */
			const XMLObject& getChild(int i) const {return ptr_->getChild(i);}

			/** Return the number of lines of character content stored in this node */
			int numContentLines() const {return ptr_->numContentLines();}

			/** Return the i-th line of character content stored in this node */
			const string& getContentLine(int i) const {return ptr_->getContentLine(i);}

			/** Represent this node and its children as a String */
			string toString() const {return ptr_->toString();}

			/** Write the header for this object to a string */
			string header() const {return ptr_->header();}

			/** Write the footer for this object to a string */
			string footer() const {return ptr_->footer();}

			/** Find out if a node is empty */
			bool isEmpty() const {return ptr_.get()==0;}

			/** Check that a tag is equal to an expected string */
			void checkTag(const string& expected) const ;
			//@}
	
			/** Methods for assembling a tree */
			//@{
			/** Add an attribute to the current node's atribute list */
			void addAttribute(const string& name, const string& value)
				{ptr_->addAttribute(name, value);}
			
			/** Add a double as an attribute */
			void addDouble(const string& name, double val)
				{addAttribute(name, Teuchos::toString(val));}

			/** Add an int as an attribute */
			void addInt(const string& name, int val)
				{addAttribute(name, Teuchos::toString(val));}

			/** Add a bool as an attribute */
			void addBool(const string& name, bool val)
				{addAttribute(name, Teuchos::toString(val));}
			
			/** Add a child node to the node */
			void addChild(const XMLObject& child)
				{ptr_->addChild(child);}
			/** Add a line of character content */
			void addContent(const string& contentLine)
				{ptr_->addContent(contentLine);}
			//@}
	
		private:
			RefCountPtr<XMLObjectImplem> ptr_;
		};

	/** \relates XMLObject write to stream */
	inline ostream& operator<<(ostream& os, const XMLObject& xml)
		{
			return os << xml.toString();
		}

	/** \relates XMLObject write to string */
	inline string toString(const XMLObject& xml)
		{
			return xml.toString();
		}

}
#endif

