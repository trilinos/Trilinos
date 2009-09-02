#include "Teuchos_ParameterList.hpp"
#include <boost/python.hpp>

using namespace boost::python;
using namespace Teuchos;


/******************************************************************/
// Type method: return python type of requested parameter.  Replaces
// templated C++ isType( ParameterList& self) methods.
PyObject * type( ParameterList& self,const std::string & name);
//////////////////////////////////////////////
// ** The following methods are added to ** //
// ** give ParameterList a PyDict "feel" ** //
//////////////////////////////////////////////

/******************************************************************/
// Comparison operators
int __cmp__( ParameterList& self,PyObject * obj);
// int __cmp__( ParameterList& self,const ParameterList & ParameterList) ;

/******************************************************************/
// Contains operator
int __contains__( ParameterList& self,const std::string & name) ;

/******************************************************************/
// Equals operators
PyObject * __eq__( ParameterList& self,PyObject * obj) ;

// PyObject * __eq__( ParameterList& self,const ParameterList & ParameterList) ;

/******************************************************************/
// __iter__ method
PyObject * __iter__( ParameterList& self) ;
/******************************************************************/
// Length operator
int __len__( ParameterList& self);
/******************************************************************/
// Not equals operators
PyObject * __ne__( ParameterList& self, PyObject * obj) ;
/******************************************************************/
// String representation method
PyObject * __repr__( ParameterList& self) ; 

/******************************************************************/
// String conversion method
PyObject * __str__( ParameterList& self) ;

/******************************************************************/
// Has_key method
int has_key( ParameterList& self,const std::string & name) ;
/******************************************************************/
// Items method
PyObject * items( ParameterList& self);

/******************************************************************/
// Iteritems method
PyObject * iteritems( ParameterList& self);
/******************************************************************/
// Iterkeys method
PyObject * iterkeys( ParameterList& self) ;
/******************************************************************/
// Itervalues method
PyObject * itervalues( ParameterList& self) ;
/******************************************************************/
// Keys method
PyObject * keys( ParameterList& self);
/******************************************************************/
// Update methods
void update( ParameterList& self,PyObject * dict, bool strict);

// void update( ParameterList& self,const ParameterList & ParameterList);

/******************************************************************/
// Values method
PyObject * values( ParameterList& self);
/******************************************************************/
// AsDict method: return a dictionary equivalent to the
// ParameterList
PyObject * asDict( ParameterList& self);


void _myprint(ParameterList& self,PyObject * pf, int indent, bool showTypes,
	      bool showFlags);


