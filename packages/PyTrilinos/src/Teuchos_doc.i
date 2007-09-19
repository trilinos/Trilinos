// -*- C++ -*-
// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
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
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

// The purpose of this file is to provide manually-written
// documentation strings for PyTrilinos.Teuchos functions, classes and
// methods, when the doxygen- or swig-generated documentation is
// either insufficient or inaccurate.

// Teuchos.ParameterList class

%feature("docstring") Teuchos::ParameterList
"The ``ParameterList`` class is an important utility class that is used
by several Trilinos packages for communicating arbitrary-type
parameters between users and packages.

Often, the ``ParameterList`` class is invisible to the user.  It is
analagous to the python dictionary (with the restriction that the
dictionary keys must be strings), and python programmers can provide a
python dictionary wherever a ``ParameterList`` is expected.
``Teuchos`` is imported by other packages that use the
``ParameterList`` class and converts between dictionaries and
``ParameterList`` objects automatically.

The user can create a ``Teuchos.ParameterList`` directly, using the
constructor, ``set`` and ``sublist`` methods, if he so chooses, and
methods that accept ``ParameterList`` objects will work as expected.
It is really just a question of verbosity and elegance that argues in
favor of using a python dictionary.

The python implementation of the ``ParameterList`` class has been
expanded extensively.  Its constructor can accept a python dictionary,
and several methods and operators have been added to the class so that
it behaves somewhat like a dictionary.

C++ ``ParameterList`` objects are designed to support parameters of
arbitrary type.  The python implementation supports a subset of types
*a priori* :

  +-------------------------+-----+-------------------+
  |       Python type       | Dir |     C/C++ type    |
  +-------------------------+-----+-------------------+
  | ``bool``                | <-> | ``bool``          |
  +-------------------------+-----+-------------------+
  | ``int``                 | <-> | ``int``           |
  +-------------------------+-----+-------------------+
  | ``float``               | <-> | ``double``        |
  +-------------------------+-----+-------------------+
  | ``str``                 | <-- | ``char *``        |
  +-------------------------+-----+-------------------+
  | ``str``                 | <-> | ``std::string``   |
  +-------------------------+-----+-------------------+
  | ``dict``                | --> | ``ParameterList`` |
  +-------------------------+-----+-------------------+
  | ``ParameterList``       | <-> | ``ParameterList`` |
  +-------------------------+-----+-------------------+

The C++ ``ParameterList`` class supports ``begin()`` and ``end()``
methods for iterating over the parameters.  These methods are disabled
in the python implementation, in favor of the dictionary iterator
methods: ``__iter__()``, ``iteritems()``, ``iterkeys()`` and
``itervalues()``.

Note that the C++ implementation of the ``ParameterList`` class does
not support parameter deletion.  Therefore, python dictionary methods
that delete items, such as ``pop()`` or ``__delitem__()``, have not
been added to the ``ParameterList`` class."

%feature("docstring")
Teuchos::ParameterList::ParameterList
"If ``dict`` is provided, it must be a dictionary whose keys are all
strings and whose values are all of supported types.  The string name
argument is optional and defaults to ``ANONYMOUS``."

%feature("docstring")
Teuchos::ParameterList::set
"The templated C++ ``set()`` method is replaced in python with a method
that takes a string name and a python object of supported type.  For
example::

    plist = Teuchos.ParameterList()
    plist.set('b',True)
    plist.set('i',10)
    plist.set('f',2.718)
    plist.set('s','Trilinos')
    plist.set('d',{'a':1, 'b':2})"

%feature("docstring")
Teuchos::ParameterList::get
"The templated C++ ``get()`` method is replaced in python with a method
that returns a python object.  For example::

    plist = Teuchos.ParameterList({'f':2.718, 'd':{'a':1. 'b':2}})
    print plist.get('f')
    print plist.get('d')

  will output::

    2.718
    {'a': 1, 'b': 2}"

%feature("docstring")
Teuchos::ParameterList::setParameters
"The ``setParameters()`` method can take either a ``ParameterList`` or
a python dictionary as its argument.  The ``ParameterList`` is updated
to contain all of the entries of the argument."

%feature("docstring")
Teuchos::ParameterList::_print
"Since ``print`` is a python keyword, the ``print()`` C++ method has
been renamed ``_print`` in python.  It takes an optional file argument
that defaults to standard output.  Its output is the same as the C++
implementation."

%feature("autodoc",
"__str__(self) -> string

The ``__str__()`` method returns a string representation of the
``ParameterList`` as though it were a python dictionary.  The python
``eval`` function applied to the output of ``__str__()`` will produce
an equivalent dictionary.")
Teuchos::ParameterList::__str__;

%feature("autodoc",
"__repr__(self) -> string

The ``__repr__()`` method returns the ``__str__()`` output
encapsulated by ``ParameterList(...)``.  The python ``eval`` function
applied to the output of ``__repr__()`` will produce an equivalent
``ParameterList``.")
Teuchos::ParameterList::__repr__;

%feature("autodoc",
"unused(self, pf=None)

The ``unused()`` method in python takes an optional python file object
as its argument, defaulting to standard output.  This specifies where
the output should go.")
Teuchos::ParameterList::unused;

%feature("docstring") Teuchos::ParameterList::type
"Parameter type determination.  With the templated ``isType()`` methods
disabled, type determination of python ``ParameterList`` entries is
accomplished with the ``type()`` method, which returns python type
objects.  For example::

  plist = Teuchos.ParameterList({'b':True, 's':'Trilinos'})
  print plist.type('b')
  print plist.type('s')

results in::

  <type 'bool'>
  <type 'str'>

A non-existent key given as the argument will raise a ``KeyError``
exception."

%feature("docstring")
Teuchos::ParameterList::asDict
"The ``asDict()`` method has been added to ``ParameterList``, which
returns the contents of the ``ParameterList`` converted to a python
dictionary.
"

%feature("docstring")
Teuchos::ParameterList::__eq__
"The ``ParameterList`` equals operator (==)"

%feature("docstring")
Teuchos::ParameterList::__ne__
"The ``ParameterList`` not equals operator (!=)"

%feature("docstring")
Teuchos::ParameterList::__contains__
"The python ``in`` operator works for ``ParameterList`` objects,
searching the parameter names::

  plist = Teuchos.ParameterList({'b':False})
  print 'a' in plist
  print 'b' in plist

produces::

  False
  True"

%feature("docstring")
Teuchos::ParameterList::__len__
"The python ``len()`` function works on ``ParameterList`` objects just
as on python dictionaries::

  plist = Teuchos.ParameterList({'b':True,
                                 'i':10,
                                 'f':2.718,
                                 's':'Trilinos',
                                 'd':{'a':1, 'b':2}})
  print len(plist)

gives::

  5"

%feature("docstring")
Teuchos::ParameterList::__iter__
"To iterate over the parameters in a ``ParameterList``, treat it like a
dictionary::

  plist = Teuchos.ParameterList({'b':True,
                                 'i':10,
                                 'f':2.718,
                                 's':'Trilinos',
                                 'd':{'a':1, 'b':2}})
  for key in plist:
    print key, ':', plist[key]

will result in the output::

  b : True
  d : {'a': 1, 'b': 2}
  f : 2.718
  i : 10
  s : Trilinos

Note that the order of the parameters is somewhat indeterminant, as
with dictionaries, because the iteration object is obtained from an
equivalent dictionary, and dictionaries are ordered by hash
function."

%feature("docstring")
Teuchos::ParameterList::__setitem__
"Like dictionaries, parameters can be set using square brackets::

  plist = Teuchos.ParameterList()
  plist['zero'] = 0"

%feature("docstring")
Teuchos::ParameterList::__getitem__
"Like dictionaries, parameters can be gotten using square brackets::

  plist = Teuchos.ParameterList()
  plist['f'] = 2.718
  e = plist['f']"

%feature("docstring")
Teuchos::parameterList::update
"An ``update()`` method has been added to the ``ParameterList`` class
that can accept either a ``ParameterList`` or a python dictionary.
Otherwise, it behaves just as the dictionary method, which is
functionally equivalent to the ``setParameters()`` method."

%feature("docstring")
Teuchos::ParameterList::has_key
"Equivalent to the python dictionary has_key() method"

%feature("docstring")
Teuchos::ParameterList::items
"Equivalent to the python dictionary items() method"

%feature("docstring")
Teuchos::ParameterList::iteritems
"Equivalent to the python dictionary iteritems() method"

%feature("docstring")
Teuchos::ParameterList::iterkeys
"Equivalent to the python dictionary iterkeys() method"

%feature("docstring")
Teuchos::ParameterList::itervalues
"Equivalent to the python dictionary itervalues() method"

%feature("docstring")
Teuchos::ParameterList::keys
"Equivalent to the python dictionary keys() method"

%feature("docstring")
Teuchos::ParameterList::values
"Equivalent to the python dictionary values() method"


// Teuchos.XMLObject class

%feature("docstring")
Teuchos::XMLObject
"Representation of an XML data tree."

%feature("docstring")
Teuchos::XMLObject::XMLObject
"The constructor that takes an ``XMLObjectImplem*`` argument has been
removed.  The ``XMLObjectImplem`` class is hidden from the python
user."

%feature("docstring")
Teuchos::XMLObject::__str__
"The ``__str__()`` method is provided so that it is possible to
``print`` an ``XMLObject`` object.  It returns the same string as the
``toString()`` method, but if ``toString()`` raises an exception (such
as when the ``XMLObject`` is empty), the ``__str__()`` method returns
the empty string."
