# @HEADER
# ************************************************************************
# 
#              PyTrilinos: Python Interface to Trilinos
#                 Copyright (2004) Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Bill Spotz (wfspotz@sandia.gov) 
# 
# ************************************************************************
# @HEADER
import sys

"""
PyTrilinos.Teuchos is the python interface to Trilinos package
Teuchos:

    http://trilinos.sandia.gov/packages/teuchos

The purpose of Teuchos is to provide a number of utilities often
needed by numerical applications, but that are not necessarily
numerical by nature.  The python version of the Teuchos package
supports the following classes:

    * ParameterList           - List of arbitrarily-typed values,
                                keyed by strings
    * XMLObject               - Object-oriented interface to XML
                                objects
    * XMLParameterListReader  - ParameterList input from XML
    * XMLParameterListWriter  - ParameterList output to XML
    * XMLInputSource          - Base class for converting a stream
                                to XML
    * FileInputSource         - Class for converting file contents
                                to XML
    * StringInputSource       - Class for converting string contents
                                to XML
    * ScalarTraits            - Function factory for ScalarTraits<...>
                                classes
    * Time                    - Wall-clock timer class

The ParameterList class matches string keys to arbitrarily-typed
values.  In python, the Teuchos.ParameterList is tightly integrated
with python dictionaries -- PyTrilinos methods that expect a
ParameterList will accept a python dictionary.
"""

from _Teuchos import *

__version__ = Teuchos_Version().split()[2]


def ScalarTraits(typecode):
    if typecode == 'd':
        return ScalarTraitsDouble()
    if typecode == 'f':
        return ScalarTraitsFloat()
    else:
        raise NotImplementedError("type '%s' is not implemented")
    

class ParameterList(dict):
    
    def __init__(self,name="ANONYMOUS"):
        if isinstance(name, dict):
            dict.__init__(self,**name)
        elif isinstance(name, str):
            self._name = name
        else:
            raise ValueError("bad constructor")
        dict.__init__(self)
        
    def name(self):
        return self._name
    def setName(self,name):
        self._name = name
    
    def __repr__(self):
        return "ParameterList(%s)" %dict.__repr__(self)
    def __str__(self):
        return "Teuchos.ParameterList(%s)" %dict.__str__(self)
    
    def isParameter(self,item):
        return self.has_key(item)
    
    def isSublist(self,item):
        if  not self.has_key(item):
            return False
        if isinstance(self[item], dict):
            return True
        else:
            return False
    
    def sublist(self,name):
        return self.setdefault(name,ParameterList(name))
    
    def setParameters(self,params):
        return self.update(params)
        
    def __getitem__(self,item):
        if not isinstance(item , str):
            raise TypeError("nsuported type for plist")
        return dict.__getitem__(self,item)
    def __setitem__(self,item,val):
        
        if not isinstance(item, str):
            raise TypeError("unsuported type for plist, item must be a string")
        inst = lambda x: isinstance(val, x)
        or_ = lambda x,y: x or y
        any  = lambda x : bool(reduce(or_, x,False))
        
        valid = [int,float,str,dict]
        bool_list = map( inst,valid )
        
        if not any(bool_list):
            raise TypeError("unsuported type for plist")
            
        return dict.__setitem__(self,item, val)

    def type(self,name):
        return type( self[name] )
    
    def asDict(self):
        return dict(**self)
        
    def currentParametersString(self):
        return str(self)
    
    def set(self,item,val):
        return self.__setitem__(item, val)
    
    def _print(self,f=sys.stdout,indent=0,showTypes=False,showFlags=True):
        print >> f, str(self)
    
    def unused(self,f):
        return []
