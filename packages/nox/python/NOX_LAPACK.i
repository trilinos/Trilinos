// -*- c++ -*-

%module(package="NOX") LAPACK

// This is a trick to get swig to wrap LOCA correctly.  Swig apparently gets
// confused when creating a module that imports a module of the same name
// (here the LAPACK module) and drops the namespace qualifications in the 
// python shadow classes.  To get around this, we put all of the interface
// definition in NOX_LAPACK.swi except the module name and include it here.
// We then define another interface file, NOX_LAPACK_import.i which also
// includes NOX_LAPACK.swi, but changes the module name to NOX.LAPACK.
// NOX_LAPACK_import.i is the inteface file that should be imported in LOCA

%include "NOX_LAPACK.swi"
