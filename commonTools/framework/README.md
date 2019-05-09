# Trilinos Framework Tools Tests

This is a collection of tools used by the Trilinos framework.  This directory
is also a TriBITS package called `TrilinosFrameworkTests`.  The purpose of
this package is to test framework tools.  No other Trilinos TriBITS packages
should have a dependency on this TribitsPackage.  Tools that are used by
Trilinos packages should be put into another Trilinos TriBITS packages, like
an existing or new Teuchos subpackage.
