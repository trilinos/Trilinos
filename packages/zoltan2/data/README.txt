The file parameters.xml contains the XML definition of all the
Zoltan2 parameters and validators.  New parameters should be added
to this data file.  When CMake runs, directives in
zoltan2/src/CMakeLists.txt will create a header file which will
be used at compile time to define an XML string in the
library containing the parameters.

See zoltan2/doc/developer.dox for information on how to add
a parameter.

You can use the python script in this directory to do some basic
validity checks on paramters.xml after you have edited it.
