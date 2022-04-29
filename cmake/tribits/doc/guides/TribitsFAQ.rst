FAQ
===

**Q:** Why does not TriBITS just use the standard CMake ``Find<PACKAGE_NAME>.cmake``
modules and the standard ``find_package()`` function to find TPLs?

**A:** The different "standard" CMake ``Find<PACKAGE_NAME>.cmake`` modules do not
have a standard set of outputs and therefore, can't be handled in a uniform
way.  For example, 
