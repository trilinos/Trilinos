Why a TriBITS Package is not a CMake Package
--------------------------------------------

Note that a `TriBITS Package`_ is not the same thing as a "Package" in raw
CMake terminology.  In raw CMake, a "Package" is some externally provided bit
of software or other utility for which the current CMake project has an
optional or required dependency (see `CMake: How to Find Libraries
<http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries>`_).  Therefore, a raw
CMake "Package" actually maps to a `TriBITS TPL`_.  A raw CMake "Package"
(e.g. Boost, CUDA, etc.)  can be found using a standard CMake find module
``Find<rawPackageName>.cmake`` using the built-in CMake command
``find_package(<rawPackageName>)``.  It is unfortunate that the TriBITS and
the raw CMake definitions of the term "Package" are not exactly the same.
However, the term "Package" was coined by the Trilinos project long ago before
CMake was adopted as the Trilinos build system and Trilinos' definition of
"Package" (going back to 1998) pre-dates the development of CMake (see
`History of CMake <http://en.wikipedia.org/wiki/CMake#History>`_) and
therefore Trilinos dictated the terminology of TriBITS and the definition of
the term "Package" in the TriBITS system.  However, note that both meanings of
the term "Package" are consistent with the more general software engineering
definition of a "Package" according to `Software Engineering Packaging
Principles`_.
