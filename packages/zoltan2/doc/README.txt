NOTES:
-------------

To build the doxygen documentation run:

build_docs

You may want to edit this Doxyfile line:
HAVE_DOT               = YES
to "NO" if you don't have the dot tool.  But getting it
is a good idea, because it is used to generate inheritance
diagrams.
-------------

To read the doxygen documentation, begin at:

zoltan2/doc/html/index.html

-------------
Our Doxygen input files have the suffix "dox" rather than "doc".  This
is because of complaints for users of software that thinks all .doc
files are MS Word files and opens them with Word.



