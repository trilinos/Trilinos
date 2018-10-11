Creating a CrsMatrix
####################

.. toctree::
   :maxdepth: 2
   :hidden:

   Construct
   FillMethods
   FillComplete

Creating and adding entries ("filling") of a ``Tpetra::CrsMatrix`` involves the
following steps:

* :ref:`Construct the CrsMatrix <crsmatrix_construct>` (by calling one of its constructors)
* :ref:`Call methods to add entries <fill_methods>` to the sparse matrix
* :ref:`Call the matrix's <fill_complete>` ``fillComplete()`` method
