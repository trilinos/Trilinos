.. _index:

Tpetra: Parallel Sparse Linear Algebra
######################################

Petra, derived from the Greek word "petros" meaning "stone, rock", is an object
model describing distributed sparse linear algebra objects, such as sparse
matrices, vectors, along with computational kernels, independent of language and
implementation.

.. image:: /Images/petra.jpg

Image credit: [#pic]_.

There are two implementations of the Petra object model in Trilinos_: Epetra_, for "Essential Petra", and Tpetra_, for "Templated Petra".  Epetra was the first implementation of the Petra object model and is still the most heavily used.  Epetra is written in :math:`\text{C++}\leq 1998` and provides a stable interface accessible to C & Fortran users. Tpetra_ is the newest implementation of the Petra object model and was born of the desire to take advantage of new software and hardware features that would break backwards compatibility in Epetra.  Among the features offered by Tpetra are its support for problems greater than 2 billion unknowns with arbitrary and mixed precision,  and its adoption of hybrid parallelism (MPI+X) through the Kokkos_ package.  Many Trilinos packages and applications produce, modify, or consume Tpetra’s linear algebra objects, or depend on Tpetra’s parallel data redistribution facilities.

Indices and tables
==================

* :ref:`genindex`

* :ref:`search`

.. toctree::
   :maxdepth: 3
   :hidden:
   :numbered: 4

   FrontMatter/Index
   Obtaining/Index
   Initializing/Index
   Maps/Index
   Vectors/Index
   SparseMatrices/Index
   DataRedist/Index
   Examples/Index

.. include:: /links.txt

.. [#pic] http://www.nationalgeographic.com/archaeology-and-history/archaeology/lost-city-petra/
