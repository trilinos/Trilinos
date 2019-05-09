.. _vectors_overview:

Overview
########

|tpetra_mv|_ implements a finite-dimensional vector distributed over processes. ``Vector`` inherits from Tpetra's ``MultiVector`` class, which represents a collection of one or more vectors with the same ``Map``. Tpetra favors block algorithms, so it favors ``MultiVector``\s over single ``Vector``\s. A single Vector is just a ``MultiVector`` containing one vector, with a few convenience methods.

.. include:: /links.txt
