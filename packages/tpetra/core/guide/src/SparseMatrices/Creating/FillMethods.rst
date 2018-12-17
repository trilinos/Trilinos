.. _fill_methods:

Adding Entries to the CrsMatrix
###############################

Overview
========

Methods of ``CrsMatrix`` that start with "``insert``" actually change the
structure of the sparse matrix. Methods that start with "``replace``" or
"``sumInto``" only modify existing values.

Insertion into nonowned rows
============================

All methods (except for ``insertGlobalValues()`` and ``sumIntoGlobalValues();``)
that work with global indices only allow operations on indices owned by the
calling process. For example, methods that take a global row index expect that
row to be owned by the calling process. Access to nonowned rows, that is, rows
not owned by the calling process, requires performing an explicit communication
via the Import / Export capabilities of the ``CrsMatrix`` object. See the
documentation of DistObject_ for more details.

The methods ``insertGlobalValues()`` and ``sumIntoGlobalValues()`` are
exceptions to this rule. They both allows you to add data to nonowned rows.
These data are stored locally and communicated to the appropriate process on the
next call to ``globalAssemble()`` or ``fillComplete()``. This means that
``CrsMatrix`` provides the same nonowned insertion functionality that Epetra
provides via ``Epetra_FECrsMatrix``.

.. _DistObject: https://trilinos.org/docs/dev/packages/tpetra/doc/html/classTpetra_1_1DistObject.html
