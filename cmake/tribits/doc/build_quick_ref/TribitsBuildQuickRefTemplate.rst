=======================================================================
Generic TriBITS Project, Build, Test, and Install Quick Reference Guide
=======================================================================

:Author: Roscoe A. Bartlett
:Contact: bartlett.roscoe@gmail.com
:Date: |date|
:Version: .. include:: TribitsGitVersion.txt

.. |date| date::

:Abstract: This document is generated from the generic template body document ``TribitsBuildQickRefBody.rst`` and provides a general project-independent quick reference on how to configure, build, test, and install a project that uses the TriBITS CMake build system.  The primary audience of this particular build of this document are TriBITS project developers themselves.  A project-specific version of this document should be created and accessed by users of a particular TriBITS-based project.

.. sectnum::

.. contents::

Introduction
============

This document is created using the script ``create-build-quickref.sh`` in this
directory which just runs::

  $ ./create-project-build-quickref.py \
    --project-name="<Project>" \
    --project-template-file=TribitsBuildQuickRefTemplate.rst \
    --file-base=TribitsBuildQuickRef

to generate this document.  In a project-specific version, ``<Project>`` is
replaced with the actual project name (e.g. ``Trilinos``, or ``TriBITS``,
etc.).  This version of the generated document is referred to by the general
TribitsDeveloperGuide.[rst,html,pdf] document.

Below are given genetic versions of the sections that show up in every
project-specific build of this document.
