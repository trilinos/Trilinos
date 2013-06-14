==================================================================
Trilinos Configure, Build, Test, and Install Quick Reference Guide
==================================================================

:Author: Roscoe A. Bartlett
:Contact: bartlett.roscoe@gmail.com

:Abstract: This document contains quick reference information on how to configure, build, test, and install Trilinos using the TriBITS CMake build system.  The primary audience are users of Trilinos that need to configure and build the software.  The secondary audience are actual developers of Trilinos.

.. sectnum::

.. contents::

Introduction
============

Trilinos contains a large number of packages that can be enabled and there is a fairly complex dependency tree of required and optional package enables.  The following sections contain fairly generic information on how to configure, build, test, and install Trilinos that addresses a wide range of issues.

This is not the first document that a user should read when trying to set up to install Trilinos.  For that, see the INSTALL.* file.  There is a lot of information and activities mentioned in this quickref that most users (and even some Trilinos developers) will never need to know about.

Also, this particular quick reference has no information at all on what is actually in Trilinos.  For that, go to:

  http://trilinos.org

to get started.
