=======================================
TriBITS Maintainers Guide and Reference
=======================================

:Author: Roscoe A. Bartlett (rabartl@sandia.gov)
:Date: |date|
:Version: .. include:: ../TribitsGitVersion.txt

.. |date| date::

:Abstract: This document describes the internal implementation and the maintenance of the TriBITS project itself.  The primary audience are those individuals who will make changes and contributions to the TriBITS project or just want to understand its implementation details.  Included is all of the same information as the TriBITS Users Guide but also include file names and line numbers for all of the documented TriBITS macros and functions.

.. sectnum::
   :depth: 2

.. Above, the depth of the TOC is set to just 2 because I don't want the
.. TriBITS function/macro names to have section numbers appearing before them.
.. Also, some of them are long and I don't want them to go off the page of the
.. PDF document.

.. Sections in this document use the underlines:
..
.. Level-1 ==================
.. Level-2 ------------------
.. Level-3 ++++++++++++++++++
.. Level-4 ..................

.. contents::


Introduction
=============

This document describes the usage and maintenance of the TriBITS (Tribal
Build, Integration, Test System) package itself.  This document includes a
super-set the material from the `TriBITS Users Guide and Reference`_ document.
In addition, all of the detailed function and macro documentation blocks
include the file names and line numbers where they are implemented in the
TriBITS source repository. This makes it easier to navigate around the TriBITS
source code when developing on TriBITS itself or just trying to understand its
implementation.


.. include:: ../TribitsGuidesBody.rst


.. include:: TribitsCoreDetailedReference.rst


TriBITS System Maintainers Documentation
========================================

This section contains more detailed information about the internal
implementation of TriBITS.  This information is meant to make it easier to
understand and manipulate the data-structures and macros/functions that make
up internal implementation of TriBITS and is important for `TriBITS System
Developers`_ and `TriBITS System Architects`_.

.. include:: ../../../core/package_arch/TribitsSystemDataStructuresMacrosFunctions.rst


.. include:: TribitsSystemMacroFunctionDoc.rst


.. include:: ../TribitsGuidesReferences.rst


.. include:: ../TribitsFAQ.rst


Appendix
========

.. include:: ../TribitsCMakeLanguageOverviewAndGotchas.rst

.. include:: ../TribitsHistory.rst

.. include:: ../TribitsPackageNotCMakePackage.rst

.. include:: ../TribitsDesignConsiderations.rst

.. include:: ../TribitsToolsDocumentation.rst
