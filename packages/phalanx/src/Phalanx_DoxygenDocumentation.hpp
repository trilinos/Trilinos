#ifndef PHX_DOXYGEN_DOCUMENTATION_HPP
#define PHX_DOXYGEN_DOCUMENTATION_HPP

/*!

\mainpage

\section xpr_contents Contents

  - \ref xpr_introduction

  - \ref xpr_overview

  - \ref xpr_user_guide

  - \ref xpr_developer_guide

  - \ref xpr_faq

  - \ref xpr_authors 

\section xpr_introduction Introduction

\section xpr_overview Overview

The ??? is a package to ...

\section xpr_user_guide User's Guide

 - Requires the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> utilities library, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.
 
 - Requires the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.

 - Requires the <a href="http://www.boost.org">Boost Template Metaprogramming (MPL) Library</a>


\section xpr_developer_guide Developer's Guide

\section xpr_faq Frequently Asked Questions

  - Why is the scalar type not embedded in the FieldTag?  We designed the FieldTag to be scalar type independent.  The reason is that we could then use FieldTags as arguments for FieldEvaluator objects of ANY scalar type.  We have a factory that can automatically build FieldEvaluators for all scalar types.  This automation requires constructor arguments that are not dependent on the scalar type.  

\section xpr_authors Authors

  - Roger Pawlowski (PI), SNL 01414

*/
#endif
