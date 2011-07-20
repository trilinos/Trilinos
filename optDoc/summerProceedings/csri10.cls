%%
%% This is file `csri10.cls'
%% CSRI Summer Papers Macros.  For use with siamltex
%%
%% By Chris Siefert
%% (Modified by Mike Parks)
%% (Modified by Denis Ridzal)
%% (Modified by Eric Cyr)
%% Sandia National Labs
%%
%% May 12, 2010
%%
%% Version 1.1
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{csri10}[2010/05/12 CSRI 2010 Summer Report Class]

% Options for siamltex
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{siamltex}}

% Process Options
\ProcessOptions

% Load siamltex
\LoadClass{siamltex}

% Table of contents
\renewcommand\contentsname{Table of Contents}
\setcounter{tocdepth}{0}
\newcommand\tableofcontents{%
\noindent{\bf\contentsname}
%    \section*{\contentsname
%        \@mkboth{%
%           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
%           ~%
           \bigskip%
    \@starttoc{toc}%
    }
\newcommand*\l@chapter{\@dottedtocline{0}{1.5em}{2.3em}}
\newcommand*\l@section{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand*\l@subsection{\@dottedtocline{2}{3.8em}{3.2em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{10em}{5em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{12em}{6em}}

% CSRI Title Page: modified from siamltex.cls
\def\ps@csriheadings{%
  \let\@oddfoot\@empty\let\@evenfoot\@empty
  \def\@defaulthead{ 
    \rlap{\footnotesize CSRI Summer Proceedings 2010}\hfil\llap{\thepage}
  }% end default headings   
  \def\@oddhead{\@defaulthead}
  \def\@evenhead{\@defaulthead} 
  \let\@mkboth\@gobbletwo
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}

% New \maketitle.  Horked + slightly modified from siamltex.cls
\renewcommand\maketitle{\par
  \begingroup
  \renewcommand{\thefootnote}{\fnsymbol{footnote}}%
  \def\@makefnmark{\hbox to\z@{$\m@th^{\@thefnmark}$\hss}}%
  \long\def\@makefntext##1{\parindent 1em\noindent
    \hbox to1.8em{\hss$\m@th^{\@thefnmark}$}##1}%
  \if@twocolumn
  \ifnum \col@number=\@ne
  \@maketitle
  \else
  \twocolumn[\@maketitle]%
  \fi
  \else
  \global\@topnum\z@   % Prevents figures from going at top of page.
  \@maketitle
  \fi
  \thispagestyle{csriheadings}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%
  %\let\maketitle\relax\let\@maketitle\relax
  \gdef\@thanks{}\gdef\@author{}\gdef\@title{}
  %\let\thanks\relax
}

\endinput
