%  American Nuclear Society Transaction template
%  by Seth R. Johnson, based on the ANS transaction guidelines
%
% For further updates, see http://reference-man.com/projects.html#anstemplate
%
% Usage: \documentclass{anstrans}
%  Define \title and \author as usual, and also define \institute and \email .
%  For two-column tables, use the table* environment. Multiple authors should
%  look like 
%   {Seth R. Johnson \and Edward W. Larsen}
%
% v1.1.3: Improved header behavior
% v1.1.2: Improved spacing for standards compliance
%         now commas separate author names
%         Added support for proper equation numbering in appendices
% v1.1.1: Figure and table captions are now flush with the left as the
%         standards specify. Now uses the ``txfonts'' package for much improved
%         Times-compatible math font.
%

\NeedsTeXFormat{LaTeX2e}[1995/12/01]

\ProvidesClass{anstrans}[2011/06/07 v1.1.2 ANS Transaction class]

%%%%%%%%%%% TYPE AND GEOMETRY %%%%%%%%%%%
\LoadClass[twocolumn,10pt]{article}

\pagestyle{empty}						%no page numbering

\RequirePackage[text={7in,9in},centering]{geometry}

% ams math must be loaded before txfonts
\usepackage{amsmath}
% Use Times for math font and text font.
\RequirePackage[T1]{fontenc}
\RequirePackage{txfonts}
% bold math must be loaded after Times font
\usepackage{bm}

%%%%%%%%%%% INCLUDE PACKAGES %%%%%%%%%%%
\RequirePackage{cuted} % defines 'strip' environment for single column
\RequirePackage{indentfirst}	% indent the leading paragraph

%%%%%%%%%%% REQUIRED VARIABLES %%%%%%%%%%%
\newif\if@institutedefined
\@institutedefinedfalse
\def\@institute{\PackageWarning{anstrans}{No institute has been defined.}}
\newcommand{\institute}[1]{\def\@institute{#1}\@institutedefinedtrue}

\newif\if@emaildefined
\@emaildefinedfalse
\def\@email{\PackageWarning{anstrans}{No email address has been defined.}}
\newcommand{\email}[1]{\def\@email{#1}\@emaildefinedtrue}

%%%%%%%%%%% PARAGRAPH FORMATTING %%%%%%%%%%%
\setlength{\parskip}{0pt plus 0pt minus 0pt}
\setlength{\parindent}{0.25in}
%\let\@oldfootnote\footnote
%\renewcommand\footnote[1]{%
%\@oldfootnote{#1}%
%}

%%%%%%%%%%% SECTIONING %%%%%%%%%%%
\setcounter{secnumdepth}{-2} %hide all section heading numbers

\renewcommand\section{\@startsection
	{section}{1}{0pt}%		%name, level, indent
	{-\baselineskip}%		%beforeskip
	{.75\baselineskip}%		%afterskip
	{\reset@font\bf\MakeUppercase}%
}

\renewcommand\subsection{\@startsection
	{subsection}{1}{0pt}%		%name, level, indent
	{-\baselineskip}%		%beforeskip
	{.75\baselineskip}%		%afterskip
	{\reset@font\bf}%
}

\renewcommand\subsubsection{\@startsection
	{subsubsection}{1}{-\parindent}%		%name, level, indent
	{-.75\baselineskip}%		%beforeskip
	{.5\baselineskip}%		%afterskip
	{\reset@font\emph}%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% APPENDIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% appendix needs equations labeled A.1 etc., but shouldn't look like 
%   A   MY FIRST APPENDIX
\let\@oldappendix\appendix
\renewcommand\appendix{%
\setcounter{secnumdepth}{-1}%keep section headings
\renewcommand\theequation{\thesection.\arabic{equation}}
\@oldappendix
\renewcommand\section{%
\addtocounter{section}{1}%
\setcounter{equation}{0}%reset equation counter
\@startsection%
	{section}{1}{0pt}%		%name, level, indent
	{-\baselineskip}%		%beforeskip
	{.75\baselineskip}%		%afterskip
	{\reset@font\bf\MakeUppercase}%
}%
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIBLIOGRAPHY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\@biblabel[1]{#1.} %change from [x] to x.

%% reduce spacing in bibliography; thanks to Jose Carlos Santos.
\let\@oldthebibliography=\thebibliography
\let\@endoldthebibliography=\endthebibliography
\let\@olditemsep=\itemsep
\renewenvironment{thebibliography}[1]{%
	\begin{@oldthebibliography}{#1}%
		\setlength{\itemsep}{-\@olditemsep}%
}%
{%
	\@endoldthebibliography%
	\end{@oldthebibliography}%
}

%%%%%%%%%%% CHANGE TABLE/FIGURE STYLE %%%%%%%%%%%
\renewcommand{\figurename}{Fig.}
\renewcommand{\tablename}{TABLE}

\setlength{\abovecaptionskip}{0.25\baselineskip}
\setlength{\belowcaptionskip}{0.25\baselineskip}
\setlength{\textfloatsep}{0.5\baselineskip}
\setlength{\intextsep}{0.5\baselineskip}

\renewcommand{\bottomfraction}{0.9}
\renewcommand{\topfraction}{0.9}

% captions should be flush with the left and use a stop instead of a colon
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  {#1. #2\par}
  \vskip\belowcaptionskip}

% tables use roman numerals
\renewcommand\thetable{\@Roman\c@table}

%%%% Title section
\setlength{\stripsep}{0pt}
\renewcommand{\maketitle}{%
  \begin{strip}%
    \centering%
      %title
      {\bf \@title}
      \vspace{.75\baselineskip}%

      %author
      {\def\and{\!, } \itshape \@author}%
      \vspace{.75\baselineskip}%

      %institute (optional?)
      \if@institutedefined%
        {\let\and\par \itshape \@institute}% 
      \fi%

      %email
      \if@emaildefined%
        {\def\and{\!, } \itshape \@email}%
      \fi%
      \vspace{\baselineskip}%

  \end{strip}%
}

\AtBeginDocument{%
  \maketitle%
}


