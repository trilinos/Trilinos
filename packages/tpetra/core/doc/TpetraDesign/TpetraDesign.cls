%%-------------------------------------------------------------------------
%% Copyright Notice
%%
%% Copyright (c) 2000, Sandia Corporation, Albuquerque, NM.
%%-------------------------------------------------------------------------

%%-------------------------------------------------------------------------
%% Purpose        : LaTeX class (.cls) file for the Tpetra Design Document report.
%% Special Notes  : Under construction - a merge between the original Xyce
%%                  version hack and the SAND version by Rolf Riesen.
%% Creator        : Michael A. Heroux, Sandia National Labs
%%
%%-------------------------------------------------------------------------
% ---------------------------------------------------------------------------- %
% Started with Scott Hutchinsons Xyce cls file and modified for use with Trilinos.
%
% Michael A. Heroux January 28, 2003.
%
% We start out with the standard LaTeX document class "article" and
% modify it so it conforms to the Sandia Technical report style.
% The modifications are based on "Guide to Preparing SAND Reports",
% SAND98-0730, edited by Tamara K. Locke.
%
% Rolf Riesen, May 29, 2002.
%
%
% 06/14/2002 v1.01
%     Make distribution list wider. Spaces between area codes and phone numbers
% 06/17/2002 v1.02
%     Use \raggedright for main title
% 06/19/2002 v1.03
%     Set tocdepth to 2 to only list sections and subsections in TOC
% 08/07/2002 v1.04
%     Added a Makefile based on one Mike Heroux provided
%     Allow for distribution lists that span more than one page and
%     typeset the distribution in two columns
%


% The version of this file
\typeout{Using Sandia Trilinos Technical Report class:
  January 28, 2003, v1.00 Michael Heroux}
\newcommand\SANDreportVersion{2003/01/28 v1.00}

% We need LaTeX 2e or better
\NeedsTeXFormat{LaTeX2e}

% We need ifthenelse. \newboolean are false by default
\RequirePackage{ifthen}
\RequirePackage{graphicx}
\RequirePackage{multicol}

\newboolean{relaxedSAND}
\newboolean{strictSAND}
\newboolean{SANDmainProvided}
\newboolean{SANDnumProvided}
\newboolean{SANDprintDateProvided}
\newboolean{SANDauthorProvided}
\newboolean{SANDsupersedProvided}
\newboolean{SANDdistributionProvided}
\newboolean{inAppendix}

% We need calc to divide by 3 in setlength
\RequirePackage{calc}

% Needed for Trilinos
\RequirePackage{color}
\RequirePackage{eepic}
\RequirePackage{amssymb}
\RequirePackage{scalefnt}
\RequirePackage{fancyhdr}
\RequirePackage{epsfig}
%\RequirePackage{tabularx}
%\RequirePackage{supertabular}
\RequirePackage{alltt}
%\RequirePackage{colortbl}
%\RequirePackage{framed}
\RequirePackage{verbatim}
\RequirePackage{fancybox}
\RequirePackage{float}


% Specify what this file provides
\ProvidesPackage{TpetraDesign}[\SANDreportVersion]
\ProvidesFile{TpetraDesign.cls}[\SANDreportVersion]


% We understand two options: "relax" and "strict".
% Using "relax" make the document look more LaTeXy
% No option follows the format of SAND98-0730
% Using "strict" tries to follow the guideline even more closely
\DeclareOption{relax}{
  \typeout{"relax": Relaxing some of the SAND98-0730 requirements}
  \setboolean{relaxedSAND}{true}
}
\DeclareOption{strict}{
  \typeout{"strict": Follow the SAND98-0730 requirements more strictly}
  \setboolean{strictSAND}{true}
}
\ifthenelse{\boolean{relaxedSAND} \and \boolean{strictSAND}}   {
  \ClassWarning{SANDreport.cls}{Both, "strict" and "relax" options given}
}{
}

% We're based on the standard article class and assume certain
% options. Set those, and pass all others to the article class
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass[letterpaper,twoside,final,notitlepage]{article}


%
% Setup for letter paper (This IS a SAND report!)
%
\setlength{\textwidth}{6in}
\setlength{\oddsidemargin}{0.5in}
\setlength{\evensidemargin}{0.0in}
\setlength{\marginparwidth}{1.0in - \marginparsep}

% Page header and foooter
\pagestyle{fancyplain}
\addtolength{\headwidth}{\marginparsep}
%\addtolength{\headwidth}{0.5in}
%\addtolength{\headwidth}{\marginparsep}
%\addtolength{\headwidth}{\marginparwidth}
\renewcommand{\headrulewidth}{1.6pt}
\setlength{\headheight}{48pt}
\setlength{\headsep}{24pt}

%\setlength{\textheight}{9in - \footskip - \headheight - \headsep}
\setlength{\textheight}{10in - \footskip - \headheight - \headsep}
\setlength{\topmargin}{-0.5in}

\AtBeginDocument{

  % Do some usage checking
  \ifthenelse{\boolean{SANDnumProvided}}   {
  }{
    \ClassError{TpetraDesign.cls}{\@backslashchar SANDnum
      not provided}{Insert \@backslashchar SANDnum in the preamble
      of your document}
  }
  \ifthenelse{\boolean{SANDauthorProvided}}   {
  }{
    \ClassError{TpetraDesign.cls}{\@backslashchar SANDauthor
      not provided}{Insert \@backslashchar SANDauthor in the preamble
      of your document}
  }
  \ifthenelse{\boolean{SANDprintDateProvided}}   {
  }{
    \ClassError{TpetraDesign.cls}{\@backslashchar SANDprintDate
      not provided}{Insert \@backslashchar SANDprintDate in the
      preamble of your document}
  }

  % Save the title, it will be cleared by \maketitle, but we need it later
  \newsavebox{\SANDtitletempbox}
  \sbox{\SANDtitletempbox}{\parbox{\textwidth - 1.5in}{
      \centering \LARGE \bf \@title \par}
  }

  %
  % Some measurements to format the title page properly
  %
  \newlength{\textwidththird}
  \newlength{\topheight}
  \newlength{\topleftwidth}
  \newlength{\topmiddlewidth}
  \newlength{\toprightwidth}
  \setlength{\textwidththird}{(\textwidth / 3) - (\fboxsep * 2)}
  \setlength{\topleftwidth}{\the\textwidththird}
  \setlength{\topmiddlewidth}{\the\textwidththird}
  \setlength{\toprightwidth}{\the\textwidththird}
  \addtolength{\topleftwidth}{-3em}
  \addtolength{\topmiddlewidth}{+4em}
  \addtolength{\toprightwidth}{-3em}


  \ifthenelse{\boolean{SANDsupersedProvided}}   {
    \setlength{\topheight}{1.5in}
  }{
    \setlength{\topheight}{1.0in}
  }

  %
  % We want a blank line between paragraphs
  %
  \setlength\parskip{2ex}


  %
  % Print the cover pages
  %
  {
    \pagestyle{empty}
    \newlength{\coverindent}
    \setlength{\coverindent}{1em}
    \setlength{\parindent}{0em}

    %
    % Front Cover
    %
    \addtolength{\headheight}{-0.5in}
    \hspace*{-0.5in}
    \begin{minipage}[t]{6.0in}
      {
        \CoverFont{b}{18}{20pt}
        SAND REPORT
        \vspace{0.125in}\\
        \CoverFont{m}{12}{14pt}
        \SANDnumVar\\
        \SANDreleaseTypeVar\\
        \SANDprintDateVar
        \ifthenelse{\boolean{SANDsupersedProvided}}   {
          \vspace{0.125in}\\
          \SANDsupersedVar
        }{
          \\
        }

        \vspace{0.25in}
        \CoverFont{b}{24}{28pt}
        \parbox[t]{\textwidth}{
%          \@title
          \begin{flushleft}
            \Huge\scalefont{1.1}\@title
          \end{flushleft}
%          \vspace{-1em}
%          \Thicklines\raggedleft\line(1,0){400} \\[0.5em]
          \LARGE\bf \SANDsubtitleVar
        }


        \vspace{0.75in}
        \CoverFont{m}{12}{14pt}
%        \SANDauthorVar
        \@author

        \vspace{0.25in}
        \CoverFont{m}{8}{10pt}
        \selectfont
        Prepared by\\
        Sandia National Laboratories\\
        Albuquerque, New Mexico 87185 and Livermore, California 94550\\
        \\
	      Sandia National Laboratories is a multimission laboratory managed and operated\\
        by National Technology and Engineering Solutions of Sandia, LLC, a wholly owned\\
        subsidiary of Honeywell International Inc., for the U.S. Department of Energy’s\\
        National Nuclear Security Administration under contract DE-NA0003525.\\
        \\
        Approved for public release; further dissemination unlimited.\\

        \includegraphics{snllineblk}
%           \includegraphics{snllineblk.eps}
      }
    \end{minipage}

    %
    % Back cover
    %
    \newpage
    \hspace{1.0in}
    \begin{minipage}[t]{5.0in}
      \fontencoding{OT1}
      \fontfamily{pnc}
      \fontseries{m}
      \fontshape{n}
      \fontsize{9}{11pt}
      \selectfont

      Issued by Sandia National Laboratories, operated for the
      United States Department of Energy by Sandia Corporation.
      \vspace{\baselineskip}

      {\fontsize{11}{13pt}\fontseries{b}\selectfont NOTICE:} This report was
      prepared as an account of work sponsored by an agency of the United
      States Government. Neither the United States Government, nor any agency
      thereof, nor any of their employees, nor any of their contractors,
      subcontractors, or their employees, make any warranty, express or
      implied, or assume any legal liability or responsibility for the
      accuracy, completeness, or usefulness of any information, apparatus,
      product, or process disclosed, or represent that its use would not
      infringe privately owned rights. Reference herein to any specific
      commercial product, process, or service by trade name, trademark,
      manufacturer, or otherwise, does not necessarily constitute or imply its
      endorsement, recommendation, or favoring by the United States Government,
      any agency thereof, or any of their contractors or subcontractors. The
      views and opinions expressed herein do not necessarily state or reflect
      those of the United States Government, any agency thereof, or any of
      their contractors.
      \vspace{\baselineskip}

      Printed in the United States of America. This report has been reproduced
      directly from the best available copy.
      \vspace{\baselineskip}

      Available to DOE and DOE contractors from\\
      \begin{tabular}{ll}
        \hspace*{\coverindent} & U.S. Department of Energy \\
        & Office of Scientific and Technical Information \\
        & P.O. Box 62 \\
        & Oak Ridge, TN 37831 \\
        &  \\
      \end{tabular}
      \\
      \begin{tabular}{lll}
        \hspace*{\coverindent} & Telephone: & (865) 576-8401\\
        & Facsimile: & (865) 576-5728\\
        & E-Mail: &         reports@adonis.osti.gov\\
        & Online ordering: & http://www.doe.gov/bridge\\
      \end{tabular}
      \vspace{\baselineskip}
      \vspace{0.25in}

      Available to the public from\\
      \begin{tabular}{ll}
        \hspace*{\coverindent} & U.S. Department of Commerce \\
        & National Technical Information Service \\
        & 5285 Port Royal Rd \\
        & Springfield, VA 22161 \\
        &  \\
      \end{tabular}
      \\
      \begin{tabular}{lll}
        \hspace*{\coverindent} & Telephone: & (800) 553-6847 \\
        & Facsimile: & (703) 605-6900 \\
        & E-Mail: & orders@ntis.fedworld.gov \\
        & Online ordering: & http://www.ntis.gov/ordering.htm \\
      \end{tabular}
      \vspace{\baselineskip}
      \vspace{0.125in}

      \includegraphics[width=0.75in]{DOEbwlogo}
%               \includegraphics[width=0.75in]{DOEbwlogo.eps}
    \end{minipage}
    \newpage
    \addtolength{\headheight}{0.5in}
  }
} % End of AtBeginDocument{}


\AtEndDocument{
  % Do some more usage checking
%   \ifthenelse{\boolean{SANDmainProvided}}   {
%   }{
%     \ClassError{TpetraDesign.cls}{\@backslashchar SANDmain
%       not provided}{Add \@backslashchar SANDmain before your
%       introduction (first) section of your document}
%   }
%  \ifthenelse{\boolean{SANDdistributionProvided}}   {
%  }{
%    \ClassError{TpetraDesign.cls}{\@backslashchar SANDdistribution
%      not provided}{Add \@backslashchar SANDdistribution at the
%      end of your document}
%  }
} % End of AtEndDocument{}


% ---------------------------------------------------------------------------- %
% We need to (re)define some commands from the standard report class.
%

  % ------------------------------------------------------------------------ %
  % Font selections for the cover page \Coverfont{series}{size}{skip}
  %
\newcommand{\CoverFont}[3]{
  \fontencoding{OT1}
  \fontfamily{phv}
  \fontseries{#1}
  \fontshape{n}
  \fontsize{#2}{#3}
  \selectfont
}


    % ------------------------------------------------------------------------ %
    % Here we define the mandatory declarations
    %
\newcommand{\SANDnum}[1]{
  \newsavebox{\SANDnumtempbox}
  \newcommand{\SANDnumVar}{#1}
  \sbox{\SANDnumtempbox}{#1}
  \setboolean{SANDnumProvided}{true}
}

\newcommand{\SANDauthor}[1]{
  \newcommand{\SANDauthorVar}{#1}
  \setboolean{SANDauthorProvided}{true}
}

\newcommand{\SANDprintDate}[1]{
  \newsavebox{\SANDprintDatetempbox}
  \newcommand{\SANDprintDateVar}{#1}
  \sbox{\SANDprintDatetempbox}{#1}
  \setboolean{SANDprintDateProvided}{true}
}


    % ------------------------------------------------------------------------ %
    % Distribution page at the end
    %
\renewcommand{\descriptionlabel}[1]{
  \makebox[3em][r]{#1}
}
\newenvironment{SANDdistribution}{
  \setboolean{SANDdistributionProvided}{true}
  \clearpage
  \setlength{\columnseprule}{1pt}
  \begin{multicols}{2}[{\Large DISTRIBUTION:}]
    \normalsize
    \begin{description}
    }{
    \end{description}
  \end{multicols}
}

\newcommand{\SANDdistInternal}[4]{
  \ifthenelse{\boolean{SANDdistributionProvided}}   {
  }{
    \ClassError{TpetraDesign.cls}{\@backslashchar SANDdistribution
      not (yet) provided}{\@backslashchar SANDdistribution must
      be used before any \@backslashchar SANDdistInternal}
  }
\item[#1] MS \makebox[6ex][r]{#2}\\
  \mbox{}
  \hspace{1em} \parbox[t]{2.0in}{#3, \mbox{#4}}
}


\newcommand{\SANDdistExternal}[2]{
  \ifthenelse{\boolean{SANDdistributionProvided}}   {
  }{
    \ClassError{TpetraDesign.cls}{\@backslashchar SANDdistribution
      not (yet) provided}{\@backslashchar SANDdistribution must
      be used before any \@backslashchar SANDdistExternal}
  }
\item[#1] \parbox[t]{2.1in}{#2}
}


    % ------------------------------------------------------------------------ %
    % And now some optional declarations
    %
\newsavebox{\SANDreleaseTypetempbox}
\newcommand{\SANDreleaseTypeVar}{Unlimited Release}
\sbox{\SANDreleaseTypetempbox}{Unlimited Release}
\newcommand{\SANDreleaseType}[1]{
  \sbox{\SANDreleaseTypetempbox}{#1}
  \renewcommand{\SANDreleaseTypeVar}{#1}
}

\newcommand{\SANDsubtitle}[1]{
  \newcommand{\SANDsubtitleVar}{#1}
}

\newcommand{\SANDsupersed}[2]{
  \newsavebox{\SANDsupersedtempbox}
  \newcommand{\SANDsupersedVar}{
    Supersedes #1 \\
    dated #2
  }
  \sbox{\SANDsupersedtempbox}{
    \begin{tabular}{c}
      Supersedes #1 \\
      dated #2
    \end{tabular}
  }
  \setboolean{SANDsupersedProvided}{true}
}

\newsavebox{\SANDdistcattempbox}
\sbox{\SANDdistcattempbox}{}
\newcommand{\SANDdistcategory}[1]{
  \sbox{\SANDdistcattempbox}{
    \begin{tabular}{r}
      Distribution\\
      Category #1
    \end{tabular}
  }
}


    % Repeat the title and start on an odd page for the main part
 \newcommand{\SANDmain}{
   \cleardoublepage    % Must start on an odd page
   \begin{center}
     \pagestyle{plain}
     \SANDnumVar\\
     \SANDreleaseTypeVar\\
     \SANDprintDateVar
     \ifthenelse{\boolean{SANDsupersedProvided}}   {
       \vspace{0.125in}\\
       \SANDsupersedVar
     }{
       \\
     }
%     \vspace{0.5in}
     \usebox{\SANDtitletempbox} \\[1ex] %
     {\Large\bf \SANDsubtitleVar} \\[0.5in] %
     \SANDauthorVar \\[1ex] %
     Sandia National Laboratories \\ %
     P.O. Box 5800 \\ %
     Albuquerque, NM 87185-0316 %
   \end{center}
   \setboolean{SANDmainProvided}{true}
 }


    % ------------------------------------------------------------------------ %
    % Table of Contents
    % The names we use to label the table of contents, and lists of figures,
    % and tables. Only list Sections and Subsections.
    %
\renewcommand{\contentsname}{\bf Contents}
%\renewcommand{\contentsname}{\raggedright \normalfont\Large\bfseries
%  Contents}
\renewcommand{\listfigurename}{\bf Figures}
\renewcommand{\listtablename}{\bf Tables}
%\renewcommand\refname{\raggedright \bf References}
%\setcounter{tocdepth}{2}

    %
    % Don't enter figures, tables, and subsections of an appendix in
    % the main table of contents.
 \renewcommand{\addcontentsline}[3]{
   \ifthenelse{\boolean{strictSAND}}   {
     \ifthenelse{\boolean{inAppendix}}   {
       \ifthenelse{\equal{#2}{section}}   {
         \addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}
       }{
       }
     }{
       \addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}
     }
   }{
     \addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}
   }
 }

% ------------------------------------------------------------------------ %
% We want the dots closer together, and the same as in \dotfill
%
\renewcommand{\@dotsep}{2.27}

\renewcommand{\l@section}{\@dottedtocline{1}{1.5em}{1em}}


% ------------------------------------------------------------------------ %
% We don't want blank lines between entries and no boldface.
% This is the l@section definition from article.sty. We get rid of
% \addvspace{1.0em... and \bfseries (if relax is not set in the main doc)
%
\renewcommand*\l@section[2]{
  \setlength\parskip{0ex}
  \ifnum \c@tocdepth >\z@
  \addpenalty\@secpenalty
  \ifthenelse{\boolean{relaxedSAND}}   {
    \addvspace{1.0em \@plus\p@}
  }{
  }

% Changed this from 1.5em to 2.5em since if the number of sections went to
% double digits the TOC entry was screwed up - SAH, 20 Sep. 2002.
%  \setlength\@tempdima{1.5em}
  \setlength\@tempdima{2.5em}
  \begingroup
  \parindent \z@ \rightskip \@pnumwidth
  \parfillskip -\@pnumwidth
  \leavevmode
  \ifthenelse{\boolean{relaxedSAND}}   {
    \bfseries
  }{
  }
  \advance\leftskip\@tempdima
  \hskip -\leftskip
  #1\nobreak\dotfill \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
  \endgroup
  \fi}

\renewcommand*\l@figure{\setlength\parskip{0ex}
  \@dottedtocline{1}{1.5em}{2.3em}
}
\let\l@table\l@figure


% ------------------------------------------------------------------------ %
% If there is an appendix, set it apart in the table of contents
%
\renewcommand\appendix{
  \setboolean{inAppendix}{true}
  \par
  %\setcounter{chapter}{0}
  \setcounter{section}{0}
  \setcounter{subsection}{0}
  \setcounter{subsubsection}{0}
  \setcounter{paragraph}{0}
  \setcounter{subparagraph}{0}
  \setcounter{figure}{0}
  \setcounter{table}{0}

  \renewcommand{\thefigure}{\thesection.\arabic{figure}}
  \renewcommand{\thetable}{\thesection.\arabic{table}}

  \gdef\thesection{\@Alph\c@section}
  \addtocontents{toc}{\protect\subsubsection*{\appendixname}}
  \setcounter{secnumdepth}{1}

  % We want each appendix starting on a new page
  % Add \cleardoublepage to the standard definition
  \renewcommand\section{
    \cleardoublepage
    \@startsection {section}{2}{-0.25in}%
%    \@startsection {section}{1}{\z@}%
    {-3.5ex \@plus -1ex \@minus -.2ex}%
    {2.3ex \@plus.2ex}%
    {\normalfont\Huge\color{TrilinosTextBlue}}
  }

  \renewcommand{\sectionmark}{\markboth%
%    {\appendixname~\thesection}%
    {\appendixname~\thesection}%
  }
}

% ------------------------------------------------------------------------ %
% Sections
% We want centered and unnumbered section headings (subsection and
% subsubsection start on the left hand side)
%
\ifthenelse{\boolean{strictSAND}}   {
        % The samples in SAND98-0730 don't number sections.
  \setcounter{secnumdepth}{0}
}{
        % The report SAND98-0730 itself does, though ;-)
  \setcounter{secnumdepth}{2}
}

% Instead of this, use the Trilinos definition below.
%\renewcommand\section{\@startsection {section}{1}{0mm}
%  {-3.5ex \@plus -1ex \@minus -.2ex}
%  {2.3ex \@plus.2ex}
%  {\normalfont\Large\bfseries\centering}}


    % ------------------------------------------------------------------------ %
    % Captions
    % We want the label bold face, and the whole captionwidth about 4",
    % and the text a little smaller than the main text.
    %
\newlength{\SANDcaptionwidth}
\setlength{\SANDcaptionwidth}{4in}
\newsavebox{\SANDtempbox}
\renewcommand{\@makecaption}[2]{
  \par\vspace{\abovecaptionskip}
  \sbox{\SANDtempbox}{{\bf #1. }\ #2}
  \ifthenelse{\lengthtest{\wd\SANDtempbox < \SANDcaptionwidth}}   {
    \begin{center}
      {\bf #1. }\small #2
    \end{center}
  }{
    \begin{center}
      \parbox[t]{\SANDcaptionwidth}{{\bf #1. }\small #2}
    \end{center}
  }
}


    % ------------------------------------------------------------------------ %
    % The basic report class title page is almost what we want. We
    % need to make sure it is on page three and has the SAND number on top
    %
\renewcommand{\@maketitle}{
  \newpage
  {
    % use framebox instead of makebox for debugging
    \makebox[\topleftwidth][t]{
      \begin{minipage}[c][\topheight][t]{\topleftwidth}
                    % There is nothing at top left
      \end{minipage}
    }
    \makebox[\topmiddlewidth][c]{
      \begin{minipage}[c][\topheight][t]{\topmiddlewidth}
        \begin{center}
          \usebox{\SANDnumtempbox}\\
          \usebox{\SANDreleaseTypetempbox}\\
          \usebox{\SANDprintDatetempbox}\\
          \ifthenelse{\boolean{SANDsupersedProvided}}   {
            \bigskip
            \usebox{\SANDsupersedtempbox}\\
          }{
          }
        \end{center}
      \end{minipage}
    }
    \makebox[\toprightwidth][r]{
      \begin{minipage}[c][\topheight][t]{\toprightwidth}
        \raggedleft
        \usebox{\SANDdistcattempbox}
      \end{minipage}
    }
  }

  \begin{center}
    \let \footnote \thanks
    {\LARGE \bf \@title\par%
    \Large\bf \SANDsubtitleVar\par}
%    {\@title\par}
    \vskip 1.5em
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}
%        \@author
        \SANDauthorVar
      \end{tabular}
      \par
    }
    \vskip 1em
    {\large \@date}
  \end{center}
  \par
  \vskip 1.5em
}

%------------------------------------------------------------------------------
% Trilinos added stuff
%------------------------------------------------------------------------------

\setlength{\parindent}{0pt}
\setlength{\parskip}{1ex plus0.5ex minus0.2ex}

% Tabularx environment changes
%\renewcommand{\tabularxcolumn}[1]{>{\vspace{0.5em}}m{#1}}
%\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
%\newcolumntype{Z}{>{\raggedright\arraybackslash}X}
%\newcolumntype{W}{>{\centering\arraybackslash}X}

% Trilinos Colors remapped to BW
\definecolor{TrilinosBlue}{rgb}{0,0,0}
\definecolor{TrilinosDeepBlue}{rgb}{0,0,0}
 \definecolor{TrilinosDarkBlue}{rgb}{0,0,0}
 \definecolor{TrilinosTextBlue}{rgb}{0,0,0}
 \definecolor{TrilinosRed}{rgb}{0,0,0}
 \definecolor{TrilinosDeepRed}{rgb}{0,0,0}
 \definecolor{turquoise}{rgb}{0,0,0}
 \definecolor{grey}{rgb}{0.5,0.5,0.5}
 \definecolor{shadecolor}{rgb}{0.753,0.753,0.753}
 \definecolor{pink}{rgb}{0.5,0.5,0.5}

% Trilinos Colors
% \definecolor{TrilinosBlue}{rgb}{0.498,0.749,1}
% \definecolor{TrilinosDeepBlue}{rgb}{0.4,0.4,0.6}
% \definecolor{TrilinosDarkBlue}{rgb}{0,0,0.5}
% \definecolor{TrilinosTextBlue}{rgb}{0.2,0,0.6}
% \definecolor{TrilinosRed}{rgb}{1,0,0}
% \definecolor{TrilinosDeepRed}{rgb}{0.502,0,0}
% \definecolor{turquoise}{rgb}{0.5,1,1}
% \definecolor{grey}{rgb}{0.5,0.5,0.5}
% \definecolor{shadecolor}{rgb}{0.753,0.753,0.753}
% \definecolor{pink}{rgb}{1,.6,.8}

% Define strings
%\def\Trilinos{{\color{TrilinosBlue}{\bf \textsf{Trilinos}}}
\def\Trilinos{\textsf{\bf Trilinos}}
\def\TrilinosTM{{\color{TrilinosBlue} \bf
      \textsf{Trilinos}}\textsuperscript{\tiny{\texttrademark}}}
\def\TrilinosTitle{\bf
%  {\color{TrilinosBlue}{\Huge\scalefont{1.45}\bf
%  {\color{TrilinosBlue}{\Huge\scalefont{1.00}\bf
  {
    \linethickness{6pt}
    \raggedright{\color{TrilinosDeepBlue}\line(1,0){450}} \\
    \vspace{0.125in}
    \color{TrilinosBlue}
    {
%      \textsf{Trilinos}}}\textsuperscript{\small{\texttrademark}}}
      \textsf{Trilinos}}}\textsuperscript{\small\scalefont{1.1}{\texttrademark}}}

% New commands
\newcommand{\shadebg}[1]{\begin{shaded}{#1}\end{shaded}}
\newcommand{\textrmb}[1]{\textrm{\textbf{#1}}}

% DRAFT
\def\draft{%
\special{!userdict begin /bop-hook{gsave
200 30 translate 65 rotate
/Times-Roman findfont 216 scalefont setfont
0 0 moveto 0.9 setgray (DRAFT) show grestore}def end}
}

% Footer
%\cfoot{\bfseries\thepage \\ DRAFT - Not Approved for External Release}
\cfoot{\bfseries\thepage}

\renewcommand{\familydefault}{phv}
\renewcommand{\sfdefault}{phv}

% Heading
\newcommand{\HeadingC}[1]{%
  \par\vspace{\baselineskip}%
  {\Large\color{TrilinosTextBlue}#1}%
  \par\vspace{\baselineskip}}

% Examples
\newcommand{\Example}[1]{%
  \par\vspace{1ex}%
  {\bf Example:}%
  {\hspace{0.2 in}} {#1} {\par\vspace{1ex}}}

% DisplayCommand
\newcommand{\DisplayCommand}[1]{%
  \par\vspace{1ex}%
  {\bf Command:}%
  {\hspace{0.2 in}} {\tt #1} {\par\vspace{1ex}}}

% InlineCommand
\newcommand{\InlineCommand}[1]{
  {\hspace{0.01 in}} {\tt #1} {\hspace{0.01 in}}}

% InlineCommand
\newcommand{\InlineDirectory}[1]{
  {\hspace{0.01 in}} {\tt #1} {\hspace{0.01 in}}}

% Verbatim-quote environment
\newenvironment{vquote}{\begin{quote}\begin{alltt}}
    {\end{alltt}\end{quote}}

% Trilinos itemize environment
\newenvironment{TrilinosItemize}{%
  \renewcommand{\labelitemi}{{\color{TrilinosRed}$\blacksquare$}}%
  \begin{itemize}}{\end{itemize}}

% Trilinos chapter style
\def\chaptermark#1{}
\def\sec@hangfrom#1{\setbox\@tempboxa\hbox{#1}%
      \hangindent \z@\noindent\box\@tempboxa}

\def\@chapapp{\chaptername}

%\newcommand{\chapterrulecolor}{TrilinosDeepRed}
%\newcommand{\chapterrulewidth}{1.6pt}
%\newcommand{\plainchapterrulewidth}{\z@skip}
%\def\chapterrule{{\if@fancyplain\let\chapterrulewidth\plainchapterrulewidth\fi
%    \color{\chapterrulecolor}\hrule\@height\chapterrulewidth}}

\def\@makechapterhead#1{{\parindent0pt\raggedright
    \hyphenpenalty \@M
    \Huge\scalefont{1.45}\color{TrilinosTextBlue}
    \sec@hangfrom{\thechapter.\hspace{1em}}%
    \ignorespaces#1\par
    \color{TrilinosDeepRed}\vspace{0.5ex}\headrule
    \ifdim\pagetotal>118pt
      \vskip 24pt
    \else
      \@tempdima=118pt\advance\@tempdima by-\pagetotal
      \vskip\@tempdima
    \fi
  }
}

\def\@makeschapterhead#1{{\parindent0pt\raggedright
    \hyphenpenalty \@M
    \Large\bf\boldmath
    \ignorespaces#1\par
    \ifdim\pagetotal>118pt
      \vskip 24pt
    \else
      \@tempdima=118pt\advance\@tempdima by-\pagetotal
      \vskip\@tempdima
    \fi}
}

\def\chapter{\cleardoublepage\thispagestyle{empty}
  \global\@topnum\z@\@afterindentfalse
  \secdef\@chapter\@schapter}

\def\@chapter[#1]#2{%
  \ifnum\c@secnumdepth>\m@ne
    \refstepcounter{chapter}
    \typeout{\@chapapp\space\thechapter}
    \addcontentsline{toc}{chapter}{\protect
      \numberline{\thechapter.}#1}
  \else
    \addcontentsline{toc}{chapter}{#1}
  \fi
  \chaptermark{#1}
  \addtocontents{lof}{\protect\addvspace{10pt}}
  \addtocontents{lot}{\protect\addvspace{10pt}}
  \if@twocolumn
    \@topnewpage[\@makechapterhead{#2}]
  \else \@makechapterhead{#2}
    \@afterheading
  \fi
}

\def\@schapter#1{
  \if@twocolumn\@topnewpage[\@makeschapterhead{#1}]
  \else \@makeschapterhead{#1}
  \@afterheading
  \fi
}

% Chapter Overview style
\newcommand{\chapteroverview}[2]{%
  \vspace{\baselineskip}%
  {\Huge\color{TrilinosTextBlue}#1}%
  \vspace{\baselineskip}\\ #2 \clearpage}

% Trilinos Abstract style
\renewcommand{\abstractname}%
{\raggedright\Large Abstract}
%{\raggedright\Large\color{TrilinosTextBlue} Abstract}

% Trilinos Contents heading style
%\renewcommand{\contentsname}{\raggedright\Huge\color{TrilinosTextBlue} Contents}

% Trilinos List of Figures heading style
%\renewcommand{\listfigurename}{\raggedright\normalfont\Large\bfseries\color{TrilinosTextBlue}
%List of Figures}

% Trilinos section style
\renewcommand{\section}{\@startsection
  {section}
  {1}
  {-0.25in}
  {1ex plus 0.5ex minus 0.5ex}
  {1ex plus 1ex minus 0ex}
  {\raggedright\Huge\color{TrilinosTextBlue}}
%  {\clearpage\raggedright\Huge\color{TrilinosTextBlue}}
}

% Trilinos subsection style
\renewcommand{\subsection}{\@startsection
  {subsection}
  {2}
  {-0.25in}
  {1ex plus 0.5ex minus 0.5ex}
  {0.5ex plus 1ex minus 0ex}
  {\raggedright\LARGE\color{TrilinosTextBlue}}
}

% Trilinos subsubsection style
\renewcommand{\subsubsection}{\@startsection
  {subsubsection}
  {3}
  {-0.25in}
  {0.75ex plus 0.5ex minus 0.5ex}
  {0.5ex plus 0.5ex minus 0ex}
  {\raggedright\large\color{TrilinosTextBlue}}
}

% Trilinos paragraph style
\renewcommand{\paragraph}{\@startsection
  {paragraph}
  {4}
  {\z@}
  {0ex plus 0.1ex minus 0.1ex}
  {0.3ex plus 0.2ex minus 0ex}
  {\raggedright\bfseries\large\color{TrilinosTextBlue}}
}

% Trilinos subparagraph style
\renewcommand{\subparagraph}{\@startsection
  {subparagraph}
  {5}
  {\z@}
  {0ex plus 0.1ex minus 0.1ex}
  {0.3ex plus 0.2ex minus 0ex}
%  {\color{TrilinosTextBlue}\raggedright\normalsize\underline}
  {\raggedright\normalsize\underline}
}

% Trilinos Header style
\renewcommand{\chaptermark}[1]{\markboth{#1}{}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}
\lhead[\fancyplain{}{The Design and Evolution of Tpetra}]%
{\fancyplain{}{\rightmark}}
\rhead[\fancyplain{}{\rightmark}]
{\fancyplain{}{The Design and Evolution of Tpetra}}
