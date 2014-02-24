%%---------------------------------------------------------------------------%%
%% TechNote.cls
%%
%% OPTIONS:
%%     memo -> gives memorandum
%%     note -> gives technical note (default)
%%     all other options accepted by article
%% 
%% COMMANDS:
%%     \memotoc -> tableofcontents blocked by single lines
%%     \classification -> puts classification markers at top and
%%                        bottom of page
%%     \revisionnum -> revision number that goes after subject
%%     \caution -> adds ORNL caution box on its own page at the end
%%---------------------------------------------------------------------------%%
%% $Id: ResearchNote.cls,v 1.3 2009/09/03 01:59:52 9te Exp $
%%---------------------------------------------------------------------------%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{TechNote}
[2010/05/14 ORNL Technical Note Class, T.M. Evans]

%% OPTIONS
\DeclareOption{memo}{%
  \def\memotitle{\rule{0pt}{21pt}Memo\kern.02em randum}}
\DeclareOption{note}{%
  \def\memotitle{\rule{0pt}{21pt}Technical \kern.02em Note}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{note}
\ProcessOptions

%% This is based on the article class
\LoadClass{article}
\RequirePackage{graphicx}
\RequirePackage{color}
\RequirePackage[colorlinks]{hyperref}
\RequirePackage{array}

%% Modify the normal text area on the page
\setlength{\textwidth}{6.5in}
\setlength{\textheight}{8.5in}
\setlength{\oddsidemargin}{0in}
\setlength{\topmargin}{-3\baselineskip}
\setlength{\headheight}{5\baselineskip}
\setlength{\headsep}{18pt}
\setlength{\footskip}{3\baselineskip}

\definecolor{dgreen}{rgb}{0.133333, 0.545098, 0.133333}

%% New names
\renewcommand{\abstractname}{Executive Summary}

%% Memo tableofcontents
\newcommand{\lines}{%
  \newline\noindent
  \begin{picture}(0,0)
    \setlength{\unitlength}{1in}
    \put(0,0){\line(1,0){6.5}}
  \end{picture}\newline
}

\newcommand{\doublelines}{%
  \newline\noindent
  \begin{picture}(0,0)
    \setlength{\unitlength}{1in}
    \put(0,0){\line(1,0){6.5}}
    \put(0,.03125){\line(1,0){6.5}}
  \end{picture}\newline
}

\newcommand{\memotoc}{%
  \noindent\lines
  \vspace{-2\baselineskip}
  \tableofcontents
  \noindent\lines
}

\newcommand{\email}[2]{\href{mailto:#2}{#1 $<$#2$>$}}

%% Define a bunch of stuff for the header block
\def\revisionNumber{1}
\def\srd{}
\def\divname{}
\def\grpname{}
\def\fromwho{}
\def\refnum{}
\def\subjectline{}
\def\dist{}

\newcommand*{\revisionnum}[1]{\def\revisionNumber{#1}}
\newcommand*{\classification}[1]{\def\srd{#1}}
\newcommand*{\TIname}[1]{\def\divname{#1}}
\newcommand*{\groupname}[1]{\def\grpname{#1}}
\newcommand*{\from}[1]{\def\fromwho{#1}}
\newcommand*{\refno}[1]{\def\refnum{#1}}
\newcommand*{\subject}[1]{\def\subjectline{#1}}
\newcommand*{\audience}[1]{\def\dist{#1}}

%% headings (regular pages)
\renewcommand{\ps@headings}{%
  \renewcommand{\@oddhead}{
    \parbox[b]{1.5in}{\refnum}
    \hfil
    \parbox[b]{2.5in}{\centerline{\Large\srd}\vspace{3\baselineskip}
      \centerline{-\textsl{\thepage}-}}
    \hfil
    \parbox[b]{1.5in}{\raggedleft\@date}}
  \renewcommand{\@oddfoot}{\centerline{\Large\srd}}%
  \renewcommand{\@evenhead}{}%
  \renewcommand{\@evenfoot}{}%
}
\pagestyle{headings}

%% title page heading
\renewcommand{\ps@empty}{%
  \renewcommand{\@oddhead}{
    \hfil
    \parbox[b]{2.5in}{\centerline{\Large\srd}\vspace{3\baselineskip}
    \centerline{$ $}}
    \hfil}
  \renewcommand{\@oddfoot}{\centerline{\Large\srd}}%
  \renewcommand{\@evenhead}{}%
  \renewcommand{\@evenfoot}{}%
}

%% title page -> opening
\newcommand{\opening}{%
  \thispagestyle{empty}
  \vspace*{-.5in}
  \begin{center}
    \includegraphics[width=6.5in]{WordMarkLeaf}
  \end{center}
  \begin{tabular}{p{2.5in}p{3.5in}}
    \begin{minipage}{2.5in}
      {\sffamily\bfseries\fontsize{21}{21}\selectfont \memotitle} \\
      {\rmfamily\bfseries\slshape\fontsize{8}{9}\selectfont \divname} \\
      {\rmfamily\itshape\fontsize{8}{9}\selectfont \grpname}
    \end{minipage}
    &
    \noindent
    \begin{tabular}{p{0.35in}p{3.15in}}
      {\scriptsize From:} & \fromwho \\
      {\scriptsize Number:} & \refnum \\
      {\scriptsize Date:} & \@date
    \end{tabular}\\\\
    \multicolumn{2}{p{6in}}{
      \textbf{\large Subject: \subjectline\ (Rev. \revisionNumber)}}\\\\
  \end{tabular}
  \doublelines
}

%% end -> closing
\newcommand{\closing}{%
  \section*{Distribution}
  \dist
}

%% caution for publishing
\newcommand{\caution}{%
  \setlength{\parindent}{0in}
  \pagebreak
  \vspace*{\fill}
  \begin{center}
    \fbox{
      \parbox{3in}{\small\vspace{-1\baselineskip}
        \begin{center}\bf CAUTION\end{center}
        This document has not been given final patent clearance and is for
        internal use only. If this document is to be given public release, it
        must be cleared through the site Technical Information Office, which
        will see that the proper patent and technical information reviews are
        completed in accordance with the policies of Oak Ridge National
        Laboratory and UT-Battelle, LLC.}}
  \end{center}
}

%%---------------------------------------------------------------------------%%
%% end of TechNote.cls
%%---------------------------------------------------------------------------%%
