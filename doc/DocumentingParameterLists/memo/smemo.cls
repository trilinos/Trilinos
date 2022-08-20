%
% Version Feb 28, 2008
%
% This is LaTeX class file to create memos that are designed to
% resemble the official Sandia memo format.
%
% I believe the class file was started ~10 years ago in 1520 by Mark
% Blanford. It was rewritten in October of 2004 by a contractor (see
% sfonts.sty for contact information). I (Rich Field) have made a few
% very minor changes since then.
%
\ProvidesClass{smemo}
              [2004/11/12 v1
 Sandia Memo LaTeX Document Class]
\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagetrue
\newif\if@openright
\newif\if@mainmatter \@mainmattertrue
\if@compatibility\else
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\fi
\if@compatibility
  \renewcommand\@ptsize{0}
\else
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\fi
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}
\if@compatibility\else
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}
\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{titlepage}{\@titlepagetrue}
\if@compatibility\else
\DeclareOption{notitlepage}{\@titlepagefalse}
\fi
\if@compatibility
\@openrighttrue
\else
\DeclareOption{openright}{\@openrighttrue}
\DeclareOption{openany}{\@openrightfalse}
\fi
\if@compatibility\else
\DeclareOption{onecolumn}{\@twocolumnfalse}
\fi
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
\ExecuteOptions{letterpaper,10pt,twoside,onecolumn,final,openright}
\ProcessOptions
\input{bk1\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{0\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{2}
\renewcommand\topfraction{.7}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.7}
\renewcommand\dblfloatpagefraction{.5}



  \if@titlepage
  \newcommand\maketitle{\begin{titlepage}%
  \let\footnotesize\small
  \let\footnoterule\relax
  \let \footnote \thanks
  \null\vfil
  \vskip 60\p@
  \begin{center}%
    {\LARGE \@title \par}%
    \vskip 3em%
    {\large
     \lineskip .75em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
      \vskip 1.5em%
    {\large \@date \par}%       % Set date in \large size.
  \end{center}\par
  \@thanks
  \vfil\null
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\else
\newcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks
  \endgroup
  \setcounter{footnote}{0}%
  \global\let\thanks\relax
  \global\let\maketitle\relax
  \global\let\@maketitle\relax
  \global\let\@thanks\@empty
  \global\let\@author\@empty
  \global\let\@date\@empty
  \global\let\@title\@empty
  \global\let\title\relax
  \global\let\author\relax
  \global\let\date\relax
  \global\let\and\relax
}
\def\@maketitle{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
    {\LARGE \@title \par}%
    \vskip 1.5em%
    {\large
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    \vskip 1em%
    {\large \@date}%
  \end{center}%
  \par
  \vskip 1.5em}
\fi
\newcommand*\chaptermark[1]{}
\setcounter{secnumdepth}{2}
\newcounter {part}
\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]
\renewcommand \thepart {\@Roman\c@part}
\renewcommand \thechapter {\@arabic\c@chapter}
\renewcommand \thesection {\thechapter.\@arabic\c@section}
\renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection .\@arabic\c@subsubsection}
\renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
\renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\newcommand\@chapapp{\chaptername}
\newcommand\frontmatter{%
    \cleardoublepage
  \@mainmatterfalse
  \pagenumbering{roman}}
\newcommand\mainmatter{%
    \cleardoublepage
  \@mainmattertrue
  \pagenumbering{arabic}}
\newcommand\backmatter{%
  \if@openright
    \cleardoublepage
  \else
    \clearpage
  \fi
  \@mainmatterfalse}
\newcommand\part{%
  \if@openright
    \cleardoublepage
  \else
    \clearpage
  \fi
  \thispagestyle{plain}%
  \if@twocolumn
    \onecolumn
    \@tempswatrue
  \else
    \@tempswafalse
  \fi
  \null\vfil
  \secdef\@part\@spart}

\def\@part[#1]#2{%
    \ifnum \c@secnumdepth >-2\relax
      \refstepcounter{part}%
      \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
    \else
      \addcontentsline{toc}{part}{#1}%
    \fi
    \markboth{}{}%
    {\centering
     \interlinepenalty \@M
     \normalfont
     \ifnum \c@secnumdepth >-2\relax
       \huge\bfseries \partname~\thepart
       \par
       \vskip 20\p@
     \fi
     \Huge \bfseries #2\par}%
    \@endpart}
\def\@spart#1{%
    {\centering
     \interlinepenalty \@M
     \normalfont
     \Huge \bfseries #1\par}%
    \@endpart}
\def\@endpart{\vfil\newpage
              \if@twoside
               \if@openright
                \null
                \thispagestyle{empty}%
                \newpage
               \fi
              \fi
              \if@tempswa
                \twocolumn
              \fi}
\newcommand\chapter{\if@openright\cleardoublepage\else\clearpage\fi
                    \thispagestyle{plain}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@chapter\@schapter}
\def\@chapter[#1]#2{\ifnum \c@secnumdepth >\m@ne
                       \if@mainmatter
                         \refstepcounter{chapter}%
                         \typeout{\@chapapp\space\thechapter.}%
                         \addcontentsline{toc}{chapter}%
                                   {\protect\numberline{\thechapter}#1}%
                       \else
                         \addcontentsline{toc}{chapter}{#1}%
                       \fi
                    \else
                      \addcontentsline{toc}{chapter}{#1}%
                    \fi
                    \chaptermark{#1}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{#2}]%
                    \else
                      \@makechapterhead{#2}%
                      \@afterheading
                    \fi}
\def\@makechapterhead#1{%
  \vspace*{50\p@}%
  {\parindent \z@ \raggedright \normalfont
    \ifnum \c@secnumdepth >\m@ne
      \if@mainmatter
        \huge\bfseries \@chapapp\space \thechapter
        \par\nobreak
        \vskip 20\p@
      \fi
    \fi
    \interlinepenalty\@M
    \Huge \bfseries #1\par\nobreak
    \vskip 40\p@
  }}
\def\@schapter#1{\if@twocolumn
                   \@topnewpage[\@makeschapterhead{#1}]%
                 \else
                   \@makeschapterhead{#1}%
                   \@afterheading
                 \fi}
\def\@makeschapterhead#1{%
  \vspace*{5\p@}%
  {\parindent \z@ \raggedright
    \normalfont
    \interlinepenalty\@M
    \Large \bfseries  #1\par\nobreak
    \vskip 40\p@
  }}
\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-2.5ex \@plus -1ex \@minus -.2ex}%
                                   {.3ex \@plus.2ex}%
                                   {\normalfont\bfseries}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-2.25ex\@plus -1ex \@minus -.2ex}%
                                     {.5ex \@plus .2ex}%
                                     {\normalfont\bfseries}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\bfseries}}
\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
\newcommand\labelenumi{\theenumi.}
\newcommand\labelenumii{(\theenumii)}
\newcommand\labelenumiii{\theenumiii.}
\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
\newcommand\labelitemi{\textbullet}
\newcommand\labelitemii{\normalfont\bfseries \textendash}
\newcommand\labelitemiii{\textasteriskcentered}
\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}
\newenvironment{verse}
               {\let\\\@centercr
                \list{}{\itemsep      \z@
                        \itemindent   -1.5em%
                        \listparindent\itemindent
                        \rightmargin  \leftmargin
                        \advance\leftmargin 1.5em}%
                \item\relax}
               {\endlist}
\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}
\if@compatibility
\newenvironment{titlepage}
    {%
      \cleardoublepage
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\z@
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
    }
\else
\newenvironment{titlepage}
    {%
      \cleardoublepage
      \if@twocolumn
        \@restonecoltrue\onecolumn
      \else
        \@restonecolfalse\newpage
      \fi
      \thispagestyle{empty}%
      \setcounter{page}\@ne
    }%
    {\if@restonecol\twocolumn \else \newpage \fi
     \if@twoside\else
        \setcounter{page}\@ne
     \fi
    }
\fi
\newcommand\appendix{\par
  \setcounter{chapter}{0}%
  \setcounter{section}{0}%
  \gdef\@chapapp{\appendixname}%
  \gdef\thechapter{\@Alph\c@chapter}}
\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\@addtoreset {equation}{chapter}
\renewcommand\theequation
  {\ifnum \c@chapter>\z@ \thechapter.\fi \@arabic\c@equation}
\newcounter{figure}[chapter]
\renewcommand \thefigure
     {\ifnum \c@chapter>\z@ \thechapter.\fi \@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}[chapter]
\renewcommand \thetable
     {\ifnum \c@chapter>\z@ \thechapter.\fi \@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{2}
\newcommand\tableofcontents{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \chapter*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    \if@restonecol\twocolumn\fi
    }
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty{-\@highpenalty}%
    \addvspace{2.25em \@plus\p@}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
    \endgroup
  \fi}
\newcommand*\l@chapter[2]{%
  \ifnum \c@tocdepth >\m@ne
    \addpenalty{-\@highpenalty}%
    \vskip 1.0em \@plus\p@
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
      \penalty\@highpenalty
    \endgroup
  \fi}
\newcommand*\l@section{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand*\l@subsection{\@dottedtocline{2}{3.8em}{3.2em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{7.0em}{4.1em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{10em}{5em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{12em}{6em}}
\newcommand\listoffigures{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \chapter*{\listfigurename
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}}%
    \@starttoc{lof}%
    \if@restonecol\twocolumn\fi
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \if@twocolumn
      \@restonecoltrue\onecolumn
    \else
      \@restonecolfalse
    \fi
    \chapter*{\listtablename
      \@mkboth{%
          \MakeUppercase\listtablename}{\MakeUppercase\listtablename}}%
    \@starttoc{lot}%
    \if@restonecol\twocolumn\fi
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newenvironment{thebibliography}[1]
     {\vskip48pt{\Large\noindent\bibname
        \@mkboth{\MakeUppercase\bibname}{\MakeUppercase\bibname}}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
\small
      \sloppy
      \clubpenalty4000
      \@clubpenalty \clubpenalty
      \widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \columnseprule \z@
                \columnsep 35\p@
                \twocolumn[\@makeschapterhead{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\@addtoreset{footnote}{chapter}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\newcommand\contentsname{Contents}
\newcommand\listfigurename{List of Figures}
\newcommand\listtablename{List of Tables}
\newcommand\bibname{Bibliography}
\newcommand\indexname{Index}
\newcommand\figurename{Figure}
\newcommand\tablename{Table}
\newcommand\partname{Part}
\newcommand\chaptername{Chapter}
\newcommand\appendixname{Appendix}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
%\pagestyle{headings}
\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Originally from Sandia.tex

%\oddsidemargin 0.3in     %   Left margin on odd-numbered pages.
%\evensidemargin 0.1in    %   Left margin on even-numbered pages.
\marginparwidth .8in     %   Width of marginal notes.

\marginparsep 11pt       %   Horizontal space between outer margin and 
                         %   marginal note

% VERTICAL SPACING:        
                         % Top of page:
\topmargin -0.35in       %    Nominal distance from top of page to top of
                         %    box containing running head.
\headheight 12pt         %    Height of box containing running head.
\headsep 0.3in           %    Space between running head and text.
% \topskip = 10pt        %    '\baselineskip' for first line of page.
                         % Bottom of page:
%\footheight 12pt         %    Height of box containing running foot.
\footskip 0.5in          %    Distance from baseline of box containing foot 
                         %    to baseline of last line of text.

% DIMENSION OF TEXT:

\textheight 8.875in % Height of text (including footnotes and figures, 
                         % excluding running head and foot).
\textwidth 6.25in        % Width of text line.
                         % For two-column mode: 
\columnsep 12pt          %    Space between columns 
\columnseprule 0pt       %    Width of rule between columns.

% FOOTNOTES:

\footnotesep 14pt    % Height of strut placed at the beginning of every
                     % footnote.
\skip\footins 8pt plus 2pt minus 4pt  % Space between last line of text and 
                                      % top of first footnote.

% FLOATS: (a float is something like a figure or table)
%
%  FOR FLOATS ON A TEXT PAGE:
%
%    ONE-COLUMN MODE OR SINGLE-COLUMN FLOATS IN TWO-COLUMN MODE:
\floatsep 14pt plus 2pt minus 4pt        % Space between adjacent floats moved
                                         % to top or bottom of text page.
\textfloatsep 20pt plus 2pt minus 4pt    % Space between main text and floats
                                         % at top or bottom of page.
\intextsep 14pt plus 4pt minus 4pt       % Space between in-text figures and 
                                         % text.

%    TWO-COLUMN FLOATS IN TWO-COLUMN MODE:
\dblfloatsep 14pt plus 2pt minus 4pt     % Same as \floatsep for double-column
                                         % figures in two-column mode.
\dbltextfloatsep 20pt plus 2pt minus 4pt % \textfloatsep for double-column 
                                         % floats.


% MARGINAL NOTES:
%
\marginparpush 7pt       % Minimum vertical separation between two marginal 
                         % notes.


%  ****************************************
%  *           PARAGRAPHING               *
%  ****************************************
%
\parskip 10pt plus2pt minus1pt     % Extra vertical space between paragraphs.
\parindent 2em                     % Width of paragraph indentation.
\topsep 10pt plus 4pt minus 6pt    % Extra vertical space, in addition to 
                                   % \parskip, added above and below list and
                                   % paragraphing environments.
\partopsep 3pt plus 2pt minus 2pt  % Extra vertical space, in addition to 
                                   % \parskip and \topsep, added when user
                                   % leaves blank line before environment.
\itemsep 5pt plus 2.5pt minus 1pt  % Extra vertical space, in addition to
                                   % \parskip, added between list items.

% The following page-breaking penalties are defined

\@lowpenalty   51      % Produced by \nopagebreak[1] or \nolinebreak[1]
\@medpenalty  151      % Produced by \nopagebreak[2] or \nolinebreak[2]
\@highpenalty 301      % Produced by \nopagebreak[3] or \nolinebreak[3]

\@beginparpenalty -\@lowpenalty    % Before a list or paragraph environment.
\@endparpenalty   -\@lowpenalty    % After a list or paragraph environment.
\@itempenalty     -\@lowpenalty    % Between list items.

% \clubpenalty         % 'Club line'  at bottom of page.
% \widowpenalty        % 'Widow line' at top of page.
% \displaywidowpenalty % Math display widow line.
% \predisplaypenalty   % Breaking before a math display.
% \postdisplaypenalty  % Breaking after a math display.
% \interlinepenalty    % Breaking at a line within a paragraph.
% \brokenpenalty       % Breaking after a hyphenated line.


%    ****************************************
%    *                LISTS                 *
%    ****************************************
%

% The following commands are used to set the default values for the list
% environment's parameters. See the LaTeX manual for an explanation of
% the meanings of the parameters.  Defaults for the list environment are
% set as follows.  First, \rightmargin, \listparindent and \itemindent
% are set to 0pt.  Then, for a Kth level list, the command \@listK is
% called, where 'K' denotes 'i', 'ii', ... , 'vi'.  (I.e., \@listiii is
% called for a third-level list.)  By convention, \@listK should set
% \leftmargin to \leftmarginK.
%
% For efficiency, level-one list's values are defined at top level, and
% \@listi is defined to be a no-op.

\leftmargini 2.5em
\leftmarginii 2.2em     % > \labelsep + width of '(m)'
\leftmarginiii 1.87em   % > \labelsep + width of 'vii.'
\leftmarginiv 1.7em     % > \labelsep + width of 'M.'
\leftmarginv 1em
\leftmarginvi 1em

\leftmargin\leftmargini
\labelwidth\leftmargini\advance\labelwidth-\labelsep
\labelsep .5em
\parsep 5pt plus 2.5pt minus 1pt
\let\@listi\relax

\def\@listii{\leftmargin\leftmarginii
   \labelwidth\leftmarginii\advance\labelwidth-\labelsep
   \topsep 5pt plus 2.5pt minus 1pt
   \parsep 2pt plus 1pt minus 1pt
   \itemsep \parsep}

\def\@listiii{\leftmargin\leftmarginiii
    \labelwidth\leftmarginiii\advance\labelwidth-\labelsep
    \topsep 2pt plus 1pt minus 1pt 
    \parsep \z@ \partopsep 1pt plus 0pt minus 1pt
    \itemsep \topsep}

\def\@listiv{\leftmargin\leftmarginiv
     \labelwidth\leftmarginiv\advance\labelwidth-\labelsep}

\def\@listv{\leftmargin\leftmarginv
     \labelwidth\leftmarginv\advance\labelwidth-\labelsep}

\def\@listvi{\leftmargin\leftmarginvi
     \labelwidth\leftmarginvi\advance\labelwidth-\labelsep}

%
% End of size-dependent definitions defaulted to 12pt font

%%%%%%%%%%%%%%%%%
%                                                                               
\mark{{}{}}   % Initializes TeX's marks

%%%%%%%% PAGE STYLES %%%%%%%%%%%%%%%%%%%%

\let\saveshipout\shipout
\long\def\shipout\vbox#1{\saveshipout\vbox to9in{\topclass#1\bottomclass}}

%% Classification hook to add to all page styles:
\def\topclass{\vbox to0pt{\vss\centerline{\@class}\vskip36pt}}
\def\bottomclass{\vbox to0pt{\vskip18pt\centerline{\@class}\vss}}


\def\ps@plain{\let\@mkboth\@gobbletwo
\def\@oddhead{}\let\evenhead\oddhead
\def\@oddfoot{\hbox to\textwidth{\reset@font\hfil-\thepage-\hfil}}
\let\@evenfoot\@oddfoot}


% Definition of 'headings' page style 
%  Note the use of ##1 for parameter of \def\sactionmark inside the
%  \def\ps@headings.
%



\def\ps@myheadings{%
    \let\@oddfoot\@empty\let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil\slshape\leftmark}%
    \def\@oddhead{{\slshape\rightmark}\hfil\thepage}%
    \let\@mkboth\@gobbletwo
    \let\chaptermark\@gobble
    \let\sectionmark\@gobble
    }

\def\ps@headings{%
   \def\@oddhead{%
\vbox{\hbox to\textwidth{\sl\hbox to 0pt{\headname\hss}
      \hfill --\thepage-- \hfill
      \ifdateon\hbox to 0pt{\hss\today}\fi
   } }}%
    \let\@evenhead\@oddhead   % Even page headers same as odd page 
    \def\@oddfoot{}
\let\@evenfoot\@oddfoot
}

\pagestyle{headings}

% Definition of 'myheadings' page style.  The user's \head command
% changes the definition of \@lhead and/or \@rhead.
%
%
% COMMAND: \ps@myheadings
% CHANGES:  The user's \head command became obsolete with LaTeX 2.09.
%           For upward compatibility, it has been reinserted, but now
%           the control sequences \@rhead and \lhead are used unless
%           they are empty definitions. When they are empty, the commands
%           \@rhead and \@lhead will be replaced with \rightmark and
%           \leftmark, respectively.
%
% \head command : Used to set page headings in 'myheadings' style.
%

\def\ps@myheadings{%
  \def\@oddhead{%
    \hbox{}\sl
    \ifx\@rhead\@empty
      \rightmark
    \else
      \@rhead
    \fi
    \hfil \rm\thepage}%
  \def\@evenhead{%
    \rm\thepage\hfil\sl
    \ifx\@lhead\@empty
      \leftmark
    \else
      \@lhead
    \fi
    \hbox{}}%
  \def\@oddfoot{}\def\@evenfoot{}\def\chaptermark##1{}\def\sectionmark##1{}%
}

%%%%%%%%%%%%%%%%%%%%%%%
\def\head{\@ifnextchar[{\@xhead}{\@yhead}\ps@myheadings}
\def\@xhead[#1]#2{\if #1l \def\@lhead{#2}\else \def\@rhead{#2}\fi}
\def\@yhead#1{\def\@lhead{#1}\def\@rhead{#1}}

% Initialization
%
\def\@lhead{}  
\def\@rhead{}

%%%
%    ****************************************
%    *            DISTRIBUTION              *
%    ****************************************
%
%  The distribution environment switches to two column mode if necessary,
%  shrinks up the \baselineskip and \parskip, allowing distribution list
%  to be set in a minimal amount of space.

\def\distribution#1{%
  \@ifnextchar[{\twocol@distrib{#1}}{\onecol@distrib{#1}}%
}

\def\twocol@distrib#1[#2]{%
  \def\@tempa{twocolumn}
  \def\@tempb{#2}
  \ifx\@tempa\@tempb
    \@restonecoltrue
    \if@twocolumn\@restonecolfalse\fi
    \columnseprule \z@
    \twocolumn[{\parindent\z@\raggedright #1\par
              \nobreak\vskip20pt}]%
    \bgroup
      \def\baselinestretch{0.9}
      \small
      \parskip=14pt plus 20pt minus 2pt\relax
      \@distrib
  \else
    \onecol@distrib{#1}%
  \fi
}
\def\onecol@distrib#1{%
  \@restonecolfalse
  {\par\goodbreak\medskip\parindent\z@
   \raggedright
\bf   #1\par
   \nobreak}%
  \bgroup
    \def\baselinestretch{0.9}
    \@normalsize
    \parskip=15pt plus 20pt minus 2pt\relax
    \@distrib
}
\def\@distrib{
    \frenchspacing
    \@rightskip\@flushglue \rightskip\@rightskip
    \raggedbottom
    \parindent\z@
    \interlinepenalty=500
    \pretolerance=5000      % Discourage hyphenation
\distribfont
    \def\attn##1{\parbox[t]{\columnwidth}{\@hangfrom{Attn: \ }\ignorespaces 
                ##1\strut}}
    \def\sandiainternal{\goodbreak\vspace{5pt}%
         {Sandia Internal:}\par
         \interlinepenalty=0
    }%
    \def\externdistrib{\goodbreak\vspace{5pt}\par
         \interlinepenalty=5000
    }%
}
\def\enddistribution{\egroup\if@restonecol\onecolumn\else\vskip 20pt\fi}




%    ************************************
%    *       SANDIA CLASSIFICATION      *
%    ************************************
% CLASSIFICATION MACROS
%
% Classification markings: Classification markings appear at the top and
%                          bottom of the page if TeX is being run from a
%                          secret or confidential account. The markings
%                          are set for the classification level of the
%                          account in the document template. 

\newif\if@classified
\let\classified\relax

\def\classification#1{\global\advance\voffset by .25in
\global\advance\textheight by-.5in \vsize=\textheight
\xdef\@class{{\LARGEtt\uppercase{#1}}}}
\def\@class{}

\def\RSnumber#1{\ifnotempty{#1}{\gdef\RS@num{RS#1}
                \if@draft\def\@evenhead{\rm\today\hfil}
                         \def\@oddhead{\hfil\rm\today}\fi}}
\def\RS@num{}\def\RS@box{\hfil\RS@num}

\def\stopbreaks{\interlinepenalty \@M
 \def\par{\@@par\nobreak}\let\\=\@nobreakcr
 \let\vspace\@nobreakvspace}
\def\@nobreakvspace{\@ifstar{\@nobreakvspacex}{\@nobreakvspacex}}
\def\@nobreakvspacex#1{\ifvmode\nobreak\vskip #1\relax\else
\@bsphack\vadjust{\nobreak\vskip #1}\@esphack\fi}
\def\@nobreakcr{\vadjust{\penalty\@M}\@ifstar{\@xnewline}{\@xnewline}}
\def\startbreaks{\let\\=\@normalcr
\interlinepenalty 0\def\par{\@@par}}


%    ****************************************
%    *             SECTIONS                 *
%    ****************************************
%
% DEFINE COUNTERS:
%
% \newcounter{NEWCTR}[OLDCTR] : Defines NEWCTR to be a counter, which is 
%                               reset to zero when counter OLDCTR is stepped.  
%                               Counter OLDCTR must already be defined.



% For any counter CTR, \theCTR is a macro that defines the printed version
% of counter CTR.  It is defined in terms of the following macros:
%
%  \arabic{COUNTER} : The value of COUNTER printed as an arabic numeral.
%  \roman{COUNTER}  : Its value printed as a lower-case roman numeral.
%  \Roman{COUNTER}  : Its value printed as an upper-case roman numeral.
%  \alph{COUNTER}   : Value of COUNTER printed as a lower-case letter:
%                         1 = a, 2 = b, etc.
%  \Alph{COUNTER}   : Value of COUNTER printed as an upper-case letter:
%                           1 = A, 2 = B, etc.
%

\def\thepart          {\Roman{part}} % Roman numeral part numbers.
\def\thesection       {\arabic{section}}
\def\thesubsection    {\thesection.\arabic{subsection}}
\def\thesubsubsection {\thesubsection.\arabic{subsubsection}}
\def\theparagraph     {\thesubsubsection.\arabic{paragraph}}
\def\thesubparagraph  {\theparagraph.\arabic{subparagraph}}

% \secdef{UNSTARCMDS}{STARCMDS} :
%    When defining a sectionioning command without using
%    \@startsection, you can use \secdef as follows:
%       \def\part { ... \secdef \CMDA \CMDB }
%       \def\CMDA  [#1]#2{ ... }     % Command to define \part[...]{...}
%       \def\CMDB  #1{ ... }         % Command to define \part*{...}

\def\part{\par               % New paragraph
   \addvspace{4ex}           % Adds vertical space above title.
   \@afterindentfalse        % Suppresses indent in first paragraph.  Change
   \secdef\@part\@spart}     % to \@afterindenttrue to have indent.

\def\@part[#1]#2{\ifnum \c@secnumdepth >\m@ne    % IF secnumdepth > -1
        \refstepcounter{part}                    %  THEN step part counter
        \addcontentsline{toc}{part}{\thepart     %       add toc line
        \hspace{1em}#1}\else                     %  ELSE add unnumbered line
      \addcontentsline{toc}{part}{#1}\fi         % FI
   { \parindent 0pt \raggedright 
    \ifnum \c@secnumdepth >\m@ne   % IF secnumdepth > -1 
      \Large \bf Part \thepart     %   THEN Print 'Part' and 
      \par\nobreak                 %          number in \Large boldface.
    \fi                            % FI
    \huge \bf                      % Select \huge boldface.
    #2\markboth{}{}\par }          % Print title and set heading marks null.
    \nobreak                       % TeX penalty to prevent page break.
    \vskip 3ex                     % Space between title and text.
   \@afterheading                  % Routine called after part and
    }                              %     section heading.

% Heading for \part* command
\def\@spart#1{{\parindent 0pt \raggedright 
    \huge \bf                       
    #1\par}                         % Title.
    \nobreak                        % TeX penalty to prevent page break.
    \vskip 3ex                      % Space between title and text.
    \@afterheading                  % Routine called after part and
  }                                 %     section heading.

% \@startsection {NAME}{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE} 
%            optional * [ALTHEADING]{HEADING}
%    Generic command to start a section.  
%    NAME       : e.g., 'subsection'
%    LEVEL      : a number, denoting depth of section -- e.g., section=1,
%                 subsection = 2, etc.  A section number will be printed if
%                 and only if LEVEL < or = the value of the secnumdepth
%                 counter.
%    INDENT     : Indentation of heading from left margin
%    BEFORESKIP : skip to leave above the heading
%    AFTERSKIP  : if positive, then skip to leave below heading,
%                       else - skip to leave to right of run-in heading.
%    STYLE      : commands to set style
%  If '*' missing, then increments the counter.  If it is present, then
%  there should be no [ALTHEADING] argument.  A sectioning command
%  is normally defined to \@startsection + its first six arguments.

\def\section{\@startsection{section}{1}{\z@}{2ex plus .5ex minus 
 .2ex}{3pt}{\large\bf}}
\def\subsection{\@startsection
 {subsection}{2}{\z@}{1.5ex plus .3ex minus .15ex}{1pt}{\normalsize\bf}}
\def\subsubsection{\@startsection
 {subsubsection}{3}{\z@}{1ex plus .2ex minus .1ex}{0.1pt}{\normalsize\bf}}
\def\paragraph{\@startsection
 {paragraph}{4}{\parindent}{\parskip}{-1em}{\normalsize\bf}}
\def\subparagraph{\@startsection
 {subparagraph}{4}{\parindent}{\parskip}{-1em}{\normalsize\bf}}

% Default initializations of \...mark commands.  (See below for their
% us in defining page styles.
%

\def\chaptermark#1{}
% \def\sectionmark#1{}           % Preloaded definitions
% \def\subsectionmark#1{}
% \def\subsubsectionmark#1{}
% \def\paragraphmark#1{}
% \def\subparagraphmark#1{}

% The value of the counter secnumdepth gives the depth of the
% highest-level sectioning command that is to produce section numbers.
%
\setcounter{secnumdepth}{2}

%%%%%%%%%%%%%%%%%%%%%
\setcounter{secnumdepth}{2}

% APPENDIX
%
% The \appendix command must do the following:
%    -- reset the section and subsection counters to zero
%    -- redefine the section counter to produce appendix numbers
%    -- redefine the \section command if appendix titles and headings are
%       to look different from section titles and headings.
\newif\if@sappendix\@sappendixfalse
\long\def\addtocontentsnow#1#2{\if@filesw \begingroup
\def\protect##1{\string\string\string##1\string\space\space}%
\edef\@tempa{\immediate\write \@auxout {\string\@writefile{#1}{#2}}}\@tempa
\endgroup\fi}
\def\appendix{\@ifstar{\@sappendixtrue\@appendix}{\@sappendixfalse\@appendix}}
\def\@appendix{\clearpage
  \setcounter{section}{0}
  \setcounter{subsection}{0}
  \if@chapter
    \setcounter{chapter}{0}
    \def\thechapter{\Alph{chapter}}
    \@addtoreset{figure}{chapter}   % Makes \chapter reset figure counter.
       \def\thefigure{\thechapter.\@arabic\c@figure}
    \@addtoreset{table}{chapter}    % Makes \chapter reset table counter.
       \def\thetable{\thechapter.\@arabic\c@table}
    \@addtoreset{equation}{chapter} % Makes \chapter reset eqn counter.
       \def\theequation{\thechapter.\@arabic\c@equation}
    \def\@chapapp{Appendix}
    \if@sappendix
      \def\@chapter[##1]##2{%
         \refstepcounter{chapter}
         \typeout{\@chapapp.}
         \addcontentsline{toc}{chapter}{\protect
            \numberline{\@chapapp.}##1}%
       \chaptermark{##1}%                              
       \if@twocolumn \@topnewpage[\@makechapterhead{##2}]
       \else \@makechapterhead{##2}\@afterheading\fi}
      \def\@makechapterhead##1{%
       \vspace*{2ex}%
        {\parindent \z@ \centering\Large\bf
         \@chapapp\\ ##1\par
         \nobreak\vskip 2ex
        }%
      }
    \addtocontentsnow{toc}{\protect\setcounter{tocdepth}{0}%
      \def\protect\l@chapter\string\string\string##1\string\string\string##2%
          {\protect\pagebreak[3]\vskip 1.0em plus 1pt
          {\bf\protect\@dottedtocline{0}{0em}{5.2em}%
          {\string\string\string##1}{\string\string\string##2}}}}
    \else
      \def\@chapter[##1]##2{%
         \refstepcounter{chapter}
         \typeout{\@chapapp\space\thechapter.}
         \addcontentsline{toc}{chapter}{\protect
            \numberline{\@chapapp\space\thechapter.}##1}%
       \chaptermark{##1}%                              
       \if@twocolumn \@topnewpage[\@makechapterhead{##2}]
       \else \@makechapterhead{##2}\@afterheading\fi}
      \def\@makechapterhead##1{%
       \vspace*{2ex}%
        {\parindent \z@ \centering\Large\bf
         \@chapapp\space\thechapter\\ ##1\par
         \nobreak\vskip 2ex
        }%
      }
    \addtocontentsnow{toc}{\protect\setcounter{tocdepth}{0}%
      \def\protect\l@chapter\string\string\string##1\string\string\string##2%
          {\protect\pagebreak[3]\vskip 1.0em plus 1pt
          {\bf\protect\@dottedtocline{0}{0em}{6.4em}%
          {\string\string\string##1}{\string\string\string##2}}}}
    \fi
    \addtocontentsnow{lof}{%
       \def\protect\l@figure{\protect\@dottedtocline{1}{0em}{2.8em}}}
    \addtocontentsnow{lot}{%
       \def\protect\l@table{\protect\@dottedtocline{1}{0em}{2.8em}}}
    \def\thebibliography##1{\section*{References}\markboth
      {REFERENCES}{REFERENCES}\list
      {[\thechapter.\arabic{enumi}]}%
      {\settowidth\labelwidth{[##1]}\leftmargin\labelwidth
        \advance\leftmargin\labelsep\usecounter{enumi}}%
        \def\newblock{\hskip .11em plus .33em minus -.07em}
        \sloppy
        \sfcode`\.=1000\relax}
    \def\@bibitem##1{\item\if@filesw \immediate\write\@auxout
           {\string\bibcite{##1}{\thechapter.\the\c@enumi}}\fi\ignorespaces}
  \else
    \def\thesubsection{\Alph{section}.\arabic{subsection}}
    \def\thefigure{\Alph{section}.\@arabic\c@figure}
    \def\thetable{\Alph{section}.\@arabic\c@table}
    \def\theequation{\Alph{section}.\@arabic\c@equation}
    \@addtoreset{figure}{section}   % Makes \section reset figure counter.
    \@addtoreset{table}{section}    % Makes \section reset table counter.
    \@addtoreset{equation}{section} % Makes \section reset eqn counter.
    \if@sappendix
      \def\thesection{Appendix.}
      \addtocontentsnow{toc}{\protect\setcounter{tocdepth}{1}%
          \def\protect\l@section{\protect\@dottedtocline{1}{0em}{5.2em}}}
    \else
      \def\thesection{Appendix \Alph{section}.}
      \addtocontentsnow{toc}{\protect\setcounter{tocdepth}{1}%
           \def\protect\l@section{\protect\@dottedtocline{1}{0em}{6.4em}}}
    \fi
    \addtocontentsnow{lof}{%
       \def\protect\l@figure{\protect\@dottedtocline{2}{0em}{2.8em}}}
    \addtocontentsnow{lot}{%
       \def\protect\l@table{\protect\@dottedtocline{2}{0em}{2.8em}}}
    \def\thebibliography##1{\subsection*{References}\markboth
      {REFERENCES}{REFERENCES}\list
      {[\Alph{section}.\arabic{enumi}]}%
      {\settowidth\labelwidth{[##1]}\leftmargin\labelwidth
        \advance\leftmargin\labelsep\usecounter{enumi}}%
        \def\newblock{\hskip .11em plus .33em minus -.07em}
        \sloppy
        \sfcode`\.=1000\relax}
    \def\@bibitem##1{\item\if@filesw \immediate\write\@auxout
           {\string\bibcite{##1}{\Alph{section}.\the\c@enumi}}\fi\ignorespaces}
  \fi}


%    ****************************************
%    *         OTHER ENVIRONMENTS           *
%    ****************************************
%

% EQUATION and EQNARRAY 
%
% \newcounter{equation}    
%  Default is for left-hand side of equations to be flushleft.
%  To make them flushright, do:
%  \let\@eqnsel = \hfil
%

\def\theequation{\@arabic\c@equation}

% The macro \@eqnnum defines how equation numbers are to appear in equations.
%
% \def\@eqnnum{(\theequation)} 
% 
% For flushleft numbers (in centered equations), use the following:
% \def\@eqnnum{\rlap{(\theequation)}\hskip \displaywidth}

% TITLEPAGE
%  In the normal environments, the titlepage environment does nothing but 
%  end a page and inhibit page numbers.  It also resets the
%  page number to zero.  In two-column style, it still makes a one-column
%  page.  It does not force a page break at the beginning, for compatibility
%  with the \classified macro.
\def\titlepage{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn\fi
     \thispagestyle{empty}\c@page\z@}

\def\endtitlepage{\if@restonecol\twocolumn \else \newpage \fi}

% ARRAY AND TABULAR
%

\arraycolsep 5pt     % Half the space between columns in an array environment.
\tabcolsep 6pt       % Half the space between columns in a tabular environment.
\arrayrulewidth .4pt % Width of rules in array and tabular environment.
\doublerulesep 2pt   % Space between adjacent rules in array or tabular env.

% TABBING
%
\tabbingsep \labelsep   % Space used by the \' command.  (See LaTeX manual.)

% MINIPAGE
%  \@minipagerestore is called upon entry to a minipage environment to
%  set up things that are to be handled differently inside a minipage
%  environment. In the current styles, the only thing it does is change
%  the way footnotes are numbered.  This is done by \let'ing \@mpfnnumber
%  to the '@' version of the appropriate numbering command -- e.g., to
%  have footnotes numbered with upper-case roman numerals, \@mpfnnumber
%  would be \let equal to \@Roman.

\def\@minipagerestore{\let\@mpfnnumber\@alph  % Numbers minipage footnotes 
  }                                           % by \alph .

% FRAMEBOX
%
\fboxsep = 3pt    % Space left between box and text by \fbox and \framebox.
\fboxrule = .4pt  % Width of rules in box made by \fbox and \framebox.


%    ****************************************
%    *             FOOTNOTES                *
%    ****************************************
%
% \footnoterule is a macro to draw the rule separating the footnotes from
% the text.  It should take zero vertical space, so it needs a negative
% skip to compensate for any positive space taken by the rule.  (See
% PLAIN.TEX.)

\def\footnoterule{\kern-1\p@   
  \hrule width .4\columnwidth 
  \kern .6\p@}                 % The \hrule has default height of .4pt .

% The numbering style (arabic, alph, etc.) for ordinary footnotes
% is defined by the macro \@footnotenumber.  It is \let equal
% to the '@' version of the appropriate numbering command -- e.g., to
% have footnotes numbered with upper-case roman numerals, \@footnotenumber
% would be \let equal to \@Roman.
%
\let\@footnotenumber\@arabic     % Defines footnote numbering to be \arabic
%                                                                               
%   \@makefntext{NOTE} :                                                        
%        Must produce the actual footnote, using \@thefnmark as the mark        
%        of the footnote and NOTE as the text.  It is called when effectively   
%        inside a \parbox of width \columnwidth (i.e., with \hsize =            
%        \columnwidth).                                                         
%                                                                               
%        The following macro indents all lines of the footnote by 10pt, and     
%        indents the first line of a new paragraph by 1em.  To change these     
%        dimensions, just substitute the desired value for '10pt' [in both      
%        places] or '1em'.  The mark is flushright against the footnote.        
%          \long\def\@makefntext#1{\@setpar{\@@par\@tempdima \hsize             
%             \advance\@tempdima-10pt\parshape \@ne 10pt \@tempdima}\par        
%             \parindent 1em\noindent \hbox to \z@{\hss$^{\@thefnmark}$}#1}     
%                                                                               
%        A simpler macro is used, in which the footnote text is                 
%        set like an ordinary text paragraph, with no indentation except        
%        on the first line of a paragraph, and the first line of the            
%        footnote.  Thus, all the macro must do is set \parindent               
%        to the appropriate value for succeeding paragraphs and put the         
%        proper indentation before mark.                                        
%
\long\def\@makefntext#1{\noindent % Macro to make the text of a footnote
    \hangindent 0.5em\hbox to 0.5em{\hss $^{\@thefnmark}$}#1}


%    ****************************************
%    *         FIGURES AND TABLES           *
%    ****************************************
% 
% Float placement parameters.  See LaTeX manual for their definition.
%
\setcounter{topnumber}{2}
\def\topfraction{.75}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.75}
\setcounter{totalnumber}{3}
\def\textfraction{0.0}
\def\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.77}
\def\dblfloatpagefraction{.5}

% \@makecaption{NUMBER}{TEXT} : Macro to make a figure or table caption.  
%      NUMBER : Figure or table number--e.g., 'Figure 3.2' 
%      TEXT   : The caption text.
%  Macro should be called inside a \parbox of right width, with \normalsize.

% Modified by G.D.Sjaardema to produce hanging indentation e.g. 
% Figure X.x:  This is the .....
%              caption for this figure
%
\newbox\@tempboxa
\newbox\@tempboxb
\long\def\@makecaption#1#2{%
      {\edef\next{\write16{#1}}\next}%
   \def\@tempa{figure}\ifx\@tempa\@captype\vskip2ex\fi
   \setbox\@tempboxa=\hbox{\bf#1.\space\space}%
   \setbox\@tempboxb=\hbox{#2}%
   \@tempdima\hsize
   \advance\@tempdima by -\wd\@tempboxa
   \ifdim \wd\@tempboxb > \@tempdima   % IF longer than one line:
%   THEN set as ordinary paragraph.
       \hangindent 1\wd\@tempboxa\noindent\box\@tempboxa
       \raggedright\hyphenpenalty=10000\unhbox\@tempboxb\par      
     \else                        
%   ELSE  center.
       \hbox to\hsize{\hfil\box\@tempboxa\box\@tempboxb\hfil}%
   \fi
   \def\@tempa{table}\ifx\@tempa\@captype\vskip2ex\fi}

% To define a float of type TYPE (e.g., TYPE = figure), the document style 
% must define the following.
%
%  \fps@TYPE   : The default placement specifier for floats of type TYPE.
%
%  \ftype@TYPE : The type number for floats of type TYPE.  Each TYPE has
%                associated a unique positive TYPE NUMBER, which is a power 
%                of two.  E.g., figures might have type number 1, tables type 
%                number 2, programs type number 4, etc.
%
%  \ext@TYPE   : The file extension indicating the file on which the 
%                contents list for float type TYPE is stored.  For example, 
%                \ext@figure = 'lof'.
%
%  \fnum@TYPE  : A macro to generate the figure number for a caption.
%                For example, \fnum@TYPE == Figure \thefigure.
%
%  The actual float-making environment commands--e.g., the commands
%  \figure and \endfigure--are defined in terms of the macros \@float
%  and \end@float, which are described below.
%
%  \@float{TYPE}[PLACEMENT] : Macro to begin a float environment for a 
%     single-column float of type TYPE with PLACEMENT as the placement 
%     specifier.  The default value of PLACEMENT is defined by \fps@TYPE.   
%     The environment is ended by \end@float.  
%     E.g., \figure == \@float{figure}, \endfigure == \end@float.


% FIGURE


\def\thefigure{\@arabic\c@figure}

\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{Figure \thefigure}
\def\figure{\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}

% TABLE
%

\def\thetable{\@arabic\c@table}

\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{Table \thetable}
\def\table{\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@float{table}}
\@namedef{endtable*}{\end@float}


%    ****************************************
%    *                LISTS                 *
%    ****************************************
%
%
% ENUMERATE
%  In the following 'N' refers to 'i', 'ii', 'iii', 'iv', etc. (i.e.
%  \theenumN refers to \theenumi, \theenumii, \theenumiii, \theenumiv, etc.)
%
%  \labelenumN ... determines how the number is printed in the label.
%  \theenumN ... determines how the number is printed on each individual level.
%  \p@enumN\theenumN ... determines how the number is printed by a \ref command.

\let\@enumspacing\relax
%
\def\theenumi{\arabic{enumi}} 
\def\theenumii{\alph{enumii}}
\def\theenumiii{\roman{enumiii}}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumii{\theenumi}
\def\p@enumiii{\theenumi(\theenumii)}
\def\p@enumiv{\p@enumiii\thenumiii}

%
\def\labelenumi{\theenumi.} 
\def\labelenumii{(\theenumii)}
\def\labelenumiii{\theenumiii.}
\def\labelenumiv{\theenumiv.}
%
% ITEMIZE
%  \item@... determines the labels of the various itemization levels.
%
\let\@itemspacing=\relax
\def\labelitemi{$\bullet$}
\def\labelitemii{\bf --}
\def\labelitemiii{$\ast$}
\def\labelitemiv{$\cdot$}


% VERSE
%   The verse environment is defined by making clever use of the
%   list environment's parameters.  The user types \\ to end a line.
%   This is implemented by \let'in \\ equal \@centercr.
%
\def\verse{\let\\=\@centercr  
  \list{}{\itemsep\z@ \itemindent -1.5em\listparindent \itemindent 
          \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]}
\let\endverse\endlist

% QUOTATION
%   Fills lines
%   Indents paragraph
%   
\def\quotation{\list{}{\listparindent 1.5em
    \itemindent\listparindent
    \rightmargin\leftmargin}\item[]}
\let\endquotation=\endlist

% QUOTE -- same as quotation except no paragraph indentation,
%
\def\quote{\list{}{\rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

% DESCRIPTION 
%
%  The description environment produces a list with all but the labeled
%  line indented by a fixed amount. The text following the label is indented
%  from the end of the label by a fixed amount equal to \labelsep, regardless
%  of the length of the label.
%
%  The description environment under LaTeX 2.06 used a dimension 
%  \descriptionmargin to determine the minimum width for the label and
%  indented all lines after the labelled line by this amount as well
%  (i.e. \leftmargin set to \descriptionmargin and \labelwidth equal to
%  \leftmargin - \labelsep (which made the left edge of the box containing
%  the label line up with the prevailing left margin of the outer environment).
%
%  To provide compatibility with older documents that may have depended on
%  the setting of \descriptionmargin, the description environment has been
%  modified to first check \descriptionmargin. If \descriptionmargin > 0pt,
%  the description environment uses \descriptionmargin to determine the
%  minimum width of the label and the subsequent indentation just as the
%  description environment in LaTeX 2.06 behaved. If \descriptionmargin is
%  not greater than 0pt, the description environment behaves as it should
%  under LaTeX 2.09 (see description environment in LaTeX User's Guide).
%
%  The \descriptionmargin dimension is initialized to 0pt, thereby making
%  the description environment behave as described in the LaTeX User's
%  Guide. Whenever \descriptionmargin is set to something larger than 0pt,
%  the older definition (LaTeX 2.06) is used.
%
%  The description environment uses \descriptionlabel to produce the label
%  in a description. To change the formatting of the label, you must
%  redefine \descriptionlabel.  For example,
%      \def\descriptionlabel #1{\it #1:\hfil}
%  will cause the label to be set in italic, flush with the left margin,
%  followed by an (italic) colon.  Remember that \descriptionlabel must contain
%  an \hfil in the older definition.

\def\descriptionlabel#1{\hfil\bf #1}

\def\description{%
 \ifdim\descriptionmargin>0pt
   \list{}{\leftmargin\descriptionmargin \labelwidth\descriptionmargin
     \advance\leftmargin\labelsep
     \itemindent\z@
     \let\makelabel\descriptionlabel}%
 \else
   \list{}{\labelwidth\z@ \advance\leftmargin-\labelsep\itemindent-\leftmargin
     \advance\leftmargin\labelsep\let\makelabel\descriptionlabel}%
 \fi
}
\let\enddescription\endlist

\newdimen\descriptionmargin
\descriptionmargin=0pt


%
%    ****************************************
%    *          REFERENCE FORMAT            *
%    ****************************************
%
% COMMAND: \refnum
% ACTION:  Sets its argument as a superscript.
%
\gdef\refnum#1{$^{#1}$}


%    ****************************************
%    *         TABLE OF CONTENTS, ETC.      *
%    ****************************************
%
% A \subsection command writes a 
%       \contentsline{subsection}{TITLE}{PAGE}
% command on the .toc file, where TITLE contains the contents of the
% entry and PAGE is the page number.  If subsections are being numbered,
% then TITLE will be of the form
%       \numberline{NUM}{HEADING}
% where NUM is the number produced by \thesubsection.  Other sectioning
% commands work similarly.  
%
% A \caption command in a 'figure' environment writes
%    \contentsline{figure}{\numberline{NUM}{CAPTION}}{PAGE}
% on the .lof file, where NUM is the number produced by \thefigure and
% CAPTION is the figure caption.  It works similarly for a 'table' environment.
%
% The command \contentsline{NAME} expands to \l@NAME.  So, to specify
% the table of contents, we must define \l@chapter, \l@section, 
% \l@subsection, ... ; to specify the list of figures, we must define
% \l@figure; and so on.  Most of these can be defined with the
% \@dottedtocline command, which works as follows.
%
% \@dottedtocline{LEVEL}{INDENT}{NUMWIDTH}{TITLE}{PAGE}
%    LEVEL    : An entry is produced only if LEVEL < or = value of 
%               'tocdepth' counter.  Note, \chapter is level 0, \section
%               is level 1, etc.
%    INDENT   : The indentation from the outer left margin of the start of 
%               the contents line.
%    NUMWIDTH : The width of a box in which the section number is to go,
%               if TITLE includes a \numberline command.
%   
% This command uses the following three parameters, which are set
% with a \def (so em's can be used to make them depend upon the font).
%   \@pnumwidth : The width of a box in which the page number is put.
%   \@tocrmarg  : The right margin for multiple line entries.  One
%                 wants \@tocrmarg > or = \@pnumwidth
%   \@dotsep    : Separation between dots, in mu units.  Should be \def'd to
%                 a number like 2 or 1.7  (18 mu = 1 em)
%
% G.D.Sjaardema: \@pnumwidth and \@tocrmarg should be integer multiples 
%                 of \@dotsep so that dots line up correctly for two line
%                 captions.  Mess with \@tocrmarg to get perfect lineup.

\def\@pnumwidth{1.5em}
\def\@tocrmarg {3.02em}
\def\@dotsep{9}
\setcounter{tocdepth}{2}

% TABLEOFCONTENTS
%  In report style, \tableofcontents, \listoffigures, etc. are always
%  set in single-column style.  @restonecol

\def\tableofcontents{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
  \fi\section*{Contents}\markboth{CONTENTS}{CONTENTS}
  {\parskip3pt plus2pt minus 1pt\@starttoc{toc}}\if@restonecol\twocolumn\fi}

\def\l@section{\@dottedtocline{1}{0em}{1.4em}}
\def\l@subsection{\@dottedtocline{2}{1.4em}{2.3em}}
\def\l@subsubsection{\@dottedtocline{3}{3.7em}{3.2em}}
\def\l@paragraph{\@dottedtocline{4}{6.9em}{4.1em}}
\def\l@subparagraph{\@dottedtocline{5}{11em}{5em}}

% LIST OF FIGURES
%
\def\listoffigures{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
  \fi\vfil\pagebreak[3]\vfilneg\section*{Figures}\markboth{FIGURES}{FIGURES}
  {\parskip3pt plus2pt minus 1pt\@starttoc{lof}}\if@restonecol\twocolumn\fi}

\def\l@figure{\@dottedtocline{0}{0em}{2.8em}}

% TABLES
%
\def\listoftables{\@restonecolfalse\if@twocolumn\@restonecoltrue\onecolumn
  \fi\vfil\pagebreak[3]\vfilneg\section*{Tables}\markboth{TABLES}{TABLES}
  {\parskip3pt plus2pt minus 1pt\@starttoc{lot}}\if@restonecol\twocolumn\fi}

\let\l@table\l@figure


%    ****************************************
%    *             BIBLIOGRAPHY             *
%    ****************************************
%
\def\thebibliography#1{\section*{References}\markboth
  {REFERENCES}{REFERENCES}\list
  {[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}\leftmargin\labelwidth
    \advance\leftmargin\labelsep\usecounter{enumi}}%
    \def\newblock{\hskip .11em plus .33em minus -.07em}
    \sloppy
    \sfcode`\.=1000\relax
    \addcontentsline{toc}{section}{References}}

\let\endthebibliography=\endlist

\def\@up#1{\leavevmode\raise .16ex\hbox{#1}}
\let\up\@up
\def\@cite#1#2{\@up[{#1\if@tempswa , #2\fi}\@up]}
\def\@biblabel#1{\@up[#1\@up]\hfill}

%    ****************************************
%    *              THE INDEX               *
%    ****************************************
%
% THE THEINDEX ENVIRONMENT
% Produces double column format, with each paragraph a separate entry.
% The user commands \item, \subitem and \subsubitem are used to
% produce the entries, and \indexspace adds an extra vertical space
% that's the right size to put above the first entry with a new letter
% of the alphabet.

\def\theindex{\@restonecoltrue\if@twocolumn\@restonecolfalse\fi
\columnseprule \z@
\columnsep 35pt\twocolumn[\section*{Index}]
    \addcontentsline{toc}{section}{Index}
    \markboth{INDEX}{INDEX}\thispagestyle{plain}\parindent\z@
    \parskip\z@ plus .3pt\relax\let\item\@idxitem}

\def\@idxitem{\par\hangindent 40pt}
\def\subitem{\par\hangindent 40pt \hspace*{20pt}}
\def\subsubitem{\par\hangindent 40pt \hspace*{30pt}}
\def\endtheindex{\if@restonecol\onecolumn\else\clearpage\fi}
\def\indexspace{\par \vskip 10pt plus 5pt minus 3pt\relax}


%    ****************************************
%    *            MISCELLANEOUS             *
%    ****************************************
%
% DATE
%
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
%
\def\Month{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\year}
\newif\ifdateon
\dateontrue   % Default
\def\nodate{\dateonfalse}
\def\date#1{\ifnotempty{#1}{\gdef\today{#1}\dateontrue}}
%
% Redefinition of \label to print labels in the left margin if user has 
% specified \viewlabeltrue
%
\newif\ifviewlabel
\viewlabelfalse  % Default
\def\label#1{\@bsphack\if@filesw {\let\thepage\relax
   \xdef\@gtempa{\write\@auxout{\string
      \newlabel{#1}{{\@currentlabel}{\thepage}}}}}\@gtempa
   \ifviewlabel\def\@strutdepth{\dp\strutbox}%
      \ifvmode{\ifinner\baselineskip\z@\else\baselineskip\@strutdepth\fi
               \label@label{#1}}
      \else\vadjust{\label@label{#1}}\fi
   \fi
   \if@nobreak \ifvmode\nobreak\fi\fi\fi\@esphack}
\def\label@label#1{\kern-\@strutdepth
  \vtop to \@strutdepth{\baselineskip\@strutdepth\vss
  \llap{\rlap{\tentt #1}\hskip 1in}\null}}
%
% This is a different ordering of footnote symbols than standard LaTeX.
% Currently commented out until someone comes up with a good reason to use it.
%
%\def\@fnsymbol#1{\ifcase#1\or* \or{**} \or\dag \or\ddag \or\S \or\| 
%               \or{\dag\dag} \or{\ddag\ddag} \or{\S\S} \else\P \fi} 
%
% Sandia National Labs defined since typed often in References:
%
\def\SANDIA{Sandia National Laboratories, Albuquerque, NM}
%
% "BOLD" Greek --  p386 TeXbook - prints argument three times, moving it 
%                                 slightly each time.
\def\bg#1{\setbox0=\hbox{$#1$}%         
      \kern-.015em\copy0\kern-\wd0
      \kern .03em\copy0\kern-\wd0
      \kern-.015em\raise.0333em\box0}
%
% \makeqnum sets ? to be an empty space the width of a roman ``0''
% this can be used in tables where the numbers have different lengths
% to line them up by the decimal point.      See TeXbook, p. 241 top

\newdimen\@digitwd
\setbox0=\hbox{\rm0}
\@digitwd=\wd0
{\catcode`?=\active\gdef?{\kern\@digitwd\relax}}
\def\makeqnum{\catcode`?=13\relax}
\def\makeqquest{\catcode`?=12\relax}
%
% Use ``\eref{reference}'' instead of ``\ref{reference}'' to get the 
%     reference in parentheses e.g. Equation (4.3.1)
\def\eref#1{(\ref{#1})}
% Use: \ifnotempty{ARG}{ACTION}
% if ARG is null or a space or tab, do nothing; else do ACTION
\long\def\ifnotempty#1#2{\def\@tempa{#1}\ifx\@tempa\@empty\else
  \ifx\@tempa\space\else #2\fi\fi}

%    ****************************************
%    *          INITIALIZATIONS             *
%    ****************************************
%
%  Initialization of skip amounts
%
\smallskipamount=.5\parskip \medskipamount=\parskip
\bigskipamount=2\parskip


\pagenumbering{arabic}      % Arabic page numbers
\onecolumn                  % Single-column.


%% End from Sandia.tex

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\typeout{Document Class 'smemo'. Version 1.0 <1 October 2004>}
%
% SMEMO DOCUMENT Class
%
% **** this is the style file for Sandia memo format *******
% See the example file 'SMEMO.tex' for the use of this     *
% style file. Written by Amy Hendrickson, TeXnology Inc,   *
% based on 1990 .sty file  Written by Mark Blanford, 1521  *
%***********************************************************

%\@options

%endinput

\newif\ifdateon
\global\dateontrue
\def\nodate{\global\dateonfalse}

\raggedbottom

%    ****************************************
%    *             MEMO MACROS              *
%    ****************************************
%


\parindent 0pt
\def\headtext#1{\def\one{#1}\ifx\one\empty\else
{\gdef\headname{#1}}\fi}\def\headname{}
\def\to#1{\begin{list}{}{\leftmargin0pt\labelwidth5em\topsep10pt\partopsep3pt%
 \itemsep5pt\parsep5pt}
\ifdateon
 \item[\smallss date:] {\ifx\@date\empty \today\else\@date \fi}
\fi
 \item[\smallss to:] #1
 \headtext{#1}}
\def\from#1{\vglue .3in
\item[\smallss from:] #1}
\long\def\subject#1\par{\vglue .3in\item[\smallss subject:] #1 \end{list}}
\def\initials#1{\ifnotempty{#1}{\par\@hangfrom{\rm Initials:} \rm\it #1\par\goodbreak}}



\def\encl#1{\ifnotempty{#1}{\par\noindent
  \parbox[t]{\textwidth}{\@hangfrom{\rm enc: }\ignorespaces #1\strut}\par}}
\def\enc#1{\ifnotempty{#1}{\par\noindent
  \parbox[t]{\textwidth}{\@hangfrom{\rm enc: }\ignorespaces #1\strut}\par}}
\def\keywords#1{\ifnotempty{#1}{\par\noindent
  \parbox[t]{\textwidth}{\@hangfrom{Key Words: }\ignorespaces\bf #1\strut}\par}}

%******************************************************************
% MEMO environment prints the Sandia letter head and defines      *
% the subsequent page style with running heads. Note, there is no *
% bothsides option for the heads                                  *
%******************************************************************
\newdimen\movetoright
\movetoright=.35in
\def\memo{\newpage\c@page\@ne\interlinepenalty=100%
\vspace*{-20pt}%
  \ifx\@empty\fromaddress
  \hbox to \textwidth{\hskip-30pt\lower 4pt%
        \vbox{\baselineskip15pt\halign{\hfil##\cr\medss\lhname\cr
              \noalign{\vskip-2pt}\smallss\lhtitle\cr}}\hfil
        \vbox{\baselineskip15pt\halign{\hfil##\hfil\cr
% modified by RF on 2/26/08 to use logo with blue chicken
%\hbox{\includegraphics[width=3in]{SandiaNL}}\cr
\hbox{\includegraphics[width=3in]{snllineblubrd}}\cr
\noalign{\vskip-10pt}
\hskip\movetoright
\tinylogofont Operated for the U.S. Department of Energy by\cr
\noalign{\vskip-2pt}
\hskip\movetoright\smlogofont Sandia Corporation\cr
\noalign{\vskip2pt}
\hskip\movetoright\tinylogofont  Albuquerque, New Mexico 87185\cr}}}
  \else  % home address
    {\raggedleft
      \begin{tabular}{l}\ignorespaces
\hskip\movetoright        \fromaddress
      \end{tabular}%
      \par
    }%
  \fi
\begin{figure}[b]
\vskip24pt
\hbox to \textwidth{\hss\small\textsl{Exceptional Service in the National Interest}\hss}
\vskip-24pt
\end{figure}
\thispagestyle{empty}
}
\def\address#1{\def\fromaddress{#1}}
\def\fromaddress{}
\def\lhname{}\def\lhtitle{}

\def\endmemo{\@@par\pagebreak\@@par}




\renewcommand \thesection {\@arabic\c@section}


\oddsidemargin=0pt
\evensidemargin=0pt
\advance\textwidth by .25in

\let\savedocument\document
\def\document{\definefonts\savedocument}

\def\definefonts{
\expandafter\ifx\csname pscurrentversion\endcsname\relax
%% if true this means that \usepackage{sfonts}
%% has not been done, so we need to define these font
%% names using Computer Roman fonts:
%% Times-Roman
\def\timesroman{cmr10 }
\def\timesbold{cmbx10 }
\def\timesitalic{cmti10 }
\def\timesbolditalic{cmbx10 }
%% Helvetica
\def\helvetica{cmss10 }
\def\helveticaoblique{cmssi10 }
\def\helveticaboldoblique{cmssbx10 }
\def\helveticabold{cmssbx10 }
\def\helveticanarrow{cmss10 }
\def\helveticanarrowbold{cmssbx10 }
\typeout{^^J^^J Using CM Fonts^^J^^J}
\else
\typeout{^^J^^J Using PS Fonts^^J^^J}
\fi

\ifnum\@ptsize=2 %% If 12pt is used as an option, these fonts are 2pts larger:
\global\font\smallss=\helveticaoblique at 9pt
\global\font\medss=\helvetica at 12pt
\global\font\largess=\helvetica at 13pt
\global\font\distribfont=\helvetica at 10pt
\else %% Default is 10pt
\global\font\medss=\helvetica at 10pt
\global\font\smallss=\helveticaoblique at 7pt
\global\font\largess=\helvetica at 11pt
\global\font\distribfont=\helvetica at 8pt
\fi

%% These font sizes stay the same, used for Sandia logo
\global\font\logofont=\helveticanarrowbold at 14pt
\global\font\smlogofont=\helveticabold at 9pt
\global\font\tinylogofont=\helvetica at 7pt

}

%% for classification markings
\global\font\LARGEtt=cmtt10 at 24pt


\endinput
%%
