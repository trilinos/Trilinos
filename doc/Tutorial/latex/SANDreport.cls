% ---------------------------------------------------------------------------- %
% We start out with the standard LaTeX document class "article" and
% modify it so it conforms to the Sandia Technical report style.
% The modifications are based on "Guide to Preparing SAND Reports",
% SAND98-0730, edited by Tamara K. Locke, and the newer "Guide to
% Preparing SAND Reports and Other Communication Products", SAND2002-2068P.
%
% Get the latest version and more at
%    http://www.cs.sandia.gov/~rolf/SANDreport
%
% Rolf Riesen, May 29, 2002.
%
%
% 06/14/2002 v1.01
%     Make distribution list wider. Spaces between area codes and phone numbers
%
% 06/17/2002 v1.02
%     Use \raggedright for main title
%
% 06/19/2002 v1.03
%     Set tocdepth to 2 to only list sections and subsections in TOC
%
% 08/07/2002 v1.04
%     Added a Makefile based on one Mike Heroux provided
%     Allow for distribution lists that span more than one page and
%     typeset the distribution in two columns
%
% 11/25/2002 v1.05
%     Some changes from Daniel Segalman to conform to the new guide
%     (SAND2002-2068P):
%     - Adjust top margin of the title page to conform to the MS Word
%         template (this may or may not be correct depending on your
%         latex/dvips initial setup. Use texconfig to setup paper and
%         dvips).
%     - Added the word "Printed" before the date on front page
%     - Classification category is no longer needed on page 3; get
%         rid of \SANDdistcategory
%     - Added "National Nuclear Security Administration" on title page
%     - Changed NTIS URL on page 2
%
%     If "strict" is not used, we'll number sections to two decimal
%     places (1.2.3)
%
% 02/13/2003 v1.06
%     New DOEbwlogo.eps from Tim Draelos. This one is smaller and has some
%     incompatible junk removed.
%
% 10/16/2003 v1.7
%     Using \ifthenelse in the addcontentsline command interefered with
%     the subfigure class. James Gruetzner and Todd Pitts provided a fix.
%
%     Philippe Pebay reports that the color package needs to be importet
%     before the graphicx package. This class file uses the graphicx package,
%     so Philippe could not use the color package because graphicx was already
%     loaded when his \usepackage command got executed.


    % The version of this file
    \typeout{Using Sandia Technical Report class: October 16, 2003, v1.07 Rolf Riesen}
    \newcommand\SANDreportVersion{2002/10/16 v1.07}

    % We need LaTeX 2e or better
    \NeedsTeXFormat{LaTeX2e}

    % We need ifthenelse. \newboolean are false by default
    \RequirePackage{ifthen}
    \RequirePackage{color}
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


    % Specify what this file provides
    \ProvidesPackage{SANDreport}[\SANDreportVersion]
    \ProvidesFile{SANDreport.cls}[\SANDreportVersion]


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
    \LoadClass[letterpaper,twoside,final]{article}


    %
    % Setup for letter paper (This IS a SAND report!)
    %
    \setlength{\textwidth}{6.2in} %% was different
    \setlength{\oddsidemargin}{0.0in} %% was 0.5
    \setlength{\evensidemargin}{0.0in}
    \setlength{\marginparwidth}{1.0in - \marginparsep}

    \setlength{\textheight}{8.9in - \footskip - \headheight - \headsep} 
    \setlength{\topmargin}{0.2in}
 


    \AtBeginDocument{

        % Do some usage checking
        \ifthenelse{\boolean{SANDnumProvided}}   {
        }{
            \ClassError{SANDreport.cls}{\@backslashchar SANDnum
                not provided}{Insert \@backslashchar SANDnum in the preamble
                of your document}
        }
        \ifthenelse{\boolean{SANDauthorProvided}}   {
        }{
            \ClassError{SANDreport.cls}{\@backslashchar SANDauthor
                not provided}{Insert \@backslashchar SANDauthor in the preamble
                of your document}
        }
        \ifthenelse{\boolean{SANDprintDateProvided}}   {
        }{
            \ClassError{SANDreport.cls}{\@backslashchar SANDprintDate
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
%%%        \setlength\parskip{2ex}


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
            \addtolength{\headheight}{-0.625in}
            \hspace*{-0.5in}
            \begin{minipage}[t]{6.0in}
            {
            \CoverFont{b}{18}{20pt}
            SAND REPORT
            \vspace{0.125in}\\
            \CoverFont{m}{12}{14pt}
            \SANDnumVar\\
            \SANDreleaseTypeVar\\
            Printed \SANDprintDateVar
            \ifthenelse{\boolean{SANDsupersedProvided}}   {
                \vspace{0.125in}\\
                \SANDsupersedVar
            }{
                \\
            }

            \vspace{0.25in}
            \CoverFont{b}{24}{28pt}
            \parbox[t]{\textwidth}{
                \raggedright
                \@title            
            }


            \vspace{0.75in}
            \CoverFont{m}{12}{14pt}
            \SANDauthorVar

            \vspace{0.25in}
            \CoverFont{m}{8}{10pt}
            \selectfont
            Prepared by\\
            Sandia National Laboratories\\
            Albuquerque, New Mexico 87185 and Livermore, California 94550\\
            \\
            Sandia is a multiprogram laboratory operated by Sandia Corporation,\\
            a Lockheed Martin Company, for the United States Department of Energy's\\
            National Nuclear Security Administration under Contract DE-AC04-94-AL85000.\\
            \\
            Approved for public release; further dissemination unlimited.\\

            \includegraphics{snllineblk.eps}
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

                {\fontsize{11}{13pt}\fontseries{b}\selectfont NOTICE:}
                This report was prepared
                as an account of work sponsored by an agency of
                the United States Government. Neither the United
                States Government, nor any agency thereof, nor any
                of their employees, nor any of their contractors,
                subcontractors, or their employees, make any warranty,
                express or implied, or assume any legal liability or
                responsibility for the accuracy, completeness, or
                usefulness of any information, apparatus, product,
                or process disclosed, or represent that its use
                would not infringe privately owned rights. Reference
                herein to any specific commercial product, process,
                or service by trade name, trademark, manufacturer, or
                otherwise, does not necessarily constitute or imply
                its endorsement, recommendation, or favoring by the
                United States Government, any agency thereof, or any
                of their contractors or subcontractors. The views and
                opinions expressed herein do not necessarily state
                or reflect those of the United States Government,
                any agency thereof, or any of their contractors.
                \vspace{\baselineskip}

                Printed in the United States of America. This report has
                been reproduced directly from the best available copy.
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
                    & Online ordering: & http://www.ntis.gov/help/ordermethods.asp?loc=7-4-0\#online \\
                \end{tabular}
                \vspace{\baselineskip}
                \vspace{0.125in}

                \includegraphics[width=0.75in]{DOEbwlogo.eps}
            \end{minipage}
            \newpage
            \addtolength{\headheight}{0.5in}
        }
    } % End of AtBeginDocument{}



    \AtEndDocument{
        % Do some more usage checking
        \ifthenelse{\boolean{SANDmainProvided}}   {
        }{
            \ClassError{SANDreport.cls}{\@backslashchar SANDmain
                not provided}{Add \@backslashchar SANDmain before your
                introduction (first) section of your document}
        }
        \ifthenelse{\boolean{SANDdistributionProvided}}   {
        }{
%           \ClassError{SANDreport.cls}{\@backslashchar SANDdistribution
%               not provided}{Add \@backslashchar SANDdistribution at the
%               end of your document}
        }
    } % End of AtEndDocument{}



% ---------------------------------------------------------------------------- %
% We need to (re)define some commands from the standard article class.
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
            \ClassError{SANDreport.cls}{\@backslashchar SANDdistribution
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
            \ClassError{SANDreport.cls}{\@backslashchar SANDdistribution
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


    % Distribution category is no longer needed
    \newcommand{\SANDdistcategory}[1]{
        \ClassWarning{SANDreport.cls}{Distribution category is no longer needed!}
    }


    % Repeat the title and start on an odd page for the main part
    \newcommand{\SANDmain}{
        \cleardoublepage    % Must start on an odd page
        \begin{center}
            \usebox{\SANDtitletempbox}
        \end{center}
        \setboolean{SANDmainProvided}{true}
    }


    % ------------------------------------------------------------------------ %
    % Table of Contents
    % The names we use to label the table of contents, and lists of figures,
    % and tables. Only list Sections and Subsections.
    %
    \renewcommand{\contentsname}{\raggedright \normalfont\Large\bfseries
        Contents}
    \renewcommand{\listfigurename}{\raggedright \bf Figures}
    \renewcommand\listtablename{\raggedright \bf Tables}
    \renewcommand\refname{\raggedright \bf References}
    \setcounter{tocdepth}{2}

    %
    % Don't enter figures, tables, and subsections of an appendix in
    % the main table of contents.
    % Code from James Gruetzner and Todd Pitts. Avoids the ifthen package in
    % \addcontentsline which interferes with the subfigure package.
    \renewcommand{\addcontentsline}[3]{
        \ifstrictSAND%
            \ifinAppendix%
                \if{#2}{section}
                    \addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}
                \else
                    \relax
                \fi
            \else%
                \addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}
            \fi
        \else%
            \addtocontents{#1}{\protect\contentsline{#2}{#3}{\thepage}}
        \fi
    }


    % ------------------------------------------------------------------------ %
    % We want the dots closer together, and the same as in \dotfill
    %
    \renewcommand{\@dotsep}{2.27}


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
        \setlength\@tempdima{1.5em}
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
        \setcounter{section}{0}
        \setcounter{subsection}{0}
        \setcounter{subsubsection}{0}
        \setcounter{figure}{0}
        \setcounter{table}{0}

        \renewcommand{\thefigure}{\thesection.\arabic{figure}}
        \renewcommand{\thetable}{\thesection.\arabic{table}}

        \gdef\thesection{\@Alph\c@section}
        \addtocontents{toc}{\protect\subsubsection*{\appendixname}}
        \ifthenelse{\boolean{strictSAND}}   {
            \setcounter{secnumdepth}{1}
        }{
            \setcounter{secnumdepth}{3}
        }

        % We want each appendix starting on a new page
        % Add \cleardoublepage to the standard definition
        \renewcommand\section{
            \cleardoublepage
            \@startsection {section}{1}{\z@}%
            {-3.5ex \@plus -1ex \@minus -.2ex}%
            {2.3ex \@plus.2ex}%
            {\normalfont\Large\bfseries}
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
        \setcounter{secnumdepth}{3}
    }
    \renewcommand\section{\@startsection {section}{1}{0mm}
        {-3.5ex \@plus -1ex \@minus -.2ex}
        {2.3ex \@plus.2ex}
        {\normalfont\Large\bfseries\centering}}


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
    % The basic article class title page is almost what we want. We
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
                    Printed \usebox{\SANDprintDatetempbox}\\
                    \ifthenelse{\boolean{SANDsupersedProvided}}   {
                        \bigskip
                        \usebox{\SANDsupersedtempbox}\\
                    }{
                    }
                \end{center}
                \end{minipage}
            }
        }

        \begin{center}
            \let \footnote \thanks
            {\LARGE \bf \@title\par}
            \vskip 1.5em
            {\large
                \lineskip .5em%
                \begin{tabular}[t]{c}
                \@author
                \end{tabular}\par
            }
        \end{center}
        \par
        \vskip 1.5em
    }
