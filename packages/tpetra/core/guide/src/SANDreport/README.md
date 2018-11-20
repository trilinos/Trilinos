
# Steps to build pdf report

- Run `make latex` in the `Source/` directory
- Navigate to the `docs/latex` directory and modify `Tpetra.tex` as follows:
  - If the preceding command did not symbolically link the contents of
    `Source/SANDreport` to `docs/latex`, do so
  - Change document class to `\documentclass[report,10pt]{SANDreport}`
  - Add `\input{FrontMatter.txt}` before `\begin{document}`
  - Remove the date auto populated in the `\date` field
  - Add `cleardoublepage` before `\tableofcontents`
  - Wrap the first level summary with
```
\section*{Summary}
\addcontentsline{toc}{section}{Summary}
```
  - Add `\SANDmain` before the Introduction section
  - Remove the `\subsection{Tpetra License and Copyright Statement}` section
  - Remove the `\subsection{ToDo}` section
  - Remove all `\begin{align}\begin{split}` commands
