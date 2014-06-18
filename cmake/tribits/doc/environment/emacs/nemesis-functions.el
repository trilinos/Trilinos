;;-------------------------------*-emacs-list-*------------------------------;;
;; file   : nemesis-functions.el
;; author : Thomas M. Evans
;; date   : <2010-01-25 13:30:56 9te>
;; brief  : Defines useful functions that are used in the nemesis
;;          development environment
;; note   : Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;
;; Useful functions
;;---------------------------------------------------------------------------;;
;; Save and kill a buffer

(defun nemesis-save-and-kill-buffer ()
  "Save the buffer then kill it without any annoying queries."
  (interactive)
  (if (buffer-file-name (current-buffer))
      (save-buffer))
  (kill-buffer (buffer-name)))

;;---------------------------------------------------------------------------;;
;; Insert a time-stamp here

(defun nemesis-time-stamp ()
  "Inserts a time stamp."
  (interactive)
  (require 'time-stamp)
  (setq nts (time-stamp-string "%:a %:b %:d %:H:%:M:%:S %:y"))
  (insert nts))

;;---------------------------------------------------------------------------;;
;; Toggle buffers

(defun nemesis-toggle-previous-buffer ()
  "Toggle to a previous buffer"
  (interactive)
  (switch-to-buffer (other-buffer (buffer-name))))

;;---------------------------------------------------------------------------;;
;; Elisp comments
;;---------------------------------------------------------------------------;;

(defun nemesis-elisp-divider ()
"Insert an elisp comment block."
  (interactive)
  (beginning-of-line)
  (insert ";;---------------------------------------------------------------------------;;\n")
  (insert ";; \n")
  (insert ";;---------------------------------------------------------------------------;;\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun nemesis-elisp-comment-divider ()
"Insert an elisp comment block."
  (interactive)
  (beginning-of-line)
  (insert ";;---------------------------------------------------------------------------;;\n")
  (end-of-line)
)

;;---------------------------------------------------------------------------;;
;; C/C++ comments
;;---------------------------------------------------------------------------;;

(defun nemesis-cc-divider ()
"Insert an C++ comment block."
  (interactive)
  (beginning-of-line)
  (insert "//---------------------------------------------------------------------------//\n")
  (insert "// \n")
  (insert "//---------------------------------------------------------------------------//\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun nemesis-cc-comment-divider ()
"Insert an C++ comment block."
  (interactive)
  (beginning-of-line)
  (insert "//---------------------------------------------------------------------------//\n")
  (end-of-line)
)

(defun nemesis-c-divider ()
"Insert an C comment block."
  (interactive)
  (beginning-of-line)
  (insert "/*---------------------------------------------------------------------------*/\n")
  (insert "/* \n")
  (insert "/*---------------------------------------------------------------------------*/\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun nemesis-c-comment-divider ()
"Insert an C comment block."
  (interactive)
  (beginning-of-line)
  (insert "/*---------------------------------------------------------------------------*/\n")
  (end-of-line)
)

;;---------------------------------------------------------------------------;;
;; F90 comments
;;---------------------------------------------------------------------------;;

(defun nemesis-f90-divider ()
"Insert an f90 comment block."
  (interactive)
  (beginning-of-line)
  (insert "!-----------------------------------------------------------------------------!\n")
  (insert "! \n")
  (insert "!-----------------------------------------------------------------------------!\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun nemesis-f90-comment-divider ()
"Insert an f90 comment block."
  (interactive)
  (beginning-of-line)
  (insert "!-----------------------------------------------------------------------------!\n")
  (end-of-line)
)

;;---------------------------------------------------------------------------;;
;; # Comments tools
;;---------------------------------------------------------------------------;;

(defun nemesis-pound-divider ()
"Insert a pound-symbol comment block."
  (interactive)
  (beginning-of-line)
  (insert "##---------------------------------------------------------------------------##\n")
  (insert "## \n")
  (insert "##---------------------------------------------------------------------------##\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun nemesis-pound-comment-divider ()
"Insert a pound-symbol comment block."
  (interactive)
  (beginning-of-line)
  (insert "##---------------------------------------------------------------------------##\n")
  (end-of-line)
)

;;---------------------------------------------------------------------------;;
;; LaTeX/TeX
;;---------------------------------------------------------------------------;;

(defun nemesis-latex-divider ()
"Insert a LaTeX comment block."
  (interactive)
  (beginning-of-line)
  (insert "%%---------------------------------------------------------------------------%%\n")
  (insert "%% \n")
  (insert "%%---------------------------------------------------------------------------%%\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun nemesis-latex-comment-divider ()
"Insert a LaTeX divider."
  (interactive)
  (beginning-of-line)
  (insert "%%---------------------------------------------------------------------------%%\n")
  (end-of-line)
)

;;---------------------------------------------------------------------------;;

(provide 'nemesis-functions)

;;---------------------------------------------------------------------------;;
;; end of nemesis-functions.el
;;---------------------------------------------------------------------------;;
