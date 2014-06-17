;;---------------------------------------------------------------------------;;
;; file   : nemesis-modes.el
;; author : Thomas M. Evans
;; date   : Wednesday January 27 12:27:58 2010 
;; brief  : Mode modifications and additions for nemesis Emacs environment
;; note   : Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;
;; SETUP C/C++ MODE 
;;---------------------------------------------------------------------------;;

(defun nemesis-setup-c-mode()
  "Accepted Nemesis indentation style for C/C++. Added to c-mode-common-hook
  when loading nemesis environment."
  
  (c-add-style
   "nemesis" '
   (
    ; basic offset
    (c-basic-offset . 4)
    ; # goes to first column
    (c-electric-pound-behavior . 'alignleft)
    ; indentation settings
    (c-offsets-alist . ((access-label . -2)
			(block-close . 0)
			(block-open . 0)
			(case-label  . +)
			(class-close . 0)
			(class-open  . 0)
			(defun-close . 0)
			(defun-open  . 0)
			(do-while-closure  . 0)
			(else-clause       . 0)
			(extern-lang-close . +)
			(extern-lang-open  . +)
			(inline-close      . 0)
			(inline-open       . 0)
			(innamespace       . 0)
			(statement-case-intro . +)
			(statement-cont    . c-lineup-math)
			(substatement-open . 0))))))

;;---------------------------------------------------------------------------;;
;; companion files

;; C/C++ companion files
(setq nemesis-c-companion-alist
      '(
        ;; Header pairings
        ("\\(.h\\)$" . ".c")
        ("\\(.hh\\)$" . ".cc")
        ("\\(.hh\\)$" . ".t.hh")
        ("\\(.hh\\)$" . ".i.hh")
        ;; Implementation pairings
        ("\\(.c\\)$" . ".h")
        ("\\(.cc\\)$" . ".hh")
        ("\\(.t.hh\\)$" . ".hh")
        ("\\(.i.hh\\)$" . ".hh")))

;; Find a companion file
(defun nemesis-find-companion-file ()
  "Function to locate the corresponding .hh .i.hh, .t.hh or .cc file.
   When a .hh file is in the current buffer and this function is
   run, the corresponding .cc file will be loaded if it is
   available.  If it is not available, the script will look for a
   corresponding .i.hh file to load.

   The mapping between file types is stored in the emacs variable
   nemesis-c-companion-alist."
  
  (interactive)
  (let ((companion nemesis-c-companion-alist)
	(pair-file ""))
    (catch 'found
      (while companion
	(if (string-match (car (car companion)) buffer-file-name)
	    (progn
	      ;; Found a match, now check to see if its the right one.
	      (setq pair-file (replace-match (cdr (car companion))
					     t t buffer-file-name))
	      (if (file-exists-p pair-file)
		  (progn
		    (message "found matching file, throwing 'found")
		    (throw 'found t))
		(setq companion (cdr companion))))
	  (message (concat "discarding companion=" pair-file))
	  (setq companion (cdr companion)))))
    (if companion
	(find-file pair-file))))

;;---------------------------------------------------------------------------;;
;; CONFIGURE FOR DOXYMACS
;;---------------------------------------------------------------------------;;

(defun nemesis-setup-doxymacs-mode ()
  "Setup doxymacs."
  (interactive)
  (progn
    (require 'doxymacs)
    (defvar doxymacs-doxygen-style "Qt")
    (add-hook 'c-mode-common-hook 'doxymacs-mode)
    (defun nemesis-doxymacs-font-lock-hook ()
      (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))
	  (doxymacs-font-lock)))
    (add-hook 'font-lock-mode-hook 'nemesis-doxymacs-font-lock-hook)))

;;---------------------------------------------------------------------------;;
;; DBC KEYWORDS
;;---------------------------------------------------------------------------;;
;; This is heavily borrowed from doxymacs.el

(defconst nemesis-dbc-keywords
  (list
   (list
    ;; NEMESIS DBC - Match single keyword that is followed by 0 or more
    ;; spaces, followed by an opening paren.
    "\\<\\(Require\\|Ensure\\|Check\\|Remember\\|Insist\\|Assert\\)\\>\\([ ]*\\s(\\)"
    '(0 font-lock-keyword-face prepend))
    ))

(defun nemesis-font-lock ()
  "Turn on font-lock for Nemesis keywords."
  ;; FIXME How do I turn *off* font-lock for Doxygen keywords?
  (interactive)
  (if (functionp 'font-lock-add-keywords)
      ;; Use new (proper?) font-lock-add-keywords function
      (font-lock-add-keywords nil nemesis-dbc-keywords)
    ;; Use old-school way
    (let ((old (if (eq (car-safe font-lock-keywords) t)
		 (cdr font-lock-keywords)
	       font-lock-keywords)))
      (setq font-lock-keywords (append old nemesis-dbc-keywords)))))

;;---------------------------------------------------------------------------;;

(provide 'nemesis-modes)

;;---------------------------------------------------------------------------;;
;; end of nemesis-modes.el
;;---------------------------------------------------------------------------;;

