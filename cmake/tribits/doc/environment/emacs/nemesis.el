;;---------------------------------------------------------------------------;;
;; file   : nemesis.el
;; author : Thomas M. Evans
;; date   : Wednesday January 27 12:28:3 2010
;; brief  : Defines a useful development environment for Emacs with nemesis
;;          codes
;; note   : Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
;;---------------------------------------------------------------------------;;

(require 'nemesis-functions)
(require 'nemesis-modes)

;;---------------------------------------------------------------------------;;
;; INITIALIZATION VARIABLES
;;---------------------------------------------------------------------------;;

;; Unfortunately, I don't know how to make this work with a good way in Emacs
;; ---once the mode line changes it changes on all buffers, ugh. (with XEmacs
;; it appears to work ok)
(defvar nemesis-colorize-modeline t
  "Use colorized modelines to help show modes while editing.")

;; Do we want doxymacs?
(defvar nemesis-setup-doxymacs nil
  "Have nemesis setup doxymacs.  The doxymacs files (doxymacs.el
  and xml.el) must be findable in the load-path.")

;;---------------------------------------------------------------------------;;
;; GLOBAL KEYS
;;---------------------------------------------------------------------------;;

;; Compiling
(global-set-key (kbd "<f1>") 'compile)
(global-set-key (kbd "<f3>") 'previous-error)
(global-set-key (kbd "<f4>") 'next-error)

;;---------------------------------------------------------------------------;;
;; LANGUAGE AND MODE HOOKS
;;---------------------------------------------------------------------------;;
;; Elisp

(defun nemesis-elisp-mode-hook ()
  "Nemesis hook for elisp editing"
  (local-set-key (kbd "<f5>") 'nemesis-elisp-divider)
  (local-set-key (kbd "<f6>") 'nemesis-elisp-comment-divider)
  (turn-on-auto-fill))

(add-hook 'emacs-lisp-mode-hook 'nemesis-elisp-mode-hook)

;;---------------------------------------------------------------------------;;
;; C/C++

;; Common c-mode hooks
(defun nemesis-c-mode-common-hook ()
  ;; setup indentation
  (nemesis-setup-c-mode)
  (c-set-style "nemesis")
  ;; turn on auto-fill
  (turn-on-auto-fill)
  ;; last column before newline
  (set-fill-column 78)
  ;; ctrl-return is newline-and-indent
  (local-set-key [C-return] 'newline-and-indent)
  ;; find companion file
  (local-set-key (kbd "<f9>") 'nemesis-find-companion-file)
  )

;; C++ specific hooks
(defun nemesis-cc-mode-hook ()
  "Nemesis hook for C++ editing"
  (local-set-key (kbd "<f5>") 'nemesis-cc-divider)
  (local-set-key (kbd "<f6>") 'nemesis-cc-comment-divider))

;; C specific hooks
(defun nemesis-c-mode-hook ()
  "Nemesis hook for C editing"
  (local-set-key (kbd "<f5>") 'nemesis-c-divider)
  (local-set-key (kbd "<f6>") 'nemesis-c-comment-divider))

(add-hook 'c-mode-common-hook 'nemesis-c-mode-common-hook)
(add-hook 'c-mode-hook 'nemesis-c-mode-hook)
(add-hook 'c++-mode-hook 'nemesis-cc-mode-hook)

;; turn on doxymacs if requested
(if nemesis-setup-doxymacs (nemesis-setup-doxymacs-mode))

;; turn on DBC keywords
(defun nemesis-dbc-font-lock-hook ()
  (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))
      (nemesis-font-lock)))
(add-hook 'font-lock-mode-hook 'nemesis-dbc-font-lock-hook)

;;---------------------------------------------------------------------------;;
;; LaTeX/TeX

;; ReFTeX should be built with Emacs
(require 'reftex)
(add-hook 'latex-mode-hook 'turn-on-reftex) ;; Emacs latex-mode
(add-hook 'LaTeX-mode-hook 'turn-on-reftex) ;; AUCTeX LaTeX mode

(defun nemesis-tex-mode-hook ()
  (local-set-key (kbd "<f5>") 'nemesis-latex-divider)
  (local-set-key (kbd "<f6>") 'nemesis-latex-comment-divider)
  (turn-on-auto-fill))

(add-hook 'latex-mode-hook 'nemesis-tex-mode-hook)
(add-hook 'LaTeX-mode-hook 'nemesis-tex-mode-hook)

;;---------------------------------------------------------------------------;;
;; Python

(defun nemesis-python-mode-hook ()
  "Nemesis hook for Python editing"
  (local-set-key (kbd "<f5>") 'nemesis-pound-divider)
  (local-set-key (kbd "<f6>") 'nemesis-pound-comment-divider)
  (local-set-key [C-return] 'newline-and-indent))

(add-hook 'python-mode-hook 'nemesis-python-mode-hook)

;;---------------------------------------------------------------------------;;
;; Shell-Script

(defun nemesis-shell-script-mode-hook ()
  "Nemesis hook for Shell editing"
  (local-set-key (kbd "<f5>") 'nemesis-pound-divider)
  (local-set-key (kbd "<f6>") 'nemesis-pound-comment-divider)
  (local-set-key [C-return] 'newline-and-indent))

(add-hook 'sh-mode-hook 'nemesis-shell-script-mode-hook)

;;---------------------------------------------------------------------------;;
;; FORTRAN

(defun nemesis-f90-mode-hook ()
  "Nemeis hook for FORTRAN editing"
  (local-set-key (kbd "<f5>") 'nemesis-f90-divider)
  (local-set-key (kbd "<f6>") 'nemesis-f90-comment-divider)
  (set-fill-column 80)
  (turn-on-auto-fill))

(setq auto-mode-alist
      (append '(
		("\\.fm4$" . f90-mode)
		) auto-mode-alist))

(add-hook 'f90-mode-hook 'nemesis-f90-mode-hook)

;;---------------------------------------------------------------------------;;

(provide 'nemesis)

;;---------------------------------------------------------------------------;;
;; end of nemesis.el
;;---------------------------------------------------------------------------;;


