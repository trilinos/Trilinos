;;; Kevin Long's (krlong) C++ programming style. This is built on the "gnu"
;;; style (where I copied out the file cc-styles.el first) and then added
;;; Kevin's style overrides.
(defconst krlong-c-style
  '((c-basic-offset . 2)
    (c-comment-only-line-offset . (0 . 0))
    (c-offsets-alist . (
      (statement-block-intro . +)
      (knr-argdecl-intro . 5)
      (substatement-open . +)
      (label . 0)
      (statement-case-open . +)
      (statement-cont . +)
      (arglist-intro . c-lineup-arglist-intro-after-paren)
      (arglist-close . c-lineup-arglist)
      (inline-open . 0)
      (brace-list-open . +)
      ))
    (c-special-indent-hook . c-gnu-impose-minimum)
    (c-block-comment-prefix . "")
    (c-basic-indent 2)
    )
  "Kevin Long's C++ Programming Style"
  )

(defun krlong-c-mode-common-hook ()
  (c-add-style "krlong" krlong-c-style)
  (setq fill-column 78)
  (setq indent-tabs-mode nil)
  (setq tab-width 2)
  )

(add-hook 'c-mode-common-hook 'krlong-c-mode-common-hook)

;;; Default Thyra C++ programming style.  This was built on the "stroustrup"
;;; format copied for the standard file cc-styles.el and then other options
;;; where added to customize the style.  This style is designed to match the
;;; recommendations in the book "Code Complete", 2nd edition, in Chapter 31,
;;; "Layout and Style".  There are two deviations from what is recommended in
;;; "Code Complete", 2nd edition.  First, (arglist-close . +) should be
;;; (arglist-close . 0) to be consistent.  However, I feel that this looks bad
;;; so I have indented this one offset.  This could be changed however,
;;; depending on what people think about this.  Second, The braces
;;; 
(defconst thyra-c-style
  '((c-basic-offset . 2)
    (c-comment-only-line-offset . 0)
    (c-offsets-alist . (
       (innamespace . 0)  ;;; Don't indent for namespace enclusures
       (inextern-lang . 0)  ;;; As for namespaces, don't indent for extern "C" {
       (substatement-open . 0) ;;; Don't indent opening '{' from control statement
       (statement-block-intro . +) ;;; Indent code in "{' one offset
       (label . 0) ;;; Don;'t indent labels to show them off
       (statement-cont . +) ;;; Indent statement continuation lines one offset
       (case-label . +) ;;; Indent case labels one offset from opening switch
       (arglist-intro . +) ;;; Indent first line of argument list one offset
       (arglist-cont . 0) ;;; Indent arguments after blah(\n one offset
       (arglist-cont-nonempty . +) ;;; Indent args after blah( int arg1 on next line one offset
       (arglist-close . +) ;;; Align the closing ')' with arguments
       ))
    (c-label-minimum-indentation 0) ;;; Was in "stroustrup" stype
    (c-basic-indent 2) ;;; Indent 2 spaces by default
    )
  "Default C++ coding style for Thyra"
  )

;;; NOTE: I had to remove the option (c-tab-always-indent t) since it was
;;; causing string behavior in putting in tabs inside of comments.
;;; Instead, I have added this to the main emacs style ???

(defun thyra-c-mode-common-hook ()
  (c-add-style "thyra" thyra-c-style t) ;;; 't' == make default
  (setq fill-column 78)
  (setq indent-tabs-mode nil)
  (setq tab-width 2)
  (setq c-auto-align-backslashes nil)
  )

;; Note that above, we make 'thyra' the default style since we can not use
;; '(c-default-style "thyra") since only the default styles in cc-styles.el
;; are recongnised by c-default-style!
;;
;; THEREFORE: The Thyra style must be put last in order for it to be used by
;; default.

(add-hook 'c-mode-common-hook 'thyra-c-mode-common-hook)

;; Provide these styles a lisp package
(provide 'cpp-thyra-styles)
