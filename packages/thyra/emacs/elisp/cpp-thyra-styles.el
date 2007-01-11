;;; Kevin Long's (krlong) C++ programming style This is built on the "gnu"
;;; style (where I copied out the file cc-styles.el first) and then added
;;; Kevin's style overrides.
(defconst krlong-c-style
  '((c-basic-offset . 2)
    (c-comment-only-line-offset . (0 . 0))
    (c-offsets-alist . ((statement-block-intro . +)
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

;;; Default Thyra C++ programming style.  This was built
;;; on the "stroustrup" format copied for the standard file cc-styles.el and
;;; then other options where added to customize the style.
(defconst thyra-c-style
  '((c-basic-offset . 2)
    (c-comment-only-line-offset . 0)
    (c-offsets-alist . ((statement-block-intro . +)
                        (substatement-open . 0)
                        (label . 0)
                        (statement-cont . +)
                        (case-label . +)
                        (inextern-lang . 0)
                        (innamespace . 0)
                        ))
    (c-tab-always-indent t)
    (c-label-minimum-indentation 0)
    (c-basic-indent 2)
    )
  "Default C++ coding style for Thyra"
  )

(defun thyra-c-mode-common-hook ()
  (c-add-style "thyra" thyra-c-style t) ;;; 't' == make default
  (setq fill-column 78)
  (setq indent-tabs-mode nil)
  (setq tab-width 2)
  )

;; Note that above, we make 'thyra' the default style since we can not use
;; '(c-default-style "thyra") since only the default styles in cc-styles.el
;; are recongnised by c-default-style!
;;
;; THEREFORE: The Thyra style must be put last in order for it to be used by
;; default.

(add-hook 'c-mode-common-hook 'thyra-c-mode-common-hook)

;; Provide as package
(provide 'cpp-thyra-styles)
