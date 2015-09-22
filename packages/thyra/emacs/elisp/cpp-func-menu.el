;;; ********************
;;; func-menu is a package that scans your source file for function
;;; definitions and makes a menubar entry that lets you jump to any
;;; particular function definition by selecting it from the menu.  The
;;; following code turns this on for all of the recognized languages.
;;; Scanning the buffer takes some time, but not much.
;;;
;;; Send bug reports, enhancements etc to:
;;; David Hughes <ukchugd@ukpmr.cs.philips.nl>
;;;
(cond (running-xemacs
  (require 'func-menu)
  (define-key global-map 'f8 'function-menu)
  (add-hook 'find-file-hooks 'fume-add-menubar-entry)
  (define-key global-map "\C-cl" 'fume-list-functions)
  (define-key global-map "\C-cg" 'fume-prompt-function-goto)

  ;; The Hyperbole information manager package uses (shift button2) and
  ;; (shift button3) to provide context-sensitive mouse keys.  If you
  ;; use this next binding, it will conflict with Hyperbole's setup.
  ;; Choose another mouse key if you use Hyperbole.
  (define-key global-map '(shift button3) 'mouse-function-menu)

  ;; For descriptions of the following user-customizable variables,
  ;; type C-h v <variable>
  (setq fume-max-items 25
        fume-fn-window-position 3
        fume-auto-position-popup t
        fume-display-in-modeline-p t
        fume-menubar-menu-location "File"
        fume-buffer-name "*Function List*"
        fume-no-prompt-on-valid-default nil)
  )
)

;; Provide this menu stuff as a lisp package
(provide 'cpp-func-menu)
