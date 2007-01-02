;;; C++ font-lock'ing
;;
;; Copyright (C) 1992, 93, 94, 95, 96, 97, 98, 1999, 2000, 2001
;;  Free Software Foundation, Inc.
;; Copyright (C) 2002 Brad Settlemyer
;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program; if not, you can either send email to this
;; program's author (see below) or write to:
;;
;;              The Free Software Foundation, Inc.
;;              675 Mass Ave.
;;              Cambridge, MA 02139, USA.
;;
;; Send bug reports to Brad Settlemyer <bws@deepcopy.org>

;; This code is derived from the FSF emacs package, font-lock.el.  A few
;; modifications have been made to simplify, remove redundancy, improve
;; support for certain constructs, and seperate this code into a smaller
;; module.

;;; Include and initialize dependencies
;; From font-lock: font-lock-constant-face font-lock-function-name-face
;;                 font-lock-keyword-face font-lock-preprocessor-face
;;                 font-lock-string-face font-lock-type-face 
;;                 font-lock-variable-name-face font-lock-defaults
(require 'font-lock)

;; From cc-mode: c-initialize-cc-mode c-at-toplevel-p
(require 'cc-mode)
(c-initialize-cc-mode)

;; Common lisp
(require 'cl)

;;; cpp-font-lock definitions
;; Version information
(defconst cpp-font-lock-version "0.4.0" "C++ Font Lock version number.")

;; Variables used by cpp-font-lock to determine which faces to use
;; when highlight C++ code.  Initially, these are bound to the
;; cpp-font-lock face of the same name (ala' GNU Emacs).
(defvar cpp-font-lock-comment-face 'cpp-font-lock-comment-face
  "Face name to use for comments.")

(defvar cpp-font-lock-preprocessor-face 'cpp-font-lock-preprocessor-face
  "Face name to use for preprocessor directives.")

(defvar cpp-font-lock-function-name-face 'cpp-font-lock-function-name-face
  "Face name to use for function definitions.")

(defvar cpp-font-lock-keyword-face 'cpp-font-lock-keyword-face
  "Face name to use for C++ keywords.")

(defvar cpp-font-lock-literal-face 'cpp-font-lock-literal-face
  "Face name to use for C++ literals.")

(defvar cpp-font-lock-preprocessor-face 'cpp-font-lock-preprocessor-face
  "Face name to use for preprocessor directives.")

(defvar cpp-font-lock-type-name-face 'cpp-font-lock-type-name-face
  "Face name to use for type definitions and declarations.")

(defvar cpp-font-lock-variable-name-face 'cpp-font-lock-variable-name-face
  "Face name to use for type definitions and declarations.")


;; Custom faces to use for cpp-font-lock.  Initialized from the standard
;; font lock faces
(defgroup cpp-font-lock-faces nil
  "Faces used by the cpp-font-lock package."
  :prefix "cpp-font-lock"
  :group 'font-lock-faces)

(defface cpp-font-lock-comment-face 
  '((((class color)(background light)) (:foreground "FireBrick"))
    (((class color)) (:foreground "LightPink"))
    (t (:inverse-video t)))
  "Font Lock face used to highlight literals."
  :group 'cpp-font-lock-faces)

(defface cpp-font-lock-function-name-face 
  '((((class color)(background light)) (:foreground "DarkBlue"))
    (((class color)) (:foreground "NavajoWhite1"))
    (t (:inverse-video t)))
  "Font Lock face used to highlight function names."
  :group 'cpp-font-lock-faces)

(defface cpp-font-lock-keyword-face 
  '((((class color)(background light)) (:foreground "Blue" :bold t))
    (((class color)) (:foreground "SkyBlue" :bold t))
    (t (:inverse-video t)))
  "Font Lock face used to highlight keywords."
  :group 'cpp-font-lock-faces)

(defface cpp-font-lock-literal-face 
  '((((class color) (background light)) (:foreground "ForestGreen"))
    (((class color)) (:foreground "SpringGreen"))
    (t (:inverse-video t)))
  "Font Lock face used to highlight literals."
  :group 'cpp-font-lock-faces)

(defface cpp-font-lock-preprocessor-face 
  '((((class color)(background light)) (:foreground "gray40"))
    (((class color)) (:foreground "cyan2"))
    (t (:inverse-video t)))
  "Font Lock face used to highlight preprocessor directives."
  :group 'cpp-font-lock-faces)

(defface cpp-font-lock-type-name-face 
  '((((class color)(background light)) (:foreground "maroon4"))
    (((class color)) (:foreground "DarkSeaGreen1"))
    (t (:inverse-video t)))
  "Font Lock face used to highlight type names."
  :group 'cpp-font-lock-faces)

(defface cpp-font-lock-variable-name-face 
  '((((class color)(background light)) (:foreground "SaddleBrown"))
    (((class color)) (:foreground "NavajoWhite1"))
    (t (:inverse-video t)))
  "Font Lock face used to highlight variable names."
  :group 'cpp-font-lock-faces)


;; Abstracts Emacs incompatibilities for locating face objects
(defun cpp-font-lock-find-face (face-name)
  (if (featurep 'xemacs) 
      (find-face face-name)
    (internal-find-face face-name)))


;; Function exposed for setting cpp-font-lock font settings to the
;; default font-lock faces
(defun cpp-font-lock-use-font-lock-faces (toggle)
  (if toggle 
      (progn (setq cpp-font-lock-comment-face font-lock-comment-face)
             (setq cpp-font-lock-function-name-face 
                   font-lock-function-name-face)
             (setq cpp-font-lock-keyword-face font-lock-keyword-face)
             (setq cpp-font-lock-literal-face font-lock-constant-face)
             (setq cpp-font-lock-preprocessor-face font-lock-preprocessor-face)
             (setq cpp-font-lock-type-name-face font-lock-type-face)
             (setq cpp-font-lock-variable-name-face 
                   font-lock-variable-name-face))
    (progn (setq cpp-font-lock-comment-face 
                 (cpp-font-lock-find-face 'cpp-font-lock-comment-face))
           (setq cpp-font-lock-function-name-face
                 (cpp-font-lock-find-face 'cpp-font-lock-function-name-face))
           (setq cpp-font-lock-keyword-face
                 (cpp-font-lock-find-face 'cpp-font-lock-keyword-face))
           (setq cpp-font-lock-literal-face
                 (cpp-font-lock-find-face 'cpp-font-lock-literal-face))
           (setq cpp-font-lock-variable-name-face
                 (cpp-font-lock-find-face 'cpp-font-lock-variable-name-face))
           (setq cpp-font-lock-type-name-face
                 (cpp-font-lock-find-face 'cpp-font-lock-type-name-face)))))
  

;;; Binding's for C++ font-lock'ing
;; Regular expression to match c++ whitespace
(defconst cpp-font-lock-ws-re "[ \t\n]")


;; Regular expression to match identifiers (inclusion of : and ~ is probably
;; going to cause problems eventually
(defconst cpp-font-lock-word-re "[a-zA-Z_:~][a-zA-Z0-9_:\\-]*")


;; Initialized to the standard (and semi-standard) library types and 
;; types ending in _t (reserved in Standard C++ for library typedefs).
;; Standard types that have been deprecated are not included.
(defcustom cpp-font-lock-extra-types
  (concat 
   cpp-font-lock-word-re "_t\\|"
   (eval-when-compile
     (regexp-opt 
      '("auto_ptr" "bitset" "bit_vector" "complex" "const_iterator" 
        "const_reverse_iterator" "deque"
        "hash_map" "hash_multimap" "hash_set" "hash_multiset" 
        "ifstream" "iostream" "ios" "istream" "istringstream"
        "iterator" "list" "map" "multimap" "multiset" "ofstream" 
        "ostringstream" "ostream" "priority_queue" "queue" 
        "reverse_iterator" "set" "size_type" "slist" "stack" 
        "streambuf" "string" "type_info" "valarray" "vector" "wstring"))))
  "*List of extra types to fontify in C++ mode.
Each list item should be a regexp not containing word-delimiters.
For example, a value of (\"string\") means the word string is treated as a type
name.

The value of this variable is used when Font Lock mode is turned on."
  :type 'font-lock-extra-types-widget
  :group 'font-lock-extra-types)


;; Regexp matches after point:		word<word>::word (
;;						^^^^ ^^^^   ^^^^ ^
;; Where the match subexpressions are:	  1    3      5  6
;;
;; Item is delimited by (match-beginning 1) and (match-end 1).  
;; If (match-beginning 3) is non-nil, that part of the item incloses a `<>'.
;; If (match-beginning 5) is non-nil, that part of the item follows a `::'.
;; If (match-beginning 6) is non-nil, the item is followed by a `('.
(defun cpp-font-lock-match-c++-style-declaration-item-and-skip-to-next (limit)
  (when (looking-at 
         (eval-when-compile
           (concat
            ;; Skip any leading whitespace.
            "[ \t\n*&]*"
            ;; This is `c++-type-spec' from below.  (Hint hint!)
            "\\(\\sw+\\)"				; The instance?
            "\\([ \t\n]*<\\(\\(?:[^<>]\\|<[^>]+>\\)+\\)[ \t\n*&]*>\\)?"	; Or template?
            "\\([ \t\n]*::[ \t\n*~]*\\(\\sw+\\)\\)*"	; Or member?
            ;; Match any trailing parenthesis.
            "[ \t\n]*\\((\\)?")))
    (save-match-data
      (condition-case nil
	  (save-restriction
	    ;; Restrict to the end of line, currently guaranteed to be LIMIT.
	    (narrow-to-region (point-min) limit)
	    (goto-char (match-end 1))
	    ;; Move over any item value, etc., to the next item.
	    (while (not (looking-at "[ \t\n]*\\(\\(,\\)\\|;\\|\\'\\)"))
	      (goto-char (or (scan-sexps (point) 1) (point-max))))
	    (goto-char (match-end 2)))
	(error t)))))


;; FIXME
;; Matches constructors and destructors not preceded by keywords.
;; The regular expression can be logically partitioned into the following
;; sub expressions (which have been somewhat optimized):
;;  sub-ex 1 - constructor/destructor name
(defun cpp-font-lock-match-structor-declaration (limit)
  (let ((end-match-point nil)
        (regexp (concat "^\\s-+"
                        ;"\\(\\sw+\\|\\sw+\\s-+\\)(")))
                        "\\(~?" cpp-font-lock-word-re "\\s-*\\)(")))
    (while (and (setq end-match-point (re-search-forward regexp limit t))
                (and (save-excursion
                       (beginning-of-line)
                       (save-match-data
                         (not (vectorp (c-at-toplevel-p))))))))
    end-match-point))


;; Matches all member functions and constructors/destructors preceded by
;; keywords (they are matched in the return type sub-ex).
;; The regular expression can be logically partitioned into the following 
;; sub expressions (which have been somewhat optimized):
;;  sub-ex 1 - optional inline
;;  sub-ex 2 - optional function modifiers
;;  sub-ex 3 - optional const mod of return type
;;  sub-ex 4 - return type wrapper
;;  sub-ex 5 - return type name
;;  sub-ex 6 - optional return type modifiers (*, &, const)
;;  sub-ex 7 - member function name
(defun cpp-font-lock-match-member-function-declaration (buffer-size)
  (let ((end-match-point nil)
        (match-function-decl-re
         (eval-when-compile
           (concat "^\\s-+"
                   "\\(inline\\s-+\\)?"
                   "\\(static\\s-+\\|virtual\\s-+\\)?"
                   "\\(const\\s-+\\)?"
                   "\\(" cpp-font-lock-word-re "\\(*\\|&\\|\\s-\\|<\\)\\)"
                   "\\(*\\|&\\|const\\s-+\\|\\s-\\|.+>\\|<.+>\\)*"
                   "\\(" cpp-font-lock-word-re "\\)\\s-*"
                   "(\\s-*"
                   "\\(" cpp-font-lock-word-re "\\s-*" cpp-font-lock-word-re
                   "\\|)\\|\n\\)"))))
    (while (and (setq end-match-point 
                      (re-search-forward match-function-decl-re buffer-size t))
                (save-excursion
                  (beginning-of-line)
                  (save-match-data
                    (not (vectorp (c-at-toplevel-p)))))))
    end-match-point))
    ;(if (save-excursion (not (vectorp (c-at-toplevel-p))))
        ;(re-search-forward match-function-decl-re buffer-size t)
;        (re-search-forward match-function-decl-re buffer-size t)
;))


;;FIXME
;; Matches variable declarations (at any scope)
;;  sub-ex 1 - optional static modifier
;;  sub-ex 2 - optional const modifier
;;  sub-ex 3 - type wrapper
;;  sub-ex 4 - pointer or reference declarations
;;  sub-ex 5 - optional pointer const
;;  sub-ex 6 - variable name wrapper
;;  sub-ex 7 - variable name declaration
(defun cpp-font-lock-match-variable-declarations (buffer-size)
  (let ((match-var-decl-re 
         (eval-when-compile
           (concat "^\\s-+"
                   "\\(static\\s-+\\)?"
									 "\\(.*[,(]\\s-+\\)?"
                   "\\(const\\s-+\\)?"
                   "\\(" cpp-font-lock-word-re "\\(*\\|&\\|\\s-\\|<\\)\\)"
                   "\\(*\\|&\\|const\\s-+\\|\\s-\\|.+>\\|<.+>\\)*"
                   "\\(\\(" cpp-font-lock-word-re "\\)\\s-*[,;)]\\s-*\\)+"))))
    (re-search-forward match-var-decl-re buffer-size t)))

    
;; Matches all C++ literals
;; Beginning match is there to avoid case where literal is embedded in
;; an identifier
;;  sub-ex 1 - literal
(defun cpp-font-lock-match-literal (buffer-size)
  (let ((literal-re (eval-when-compile 
                      (concat
                       "[ \t\n,;='<>()&-/^~\+\*\|\"\[]\\("
                       (regexp-opt '("true" "false"))       ; bool
                       "\\|[-]?[0-9]+[eE][-]?[0-9]+[fFlL]?" ; exp lit
                       "\\|[-]?[0-9]+\\.[0-9]+[eE][-]?[0-9]+[fFlL]?" ; dec exp
                       "\\|0[xX][0-9a-fA-F]+"               ; hex
                       "\\|0[0-7]+"                         ; octal
                       "\\|[-]?[0-9]+\\.[0-9]+[fFlL]?"      ; dec
                       "\\|[-]?[1-9]?[0-9]+[uUlL]?"        ; int
                       "\\)"))))
    (re-search-forward literal-re buffer-size t)))


(let* 
    ;;
    ;; C++ Keywords from C++PL 3rd ed. (Stroustrup) Appendix A.2 
    ;; excludes builtin types and keyword operator (it should be fontified
    ;; as a function)
    ((c++-keywords
      (eval-when-compile
        (regexp-opt
         '("and" "and_eq" "asm" "auto" "bitand" "bitor" "break" "case" "catch"
           "class" "compl" "const" "const_cast" "continue" "default" "delete" 
           "do" "dynamic_cast" "else" "enum" "explicit" "export" "extern" 
           "for" "friend" "goto" "if" "inline" "mutable" "namespace" "new" 
           "not" "not_eq" "or" "or_eq" "public" "private" 
           "protected" "register" "reinterpret_cast" "return" "sizeof" 
           "static" "static_cast" "struct" "switch" "template" "this"
           "throw" "try" "typedef" "typeid" "typename" "union" "using" 
           "virtual" "volatile" "while" "xor" "xor_eq"))))
     ;;
     ;; C++ builtin types from C++PL 3rd ed. Appendix A.2
     (c++-type-names
      (eval-when-compile
        (regexp-opt
         '("bool" "char" "double" "float" "int" "long" "short" "signed" 
           "unsigned" "void" "wchar_t") t)))
     ;(c++-type-names-depth `(regexp-opt-depth ,c++-type-names))
     ;;
     ;; C++ operators from C++PL Appendix A.8.3 (excl. previous keywords)
     (c++-operators
      (eval-when-compile
        (regexp-opt
         '("+" "-" "*" "/" "%" "^" "&" "|" "~" "!" "=" "<" ">" "+=" "-="
           "*=" "/=" "%=" "^=" "&=" "|=" "<<" ">>" ">>=" "<<=" "==" "!="
           "<=" ">=" "&&" "||" "++" "--" "->*" "," "->" "[]" "()"))))
     ;;
     ;; Matches all C++ type specifications
     (c++-type-specs
      (eval-when-compile
        (regexp-opt
         '("class" "enum" "namespace" "struct" "union") t)))
     (c++-type-specs-depth (regexp-opt-depth c++-type-specs))
     ;;
     ;; A brave attempt to match templates following a type and/or match
     ;; class membership.  See and sync the above function
     ;; `font-lock-match-c++-style-declaration-item-and-skip-to-next'.
     (c++-type-suffix (concat "\\([ \t]*<\\(\\(?:[^<>\n]\\|<[^>\n]+>\\)+\\)[ \t*&]*>\\)?"
                              "\\([ \t]*::[ \t*~]*\\(\\sw+\\)\\)*"))
     (c++-type-suffix-depth (regexp-opt-depth c++-type-suffix))
     ;;
     ;; If the string is a type, it may be followed by the cruft above.
     (c++-type-spec (concat "\\(\\sw+\\)\\>" c++-type-suffix))
     ;;
     ;; Parenthesis depth of user-defined types not forgetting their cruft.
     ;(c++-type-depth `(regexp-opt-depth
     ;                  (concat ,c++-type-names ,c++-type-suffix)))
     ;;
     ;; Directives recognized by the C/C++ preprocessor
     (c-preprocessor-directives
	(eval-when-compile
	  (regexp-opt
	   '("define"  "elif" "else" "endif" "error" "file" "if" "ifdef"
	     "ifndef" "include" "line" "pragma" "undef"))))
     (c-preprocessor-directives-depth 
      (regexp-opt-depth c-preprocessor-directives))
     )
  ;;
  ;; Fontify strings, preprocessor cmds (incl macro definitions),
  ;; function names in function declarations (not within type declarations)
  (defconst cpp-font-lock-keywords-1
    (list
     ;;
     ;; Fontify filenames in #include <...> preprocessor directives as strings.
     '("^#[ \t]*\\(import\\|include\\)[ \t]*<\\([^>\"\n]*\\)>?"
       2 cpp-font-lock-literal-face)
     ;;
     ;; Fontify function macro names.
     '("^#[ \t]*define[ \t]+\\(\\sw+\\)(" 1 cpp-font-lock-function-name-face)
     ;;
     ;; Fontify symbol names in #if ... defined preprocessor directives.
     '("^#[ \t]*if\\>"
       ("\\<\\(defined\\)\\>[ \t]*(?\\(\\sw+\\)?" nil nil
        (1 cpp-font-lock-preprocessor-face) 
        (2 cpp-font-lock-variable-name-face nil t)))
     ;;
     ;; Fontify symbol names in #elif ... defined preprocessor directives.
     '("^#[ \t]*elif\\>"
       ("\\<\\(defined\\)\\>[ \t]*(?\\(\\sw+\\)?" nil nil
        (1 cpp-font-lock-preprocessor-face) 
        (2 cpp-font-lock-variable-name-face nil t)))
     ;;
     ;; Fontify preprocessor directives as preprocessor and following symbols 
     ;; as variables
     (list
      (concat "^#[ \t]*\\(" c-preprocessor-directives
              "\\)\\>[ \t!]*\\(" cpp-font-lock-word-re "\\)?")
      '(1 cpp-font-lock-preprocessor-face)
      (list (+ 2 c-preprocessor-directives-depth)
            'cpp-font-lock-variable-name-face nil t))
     ;;
     ;; Fontify function name definitions, possibly incorporating class names.
     (list (concat "^" c++-type-spec "[ \t]*(")
           '(1 (if (or (match-beginning 2) (match-beginning 4))
                   cpp-font-lock-type-name-face
                 cpp-font-lock-function-name-face))
           '(3 cpp-font-lock-type-name-face nil t)
           '(5 cpp-font-lock-function-name-face nil t))
     )
    "Subdued level highlighting for C++ mode.")

  ;; Fontify builtin types, extra-types, literals, and keywords
  (defconst cpp-font-lock-keywords-2
    (append cpp-font-lock-keywords-1
     (list
      ;;
      ;; Fontify builtin and extra types
      (cons (concat "\\<\\(" c++-type-names "\\|" 
                    cpp-font-lock-extra-types "\\)\\>")
            '(0 cpp-font-lock-type-name-face))
      ;;
      ;; Fontify macro definitions
      ;(list 'cpp-font-lock-match-macro 
      ;      '(2 font-lock-function-name-face nil t))
      ;      '(4 font-lock-variable-name-face nil t)
      ;      '(8 font-lock-variable-name-face nil t))
      ;;
      ;; Fontify literal expressions
      '(cpp-font-lock-match-literal 1 cpp-font-lock-literal-face)
      ;;
      ;; Fontify other builtin keywords.
      (cons (concat "\\<\\(" c++-keywords "\\)\\>") 
            '(0 cpp-font-lock-keyword-face))
      ;;
      ;; Fontify operator overloading methods
      (cons (concat "\\<\\(operator\\)\\>[ \t]*\\(" c++-operators "\\)?")
            '(. cpp-font-lock-function-name-face))
      ;;
      ;; Fontify some constructors and destructors
      '(cpp-font-lock-match-structor-declaration 
        1 cpp-font-lock-function-name-face)
      ;;
      ;; Fontify member functions
      '(cpp-font-lock-match-member-function-declaration 
        7 cpp-font-lock-function-name-face)
      ))
  "Medium level highlighting for C++ mode.
See also `c++-font-lock-extra-types'.")

  ;; Fontify variable declarations, class name declarations, class member
  ;; function declarations, class member data declarations, and
  ;; constructors and destructors
  (defconst cpp-font-lock-keywords-3
    (append cpp-font-lock-keywords-2
     (list
      ;;
      ;; Fontify variable declarations
      (cons 'cpp-font-lock-match-variable-declarations 
            '(8 cpp-font-lock-variable-name-face nil t))
      ;;
      ;; Fontify type (class, struct, and enum) declaration names
      `(,(concat "\\<" c++-type-specs "\\>" c++-type-suffix
                 "[ \t]*\\(" c++-type-spec "\\)?")
        ;; The name of any template type.
        (,(+ c++-type-specs-depth 2) 'cpp-font-lock-type-name-face nil t)
        ;; The name of any type.
        (,(+ c++-type-specs-depth c++-type-suffix-depth 2)
         cpp-font-lock-type-name-face nil t))
      ;;
      ;; Fontify anything at beginning of line as a type declaration or 
      ;; a function definition.
      ;; This fontifies class name declarations as variable declarations
      ;; This fontifies the class type for the additional cases
      ;;   inline Class::~Class()
      ;;   Class& Class::operator=( const Class& other )
      `(,(concat "^\\(" c++-type-spec "[ \t*&]*\\)+")
        (cpp-font-lock-match-c++-style-declaration-item-and-skip-to-next
         (prog1 (progn (skip-chars-forward "^;{}") (point))
           (goto-char (match-beginning 1)))
         (goto-char (match-end 1))
         (1 (cond ((or (match-beginning 2) (match-beginning 4))
                   cpp-font-lock-type-name-face)))
         (5 (if (match-beginning 6)
                cpp-font-lock-function-name-face
              cpp-font-lock-variable-name-face) nil t)
         ))

      ))
    "Gaudy level highlighting for C++ mode.
See also `cpp-font-lock-extra-types'.")
  )

(defvar cpp-font-lock-keywords cpp-font-lock-keywords-3
  "Default expressions to highlight in C++ mode.
See also `cpp-font-lock-extra-types'.")


;; XEmacs: Replace the provided c++ font-lock with this package
(put 'c++-mode 'font-lock-defaults
     '((cpp-font-lock-keywords
	cpp-font-lock-keywords-1 cpp-font-lock-keywords-2
	cpp-font-lock-keywords-3)
       nil nil ((?_ . "w") (?~ . "w")) beginning-of-defun))
(make-local-variable 'font-lock-defaults)

;; GNU Emacs: Replace the provided c++ font-lock with this package
(defun cpp-font-lock-defaults ()
  (setq font-lock-defaults
        '( (cpp-font-lock-keywords
            cpp-font-lock-keywords-1
            cpp-font-lock-keywords-2
            cpp-font-lock-keywords-3)
           nil
           nil
           ((95 . "w")) 
           beginning-of-defun
           (font-lock-mark-block-function . mark-defun))))
(add-hook 'c++-mode-hook 'cpp-font-lock-defaults)


;; Turn on C++ font-lock'ing for C++ mode
(add-hook 'c++-mode-hook 'font-lock-fontify-buffer)


;; Provide as package
(provide 'cpp-font-lock)
