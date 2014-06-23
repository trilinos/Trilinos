(setq custom-file "~/.emacs.d/custom.el")
(load custom-file)

(setq my-emacs "~/.emacs.d")
(add-to-list 'load-path my-emacs)
(load-library "tme-keys.el")

(setq doxymacs-path "/usr/local/src/doxymacs-1.8.0/share/emacs/site-lisp")
(add-to-list 'load-path doxymacs-path)
(setq nemesis-setup-doxymacs t)

(setq default-major-mode 'text-mode)
(add-hook 'text-mode-hook 'turn-on-auto-fill)

(setq nemesis-dir "/vendors/environment/emacs")
(add-to-list 'load-path nemesis-dir)
(load-library "nemesis")

(add-to-list 'initial-frame-alist '(width . 80))
(add-to-list 'initial-frame-alist '(height . 70))
(add-to-list 'default-frame-alist '(width . 80))
(add-to-list 'default-frame-alist '(height . 70))

(setq auctex-path "/usr/local/src/auctex-11.85/share/emacs/site-lisp")
(add-to-list 'load-path auctex-path)
(load "auctex.el" nil t t)

(defun mycompilation-mode-hook ()
  (setq truncate-lines nil))
(add-hook 'compilation-mode-hook 'mycompilation-mode-hook)

(defun mydired-mode-hook ()
  (setq truncate-lines nil))
(add-hook 'dired-mode-hook 'mydired-mode-hook)

(load-file "/usr/local/src/cedet-1.0pre7/common/cedet.el")
(global-ede-mode 1)
(semantic-load-enable-code-helpers)
(global-srecode-minor-mode 1)

(add-to-list 'load-path "/usr/local/src/ecb-2.40")
(require 'ecb)
