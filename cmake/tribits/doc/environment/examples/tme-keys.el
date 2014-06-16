;;---------------------------------------------------------------------------;;
;; file   : tme-keys.el
;; author : Thomas M. Evans
;; date   : Monday January 25 23:24:7 2010
;; brief  : Key mappings for Emacs environment.
;;---------------------------------------------------------------------------;;

(global-set-key (kbd "<f7>")          'nemesis-save-and-kill-buffer)
(global-set-key (kbd "<mouse-3>")     'kill-region)
(global-set-key (kbd "<ESC> <f1>")    'overwrite-mode)
(global-set-key (kbd "<f8>")          'nemesis-toggle-previous-buffer)
(global-set-key (kbd "C-*")           'start-kbd-macro)
(global-set-key (kbd "C-+")           'end-kbd-macro)
(global-set-key (kbd "C-<kp-enter>")  'call-last-kbd-macro)

(defun tme-dired-mode-hook ()
  (local-set-key (kbd "<mouse-2>") 'dired-find-file))
(add-hook 'dired-mode-hook 'tme-dired-mode-hook)

;;---------------------------------------------------------------------------;;
;; REMOVE TRAILING ^M
;;---------------------------------------------------------------------------;;

(defun remove-or-convert-trailing-ctl-M ()
  "Optionally remove or convert trailing ^M from a file."
  (interactive)
  (save-excursion
    (goto-char (point-min))
    (if (search-forward "\^M" nil t)
        (if (or (= (preceding-char) ?\^J)
                (= (following-char) ?\^J))
            (if (y-or-n-p "Remove trailing ^M ? ")
                (progn (goto-char (point-min))
                       (perform-replace "\^M" "" nil nil nil)
                       (pop-mark))
              (message "No transformation."))
          (if (y-or-n-p "Convert ^M into ^J ? ")
              (progn (goto-char (point-min))
                     (perform-replace "\^M" "\^J" nil nil nil)
                     (pop-mark))
            (message "No transformation."))))))

;;---------------------------------------------------------------------------;;
;; end of tme-keys.el
;;---------------------------------------------------------------------------;;
