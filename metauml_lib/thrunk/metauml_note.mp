% MetaUML, a MetaPost library for typesetting exquisite UML diagrams.
% Copyright (C) 2005 Ovidiu Gheorghies
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

if known _metauml_note_mp:
  expandafter endinput
fi;
_metauml_note_mp:=1;

% def metauml_note_debug text txt=log txt enddef;
def metauml_note_debug text txt= enddef;

% Sadly, this copy of the macro is needed to prevent multiple file loads being shown by MetaPost.
% The guard values (such as _metauml_mp) do ensure that the file isn't loaded multiple times, 
% but this macro makes sure that MetaPost won't try to load the file and display a message for that.
def inputonce text libraryFile=
	if not known scantokens ("_" & str libraryFile & "_mp"):
		%includeonce% show "Loading " & str libraryFile;
		scantokens ("input " & str libraryFile);
	else:
		%includeonce% show str libraryFile & " already loaded.";
	fi;
enddef;

inputonce metauml_defaults;
inputonce util_log;

inputonce util_object;
inputonce util_picture_stack;
inputonce util_shade;
inputonce util_margins;

vardef NoteInfo@#=
  attributes(@#);
  var(color)   foreColor, borderColor;
  var(numeric) horizontalFolding, verticalFolding;

  FontInfo.@#iFont(metauml_defaultFont, defaultscale);
  @#iFont.ignoreNegativeBase := 1;

  PictureInfo.@#iText  (2, 2, 2, 2)(@#iFont);
  PictureStackInfo.@#iStack (0, 0, 0, 0)(9)(@#iText);

  @#iStack.iPict.ignoreNegativeBase := 1;

  ShadeInfo.@#iShade;
  Margins.@#(2,2,2,2);

  @#foreColor   := .9white;
  @#borderColor := black;
  @#horizontalFolding := 10;
  @#verticalFolding := 10;
enddef;

vardef NoteInfoCopy@#(text src)=
  attributes(@#);
  var(color)   foreColor, borderColor;

  log "NoteInfoCopy: Copying usecase info";

  PictureStackInfoCopy.@#iStack (src.iStack);
  MarginsCopy.@#(src);
  ShadeInfoCopy.@#iShade(src.iShade);

  @#foreColor := src.foreColor;
  @#borderColor := src.borderColor;
  @#horizontalFolding := src.horizontalFolding;
  @#verticalFolding := src.verticalFolding;
enddef;

NoteInfo.iNote;

vardef Note@#(text contents)=
  ENote.@#(iNote)(contents);
enddef;

vardef ENote@#(text _info)(text contents)=
  ObjectEquations(@#);
  @#className := "Note";

  NoteInfoCopy.@#info(_info);

  %% The following code cannot work when pictures are given in the contents argument, such as in:
  %% ENote.myNote(iNote)("This is a formula", TEX("$a=b$"));
  %%
  %numeric @#nLines; @#nLines := 0;
  %string @#lines[];
  %for l=contents:
  %  @#lines[@#nLines] := l;
  %  @#nLines := @#nLines + 1;
  %endfor;
  %
  %EPictureStack.@#stack(@#info.iStack)(scantokens listArray(@#lines)(@#nLines))("vleftbase");

  EPictureStack.@#stack(@#info.iStack)(contents)("vleftbase");
enddef;

vardef Note_layout@#=
  if @#layedout = 1:
    log "Note " & (str @#) & " has already been layed out";
  else:
    @#layedout := 1;
    log "Note layout: " & (str @#);

    PictureStack_layout.@#stack;

    @#width = @#stack.width + @#info.left + @#info.right + @#info.horizontalFolding;
    @#height = @#stack.height + @#info.top + @#info.bottom;

    @#top = @#stack.top + @#info.top;
    @#left = @#stack.left - @#info.left;
  fi;
enddef;

vardef Note_draw@#=
  Note_layout@#;
  objectEnsurePositioning.@#;

  attributes(@#);
  var(path) border, corner;

  @#border := @#nw -- (@#ne - (@#info.horizontalFolding, 0)) --
                      (@#ne - (0, @#info.verticalFolding)  ) -- @#se -- @#sw -- cycle;

  @#corner := (@#ne - (@#info.horizontalFolding, 0)) --
	     (@#ne - (@#info.horizontalFolding, @#info.verticalFolding)) --
             (@#ne - (0, @#info.verticalFolding));

  drawObjectShade(@#);

  fill Note_border.@# withcolor @#info.foreColor;
  draw Note_border.@# withcolor @#info.borderColor;

  draw @#corner withcolor @#info.borderColor;
  drawObject(@#stack);
enddef;

vardef Note_border@#=
  @#border
enddef;

