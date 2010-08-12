% MetaUtil, an easier MetaPost life
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

if known _util_object_mp:
  expandafter endinput
fi;
_util_object_mp:=1;

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

inputonce util_infrastructure;
inputonce util_log;

def ObjectEquations(suffix $)=
  attributes($);
  var(pair) n, ne, e, se, s, sw, w, nw, c;
  var(pair) dim;
  var(numeric) width, height;

  var(numeric) left, right;
  var(numeric) top,  bottom;
  var(numeric) midx, midy;

  var(picture) pict;
  var(string) className;

  var(numeric) layedout, drawn;

  $.layedout := 0;
  $.drawn := 0;

  $.left  = xpart $.nw = xpart $.sw;
  $.right = xpart $.ne = xpart $.se;

  $.top    = ypart $.ne = ypart $.nw;
  $.bottom = ypart $.se = ypart $.sw;

  $.w = .5[$.nw,$.sw];
  $.s = .5[$.sw,$.se];
  $.e = .5[$.ne,$.se];
  $.n = .5[$.ne,$.nw];

  $.c = .5[$.nw, $.se];

  $.width  = $.right - $.left;
  $.height = $.top - $.bottom;

  $.midx = xpart($.c);
  $.midy = ypart($.c);

  xpart $.dim = $.width;
  ypart $.dim = $.height;

  $.className := "Object";
enddef;

def objectBox(text obj)=
  obj.nw -- obj.ne -- obj.se -- obj.sw -- cycle
enddef;

vardef objectBorder(suffix @#)=
  save methodName;
  string methodName;
  methodName := @#className & "_border";
  log "invoking " & methodName & " arg=" & str @#;

  scantokens (methodName).@#
enddef;

vardef Object_toString(text obj)=
  save ret; string ret;
  ret := "Generic object";
  ret
enddef;

vardef toString(suffix @#)=
  scantokens (@#className & "_toString").@#
enddef;

vardef layoutObject(suffix @#)=
  save methodName;
  string methodName;
  methodName := @#className & "_layout";
  log "invoking " & methodName & " arg=" & str @#;

  scantokens (methodName).@#;
enddef;

vardef drawObject(suffix @#)=
  save methodName;
  string methodName;
  methodName := @#className & "_draw";
  log "invoking " & methodName & " arg=" & str @#;

  scantokens (methodName).@#;
enddef;

vardef drawObjectShade(suffix @#)=
  fill objectBorder(@#) shifted
       (@#.info.iShade.shift, -@#.info.iShade.shift) withcolor @#.info.iShade.background;
enddef;

vardef objectEnsurePositioning@#=
    if unknown (xpart @#nw):
    log "*** WARNING: DEFAULTING OBJECT'S (" & (str @#) & ") x TO 0";
    xpart @#nw = 0;
  fi;
  if unknown (ypart @#nw):
    log "*** WARNING: DEFAULTING OBJECT'S (" & (str @#) & ") y TO 0";
    ypart @#nw = 0;
  fi;
enddef;

vardef drawObjectAt(suffix @#)(text eq)=
  eq;
  drawObject(@#);
enddef;

vardef layoutObjects(text objects)=
  log "Layout objects...";
  forsuffixes o = objects:
    if (str o) <> "":
      layoutObject(o);
    fi;
  endfor;
enddef;

vardef drawObjects(text objects)=
  % This needs to be done first in order to allow for equations which define
  % relative positioning to be appropiately solved. Otherwise, for example,
  % object0 may be positioned somewhere by default, and object1 also, leading
  % to inconsitent equations. Inconsistency would have been caused by the fact that the
  % positioning equations become meaningful only after object layout is performed
  % --- but that would have been too late, since the positioning would have been
  % done already to inappropriate values.
  %
  layoutObjects(objects);
  
  forsuffixes o = objects:
    if (str o) <> "":
      drawObject(o);
    fi;
  endfor;
enddef;

%%
%% Joins the given pictures according to the equation given
%% by the pictureDoJoin macro. Note: the pictureDoJoin macro
%% can be conveniently set by calling setPicutureJoin macro.
%%
%% Usage example:
%% joinObjects (p, q, r);
%%
def joinObjects(text pictures)=
  save skipFirstPicture;
  skipFirstPicture := 1;
  forsuffixes p=pictures:
    if skipFirstPicture=0:
      objectDoJoin(previousPic, p);
    else:
      % first picture is skipped
      skipFirstPicture := 0;
    fi;
    def previousPic=p enddef;
  endfor;
enddef;

%%
%% Sets the macro pictureDoJoin, which is used by joinObjects to join
%% one picture into another.
%%
%% The pictureDoJoin macro has two arguments: pa and pb.
%% pa represents the first picture and pb the second picture.
%% The second picture is positioned relatively to pa as specified
%% by the given equation.
%%
%% Usage examples:
%%
%% % place the second picture at the bottom of the first
%% setPictureJoin(pa.sw = pb.nw);
%%
%% % place the second picture at the right of the first
%% setPictureJoin(pa.ne = pb.nw);
%%
def setObjectJoin(text eq)=
  def objectDoJoin(suffix pa, pb)=eq enddef;
enddef;

%%
%% A shortcut for setting the join equation and actually performing the join.
% 
def joinObjectsEq(text objects)(text eq)=
    setObjectJoin(eq);
    joinObjects(objects);
enddef;


%%
%% Joins the pictures and draws them. A shorthand for
%% the corresponding two macro invocations.
%%
def joinDrawObjects(text pictures)=
  joinObjects(pictures);
  drawObjects(pictures);
enddef;
