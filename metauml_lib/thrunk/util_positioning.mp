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

if known _util_positioning_mp:
  expandafter endinput
fi;
_util_positioning_mp:=1;

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

inputonce util_log;

%% Gives a point which is below the given @point with the given @value.
%%
def below(expr point, value) =
    (point - (0, value))
enddef;

%% Gives a point which is above the given @point with the given @value.
%%
def above(expr point, value) =
    (point + (0, value))
enddef;

%% Gives a point which is at the right of the given @point with the given @value.
%%
def atright(expr point, value) =
    (point + (value, 0))
enddef;

%% Gives a point which is at the left of the given @point with the given @value.
%%
def atleft(expr point, value) =
    (point - (value, 0))
enddef;

%% Positions the given objects so that they:
%% * have their tops aligned
%% * the distance between the objects is @distance
%% * the center of gravity of the objects (taken on the top line) is at @middlePoint
%%
%%                         @distance
%%              |---------------------------|
%%              |                           |
%%              |          @middlePoint     |
%%    __________|_____________._____________|____________ the same top
%%      [objectA]             |             [objectB]
%%
vardef centered_align_top(suffix objectA, objectB)(expr distance, middlePoint)=
    log middlePoint;
    objectA.top = objectB.top;
    objectB.left - objectA.right = distance;
    middlePoint = .5[objectB.nw, objectA.ne];
enddef;

%% Deprecated, kept for backward compatibility (align is misspelled).
%%
vardef centered_allign_top(suffix objectA, objectB)(expr distance, middlePoint)=
    log middlePoint;
    objectA.top = objectB.top;
    objectB.left - objectA.right = distance;
    middlePoint = .5[objectB.nw, objectA.ne];
enddef;

%% theString means here not a sequence of characters but
%% the horizontal or vertical line on which objects are placed.
%% For example, if theString is top, then the objects are
%% "hanging" on a horizontal line, like freshly washed 
%% clothes on a string.

vardef align(suffix theString, extremityNew, extremityOld)
                    (text distanceBetweenObjects)(expr sign)(text objects)=
    string objectsAsString__;
    objectsAsString__ := enumToString(objects)("");
    
    log "Alligning '" & objectsAsString__ & "' at " & str theString;
    log sign;
     
    if (objectsAsString__ = ""):
        log "Nothing to do, bailing out.";
    else:
        string previousObject__;
        previousObject__ := "";
        objectIndex__ := 0;
        
        forsuffixes obj = objects:
            log "object: " & str obj;
            if objectIndex__ = 0:
                objectIndex__ := objectIndex__ + 1;
                previousObject__ := str obj;
            else:
                objectIndex__ := objectIndex__ + 1;
                
                log str theString;
                log str extremityOld;
                log str extremityNew;
                log distanceBetweenObjects;
                
                string eqA__, eqB__;
                eqA__ := previousObject__ & "." & str theString & " = " & str obj & "." & str theString;
                
                if (sign = "+"):
                    eqB__ := previousObject__ & "." & str extremityOld & " + " & (str distanceBetweenObjects) & " = " & str obj & "." & str extremityNew;
                    %eqB__ := previousObject__ & "." & str extremityOld & " + distanceBetweenObjects = " & str obj & "." & str extremityNew;
                else:
                    eqB__ := previousObject__ & "." & str extremityOld & " - " & (str distanceBetweenObjects) & " = " & str obj & "." & str extremityNew;
                fi;
                
                log eqA__, eqB__;
                scantokens eqA__;
                scantokens eqB__;

                previousObject__ := str obj;
            fi;
        endfor;
    fi;
enddef;

vardef leftToRight@#(text distanceBetweenObjects)(text objects)=
    if str @# = "":
        log "String is empty, aligning to midy";
        align(midy, left, right)(distanceBetweenObjects)("+")(objects);
    else:
        align(@#,   left, right)(distanceBetweenObjects)("+")(objects);
    fi;
enddef;

vardef rightToLeft@#(text distanceBetweenObjects)(text objects)=
    if str @# = "":
        log "String is empty, aligning to midy";
        align(midy, right, left)(distanceBetweenObjects)("-")(objects);
    else:
        align(@#,   right, left)(distanceBetweenObjects)("-")(objects);
    fi;
enddef;

vardef topToBottom@#(text distanceBetweenObjects)(text objects)=
    if str @# = "":
        log "String is empty, aligning to midx";
        align(midx, top, bottom)(distanceBetweenObjects)("-")(objects);
    else:
        align(@#,   top, bottom)(distanceBetweenObjects)("-")(objects);
    fi;
enddef;

vardef bottomToTop@#(text distanceBetweenObjects)(text objects)=
    if str @# = "":
        log "String is empty, aligning to midx";
        align(midx, bottom, top)(distanceBetweenObjects)("+")(objects);
    else:
        align(@#,   bottom, top)(distanceBetweenObjects)("+")(objects);
    fi;
enddef;

%%
%% Defines an equation which sets the given prefix to be equal for all the given objects.
%% For example:
%%
%% same.top(a, b, c);
%%
%% is equivalent to 
%% a.top = b.top = c.top;
%%
vardef same@#(text objects)=
    log "begin macro: same@#(text objects)=";

    string equation__, property__;
    equation__ := "";
    property__ := str @#;

    objectIndex__ := 0;
    forsuffixes obj = objects:
            log "object: " & str obj;

            if objectIndex__ = 0:
                objectIndex__ := objectIndex__ + 1;
                equation__ := equation__ & str obj & "." & property__;
            else:
                objectIndex__ := objectIndex__ + 1;
                equation__ := equation__ & "=" & str obj & "." & property__;
            fi;
    endfor;

    equation__ := equation__ & ";";

    log "Equation is: '" & equation__ & "'";

    if objectIndex__ <= 1:
        log "One or no objects, nothing to do!";
    else:
        scantokens equation__;
    fi;
enddef;

%% this macro is NOT usable (yet)
vardef leftToRightCentered@#(expr distanceBetweenObjects, middlePoint)(text objects)=
    string objectsAsString__;
    objectsAsString__ := enumToString(objects)("");

    if (objectsAsString__ = ""):
        log "Nothing to do, bailing out.";
    else:
        log "ASDFGHJ";
        leftToRight@#(distanceBetweenObjects)(objects);
    
        string firstObject__, lastObject__;
        objectIndex__ := 0;
        forsuffixes obj = objects:
            if objectIndex__ = 0:
                objectIndex__ := objectIndex__ + 1;
                firstObject__ := str obj;
            else:
                objectIndex__ := objectIndex__ + 1;
                lastObject__ := str obj;
            fi;
        endfor;
        
        log "___";
        .5[authA.nw,authB.ne] = middlePoint;
        log "___";
        
        string eqCenter__;
        eqCenter__ := ".5[" & firstObject__ & ".nw," & lastObject__ & ".ne] = middlePoint";
        
        log eqCenter__;
        log middlePoint;
        
        log "XXX";
        scantokens eqCenter__;
        log "XXX";
    fi;
enddef;
