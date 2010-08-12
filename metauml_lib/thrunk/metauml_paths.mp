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

if known _metauml_paths_mp: 
  expandafter endinput 
fi;
_metauml_paths_mp:=1;

def pathDirect(expr pA,pB) =
  pA--pB
enddef;

vardef pathManhattanX(text twoParams) =
  pair pointA__, pointB__;
  string  paramA_str__, paramB_str__;
  boolean paramA_isObject__, paramB_isObject__;
  
  paramA_isObject__ := true;
  paramB_isObject__ := true;
  
  boolean isFirstParam__;
  isFirstParam__ := true;

  %string twoParams;
  %twoParams := "";
  
  %twoParams := str paramA & ", " & str paramB;
  %%includeonce% show "the two params are: " & twoParams;
  
  forsuffixes param__ = twoParams:
    
    if isFirstParam__:
      isFirstParam__ := false;
      paramA_str__ := str param__;
      %includeonce% show "paramA_str__ <-- " & paramA_str__;
      
      if not known (scantokens (paramA_str__ & ".className")):
        %includeonce% show "paramA is not an object";
        paramA_isObject__ := false;
      fi;
    else:
      paramB_str__ := str param__;
      %includeonce% show "paramB_str__ <-- " & paramB_str__;
      
      if not known (scantokens (paramB_str__ & ".className")):
        %includeonce% show "paramB is not an object";
        paramB_isObject__ := false;
      fi;
    fi;
  endfor;

  %includeonce% show paramA_str__;
  %includeonce% show paramA_isObject__;
  
  %includeonce% show paramB_str__;
  %includeonce% show paramB_isObject__;
  
  string stringPath__;
  stringPath__ := "(0,0)--(50,50)";
  
  if (not paramA_isObject__) and (not paramB_isObject__):
    %includeonce% show "Both are points!";
    
    
    stringPath__ := paramA_str__ & "--(xpart " & paramB_str__ & ", ypart " & paramA_str__ & ") -- " & paramB_str__;
    %includeonce% show stringPath__;
  else:
    pair refPointA__, refPointB__;
    if paramA_isObject__:
      refPointA__ := scantokens (paramA_str__ & ".c");
    else:
      refPointA__ := scantokens paramA_str__;
    fi;
    if paramB_isObject__:
      refPointB__ := scantokens (paramB_str__ & ".c");
    else:
      refPointB__ := scantokens paramB_str__;
    fi;
    
    %includeonce% show "Reference point A: ";
    %includeonce% show refPointA__;
    %includeonce% show "Reference point B: ";
    %includeonce% show refPointB__;
    
    numeric dx__, dy__;
    dx__ := xpart refPointB__ - xpart refPointA__;
    dy__ := ypart refPointB__ - ypart refPointA__;
    
    string addornmentA__;
    addornmentA__ := "";
    string addornmentB__;
    addornmentB__ := "";
    
    if (paramA_isObject__):
      if (dx__ > 0):
        addornmentA__ := ".e";
      else:
        addornmentA__ := ".w";
      fi;
    fi;
    
    if (paramB_isObject__):
      if (dy__ > 0):
        addornmentB__ := ".s";
      else:
        addornmentB__ := ".n";
      fi;
    fi;
    stringPath__ := "pathManhattanX(" & paramA_str__ & addornmentA__ & ", " & paramB_str__ & addornmentB__ & ")";
    
    %includeonce% show stringPath;
    %includeonce% show "-------------------> recursive call";
  fi;
  
  scantokens stringPath__
enddef;

vardef rpathManhattanX(text twoParams) =
  reverse pathManhattanX(twoParams)
enddef;

vardef pathManhattanY(expr pA,pB) =
  pA--(xpart(pA), ypart(pB))--pB
enddef;

vardef rpathManhattanY(expr pA,pB) =
  pB--(xpart(pA), ypart(pB))--pA
enddef;

vardef pathStepX(expr pA,pB, deltax) =
  pA--(pA+(deltax, 0))--(((xpart pA)+deltax), ypart(pB))--pB
enddef;

vardef rpathStepX(expr pA,pB, deltax) =
  pB--(((xpart pA)+deltax), ypart(pB))--(pA+(deltax, 0))--pA
enddef;

vardef pathStepY(expr pA,pB, deltay) =
  pA--(pA+(0, deltay))--(xpart(pB), ((ypart pA)+deltay))--pB
enddef;

vardef rpathStepY(expr pA,pB, deltay) =
  pB--(xpart(pB), ((ypart pA)+deltay))--(pA+(0, deltay))--pA
enddef;

vardef pathHorizontal(expr pA, untilX) =
  pA--(untilX,ypart(pA))
enddef;

vardef rpathHorizontal(expr pA, untilX) =
  (untilX,ypart(pA))--pA
enddef;

vardef pathVertical(expr pA, untilY) =
  pA--(xpart(pA), untilY)
enddef;

vardef rpathVertical(expr pA, untilY) =
  (xpart(pA), untilY)--pA
enddef;

% Cuts path thePath so that it won't intersect the objects
% objectA and objectB. The border of the objects is obtained
% by invoking the method "_border".
% 
vardef pathCut(suffix objectA, objectB)(expr thePath)=
  save timeA, timeB;
  
  timeA = xpart (thePath intersectiontimes objectBorder(objectA));
  timeB = xpart (thePath intersectiontimes objectBorder(objectB));

  subpath (timeA, timeB) of thePath
enddef;
