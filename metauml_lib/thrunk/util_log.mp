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

%%% if known _util_log_mp:
%%%   expandafter endinput
%%% fi;
_util_log_mp:=1;

util_loglevel_PARANOIA:=1;
util_loglevel_DEBUG:=2;
util_loglevel_RELEASE:=3;
util_loglevel_SILENT:=4;

util_log_thresholdlevel := util_loglevel_DEBUG;
util_log_defaultlevel   := util_loglevel_RELEASE;

def elog(expr loglevel) text txt=
  if (loglevel <= util_log_thresholdlevel):
    %includeonce% show txt
  else:
    % do nothing
  fi;
enddef;

def log text txt=
  if (util_log_defaultlevel <= util_log_thresholdlevel):
    %includeonce% show txt
  else:
    % do nothing
  fi;
enddef;
