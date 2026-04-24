%% y = steppoint(x,a,b,c[,d])
%%
%% function to fit a two-straight-line model to data where a step change is expected before or after a flat plateau
%%
%% INPUTS
%% x = single column of data
%% a = step point of curve
%% b = baseline
%% c = plateau
%% d = direction of step: -1=before, 1=after the plateau; default=1
%%
%% example:
%% mfit = fit(x', y, 'steppoint(x,a,b,c,1)', 'StartPoint', [a,b,c]);
%%
function y = steppoint(x,a,b,c,d)
    if nargin==4 || isempty(d)                                              % if no direction argument is given
        d = 1;
    end
    if abs(d)~=1                                                            % if non-unit direction, convert to unit, preserving sign
        d = d./abs(d);
    end
    if d==1
        y = b + double(x>=a).*c;                                            % step (up) occurs before the plateau
    elseif d==-1
        y = b + double(x<=a).*c;                                            % step (down) occurs after the plateau
    end
end