%% y = kneepoint(x,a,b,c[,d])
%%
%% function to fit a two-straight-line model to data where a linear change is expected before or after a flat plateau
%%
%% INPUTS
%% x = single column of data
%% a = kneepoint
%% b = plateau
%% c = slope
%% d = direction of slope: -1=before, 1=after the plateau; default=1
%%
%% example:
%% mfit = fit(x', y, 'kneepoint(x,a,b,c,1)', 'StartPoint', [a,b,c]);
%%
function y = kneepoint(x,a,b,c,d)
    if nargin==4 || isempty(d)                                       % if no direction argument is given
        d = 1;
    end
    if abs(d) ~= 1                                                   % if non-unit direction, convert to unit, preserving sign
        d = d./abs(d);
    end
    if d==1
        y = b + double(x>=a).*((x-a).*c);                            % slope occurs after the plateau
    elseif d==-1
        y = b + double(x<=a).*((x-a).*c);                            % slope occurs before the plateau
    end
end