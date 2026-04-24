%% y = ramppoint(x,a,b,c[,d])
%%
%% function to fit a three-straight-line model to data where a linear change is expected between a flat baseline and flat plateau
%%
%% INPUTS
%% x = single column of data
%% a = start point of ramp
%% b = baseline
%% c = end point of ramp
%% d = plateau
%% e = direction of slope: 1 = positive, -1 = negative
%%
%% example:
%% mfit = fit(x', y, 'ramppoint(x,a,b,c,1)', 'StartPoint', [a,b,c]);
%%
function y = ramppoint(x,a,b,c,d,e)
    if nargin==5 || isempty(e)                                       % if no direction argument is given
        e = 1;
    end
    if abs(e) ~= 1                                                   % if non-unit direction, convert to unit, preserving sign
        e = e./abs(e);
    end
    %if e==1
        y = b + double(x>=a & x<c).*(x-a).*((d-b)./(c-a)) + double(x>=c).*d;% baseline + x_diff*slope + plateau
    %elseif e==-1
        %y = b + double(x<a).*((x-a).*c);                             % slope occurs before the plateau
    %end
end