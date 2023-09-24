%----------------------------------------------------------------
%       FUZZY SYSTEMS & EVOLUTIONARY COMPUTATION
%       Summer 2020-2021
%       function y=trap_MF(x,a,b,c,d)
%       This function provides the trapezoidal membership function,
%       according to the specifications given in the THIRD PROJECT
%       Author: Prof.Paris Mastorocostas
%----------------------------------------------------------------

function y=trap_MF(x,a,b,c,d)
if a > b
    error('Illegal parameter condition: a > b');
elseif b > c
    error('Illegal parameter condition: b > c');
elseif c > d
    error('Illegal parameter condition: c > d');
end

y=max(min(min(((x-a)/(b-a)),1),((d-x)/(d-c))),0);
