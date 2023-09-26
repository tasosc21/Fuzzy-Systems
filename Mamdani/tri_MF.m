%----------------------------------------------------------------
%       FUZZY SYSTEMS & EVOLUTIONARY COMPUTATION
%       function y=tri_MF(x,a,b,c)
%       This function provides the triangular membership function
%----------------------------------------------------------------

function y=tri_MF(x,a,b,c)
if a > b,
    error('Illegal parameter condition: a > b');
elseif b > c,
    error('Illegal parameter condition: b > c');
end

y=max( min(((x-a)/(b-a)),((c-x)/(c-b))),0 );


