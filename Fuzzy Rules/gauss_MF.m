%----------------------------------------------------------------
%       FUZZY AND UNCERTAINTY SYSTEMS
%       function y=gauss_MF(x,m,sigma)
%       This function provides the gaussian membership function
%----------------------------------------------------------------

function y=gauss_MF(x,m,sigma)
if sigma == 0
    error('The Gaussian MF needs a non-zero sigma.');
end

y = exp(-(x - m).^2/(2*sigma^2));
