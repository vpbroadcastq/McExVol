function dS = dSq(a, b, varargin)
    % Returns the distance squared between each pair of points in a set of
    % points a and a set of points b:
    % a = [(x1; y1; z1), (x2; y2; z2), ...]
    % It is possible to specify single precision like:
    % dSq(a, b, 'single')
    
    if (nargin < 2)
       error('Not enough input arguments');
    end
    
    if (size(a,1) ~= size(b,1))
       error('A and B should be of the same dimensionality');
    end
    
    if (length(varargin) == 1) && strcmp(varargin{1}, 'single')  % Run in single precision
        a = single(a); b=single(b);
        aa=single(sum(a.*a,1)); bb=single(sum(b.*b,1)); ab=single(a'*b);
    else  % Run in double precision
        aa=sum(a.*a,1); bb=sum(b.*b,1); ab=a'*b;
    end
    
    dS = abs(bsxfun(@plus, aa', bb) - 2*a'*b);  % dS ~ "distance squared"
end
