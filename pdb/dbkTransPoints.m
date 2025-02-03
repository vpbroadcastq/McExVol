function [transPoint] = dbkTransPoints(p, Tv)
    % Translates a set of N points p (3xN)
    % p = [v1, v2, v3, ...vN], vi = [xi; yi; zi]
    % Tv is a vector specifying the translation (the displacement, not the
    % point to which we want to translate)

    if size(p,1) ~= 3
        error('The N points to be transformed must be specified as a 3xN matrix');
    end
    if size(Tv,1) ~= 3 || size(Tv,2) ~= 1
        error('The translation vector must be 3x1');
    end
    
    transPoint = bsxfun(@plus, p, Tv);
end
