function [rotPoint] = dbkRotPoints(p, Rv, theta)
    % Rotates a set of N points p (3xN)
    % p = [v1, v2, v3, ...vN], vi = [xi; yi; zi]
    % Rv is a vector specifying the axis of rotation, theta is the angle through
    % which to rotate

    if size(p,1) ~= 3
        error('The N points to be transformed must be specified as a 3xN matrix');
    end

    rotPoint = dbkGetRM(Rv, theta)*p;
end
