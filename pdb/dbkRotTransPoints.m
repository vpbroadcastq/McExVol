function [rtP, RT] = dbkRotTransPoints(p, Tv, Rv, theta)
    % Rotates AND THEN translates a set of N points p (3xN)
    % p = [v1, v2, v3, ...vN], vi = [xi; yi; zi]
    % Tv is a vector specifying the translation (the displacement, not the
    % point to which we want to translate)
    % Rv is a vector specifying the axis of rotation, theta is the angle through
    % which to rotate

    if size(p,1) ~= 3
        error('The N points to be transformed must be specified as a 3xN matrix');
    end
    if size(Tv,1) ~= 3 || size(Tv,2) ~= 1
        error('The translation vector must be 3x1');
    end
    if size(Rv,1) ~= 3 || size(Rv,2) ~= 1
        error('The rotation vector must be 3x1');
    end
    if size(theta,1) ~= 1 || size(theta,2) ~= 1
        error('The rotation angle theta must be a scalar');
    end
    
    if norm(Rv) ~= 0  % Can't norm the vector if it's [0;0;0] => divide by 0
        Rv = Rv/norm(Rv);
    end
    Rx = Rv(1); Ry = Rv(2); Rz = Rv(3);

    c = cos(theta);
    s = sin(theta);
    t = 1 - cos(theta);

    % Build matrix for rotation about Rv, translation by Tv
    % Verified the rotation part is the same as is used in the MSL source
    % CartesianGeometry::getRotationMatrix()
    % RT ~ "rotation translation"
    RT = [t*Rx^2 + c, t*Rx*Ry - s*Rz, t*Rx*Rz + s*Ry, Tv(1)];
    RT = [RT; t*Rx*Ry + s*Rz, t*Ry^2 + c, t*Ry*Rz - s*Rx, Tv(2)];
    RT = [RT; t*Rx*Rz - s*Ry, t*Ry*Rz + s*Rx, t*Rz^2 + c, Tv(3)];
    RT = [RT; 0, 0, 0, 1];

    rtP = RT*[p; ones(1,size(p,2))];  % rtP ~ "rotated translated point"
    rtP = rtP(1:3,:);
end
