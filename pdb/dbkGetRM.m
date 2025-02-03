function R = dbkGetRM(Rv, theta)
    % Rv is a vector specifying the axis of rotation, theta is the angle through
    % which to rotate
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

    % Build matrix for rotation about Rv
    % Verified this is the same as is used in the MSL source
    % CartesianGeometry::getRotationMatrix()
    % RT ~ "rotation translation"
    R = [t*Rx^2 + c, t*Rx*Ry - s*Rz, t*Rx*Rz + s*Ry];
    R = [R; t*Rx*Ry + s*Rz, t*Ry^2 + c, t*Ry*Rz - s*Rx];
    R = [R; t*Rx*Rz - s*Ry, t*Ry*Rz + s*Rx, t*Rz^2 + c];
    
end