
p1R = [-0.1640362;-0.4328932;2.701168;]  % Fixed at 0,0,0
p2RT = [-1.2397; -0.5879; 15.7896]  % Moves
rTraj = [0.1260; -0.5739; -0.8092]
mdSqm = (1.5 + 2)^2

a = (rTraj)'*(rTraj)
b = 2*(rTraj')*(p2RT + p1R)
c = (p2RT)'*(p2RT) + (p1R)'*(p1R) - 2*(p2RT')*p1R - mdSqm

%a = (rTraj.^2)'*(p2RT.^2)
%b = -2*sum(rTraj.*p1R.*p2RT)
%c = sum(p1R.^2) - mdSqm

sigma = [(1/(2*a))*(-b + sqrt(b^2 - 4*a*c)); (1/(2*a))*(-b - sqrt(b^2 - 4*a*c))]

if sigma > 0 && isreal(sigma)
    % collision is possible
    % how many steps?
    % sigma/stepSize
    % Advance this many steps
    'yes'
else
    % no collision is possible
    'no'
end

