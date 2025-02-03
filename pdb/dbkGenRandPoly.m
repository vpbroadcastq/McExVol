function [d] = dbkGenRandPoly(nMono, r, L, nPoly)
    % nMono:  # of monomers in the chain
    % r:  radius of each monomer (spherical)
    % L:  monomer-monomer bond length
    % nPoly:  # of polymers to generate

    cPoly = 1;
    while cPoly <= nPoly
        theta = 2*pi*rand(nMono-2);
        phi = 2*pi*rand(nMono-2);
        
        cY = zeros(nMono,1); cX = zeros(nMono,1); cZ = zeros(nMono,1);
        cY(1,1) = 0; cX(1,1) = 0; cZ(1,1) = 0;
        cY(2,1) = 0; cX(2,1) = 0; cZ(2,1) = L;

        for i=3:nMono
            cY(i,1) = L*sin(theta(i-2,1))*sin(phi(i-2,1)) + cY(i-1,1);
            cX(i,1) = L*cos(phi(i-2,1))*sin(theta(i-2,1)) + cX(i-1,1);
            cZ(i,1) = L*cos(theta(i-2,1)) + cZ(i-1,1);
        end

        d = dSq([cX, cY, cZ]', [cX, cY, cZ]');

        % d is nMono*nMono
        % entries on the main diagonal should be 0
        
        % replace entries on the main diagonal of d with values >r
        % we don't care if an atom overlaps itself (2*r just to be sure)
        d = d + diag(((2*r)^2)*ones(nMono,1), 0);

        %if bsxfun(@lt, d, r)
        if bsxfun(@ge, d, (2*r)^2)
            % no collisions
            
            dbkObjectToPDB(strcat(num2str(cPoly), '.pdb'),...
                dbkPointsToObject([cX, cY, cZ]', 'atomRadius', r.*ones(nMono,1),...
                'atomName', repmat({'C'}, nMono, 1), 'resName',...
                repmat({'peg'}, nMono, 1)));

            cPoly = cPoly + 1;
        else
            % collision
        end

    end

end
% y[j] = (float)(L*sin((double)theta[j])*sin((double)phi[j]) + (double)y[j-1]);
% x[j] = (float)(L*cos((double)phi[j])*sin((double)theta[j]) + (double)x[j-1]);
% z[j] = (float)(L*cos((double)theta[j]) + (double)z[j-1]);