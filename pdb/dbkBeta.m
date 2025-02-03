function [debug] = dbkBeta(fileName1, fileName2, nRuns, stepSize, varargin)

    if ~ischar(fileName1) || ~ischar(fileName2)
        error('PDB filenames must be specified as character arrays.  ');
    end
    if ~isnumeric(nRuns)
        error('Number of iterations nP must be an integer.  ');
    end

    o1 = dbkPDBToObject(fileName1); o2 = dbkPDBToObject(fileName2);
    atomName1 = o1.atomName; atomName2 = o2.atomName;
    r1 = o1.atomRadius; r2 = o2.atomRadius;  % [r1;r1;r3;...]
    p1 = dbkObjectToPoints(o1); p2 = dbkObjectToPoints(o2);  % [[x1;y1;z1],[x2;y2;z2],...]
    N1 = size(p1,2); N2 = size(p2,2);
    
    %  Zero the center of mass of p1 and p2
    cm_p1 = mean(p1,2); cm_p2 = mean(p2,2);  % Center of mass calc
    p1 = dbkTransPoints(p1, -1*cm_p1);
    p2 = dbkTransPoints(p2, -1*cm_p2);
    
    % Calculate the minimum distance squared matrix (mdSqm)
    mdSqm = zeros(N1, N2);
    for i=1:N2
        mdSqm(:,i) = (r2(i,1) + r1(:,1)).^2;  % mdSqm(atom i p1, atom j p2)
    end
    
    % Calculate the maximum width of molecule 1 and molecule 2
    p1MaxWidth = sqrt(max(max(dSq(p1, p1))));  % scalar
    p2MaxWidth = sqrt(max(max(dSq(p2, p2))));  % scalar
    %startDistance = p1MaxWidth + p2MaxWidth + 2*(max(r1) + max(r2)) + 9;
    % startDistance is the starting distance between the CENTERS OF MASS
    
    if ~isempty(varargin) && strcmp(varargin{1}, 'single')
        r1 = single(r1); r2 = single(r2);
        p1 = single(p1); p2 = single(p2);
        mdSqm = single(mdSqm);
    end
    
    runTime = tic;
    numctot = 0; numc1 = 0; numc2 = 0; numc3 = 0; numc4 = 0; noncounted = 0;
    % numctot = numc1 + 2*numc2 + 2*numc3 + 4*numc4         NOT + noncounted
    % "noncounted" would be more accurately named "notclassified" or
    % "notrecorded"
    ctot = zeros(N1,N2);  % All types of _counted_ collisions, including
                          % multi-atom collisions, added together
    c1 = zeros(N1,N2);  % Number of times one atom on molecule 1 overlaps 1 atom on molecule 2

    if N1 > 1
        c2 = zeros((1/2)*(N1^2 - N1),N2); % number of times two atoms on molecule 1 overlap 1 atom on molecule 2
        % Build the c2 index matrix c2Index.  Associates row numbers of the c2
        % matrix with different possible collision types
        c2Index(:,1) = (1:((1/2)*(N1^2 - N1)))';
        currRow = 0;
        for i=1:N1
            for j = 1:N1
                if i == j
                    continue;
                elseif i > j
                    continue;
                elseif i < j
                    currRow = currRow + 1;
                    c2Index(currRow,2:3) = [i, j];
                end
            end
        end
        %c2Index
    end
    
    if N2 > 1
        c3 = zeros(N1,(1/2)*(N2^2 - N2)); % number of times two atoms on molecule 2 overlap 1 atom on molecule 1
        % Build the c3 index matrix c3Index.  Associates column numbers of the
        % c3 matrix with different possible collision types
        c3Index(1,:) = (1:((1/2)*(N2^2 - N2)))';
        currCol = 0;
        for i=1:N2
            for j = 1:N2
                if i == j
                    continue;
                elseif i > j
                    continue;
                elseif i < j
                    currCol = currCol + 1;
                    c3Index(2:3,currCol) = [i; j];
                end
            end
        end
        %c3Index
    end
    
    if N1 > 1 && N2 > 1
        c4 = 0; % number of times two atoms on molecule 1 overlap 2 atoms on molecule 2
        % NOT YET PROPERLY IMPLEMENTED
    end
    
    % Generate nRuns "random" trajectories by sampling randomly a circle of radius
    % "sampleRadius" (centered at x=0, y=0) chosen to be wide enough such that if
    % the center of mass of molecule 2 is placed _outside_ the circle there is no way
    % any atoms of molecule 2 can collide with molecule 1 (located at the
    % center of the circle).  
    % Molecule 2 starts at (0,0,startDistance) where startDistance is some
    % distance much larger than the radius of the circle (chosen here to be
    % 50*sampleRadius, which is arbitrary)
    sampleRadius = p1MaxWidth + p2MaxWidth + 2*(max(r1) + max(r2));
    startDistance = 50*sampleRadius;
    rTraj = zeros(3,nRuns);
    r = sampleRadius.*sqrt(rand(nRuns, 1)); theta = 2*pi.*rand(nRuns, 1);
    rTraj(1,:) = (r.*cos(theta))'; rTraj(2,:) = (r.*sin(theta))';
    rTraj(3,:) = -1*startDistance.*ones(1,nRuns);
    % ie, the coordinates of the point i want to hit are
    % [x;y;-startDistance] if the origin is at the COM of molecule 2
    % (molecule 2 moves)
    for i=1:nRuns
        rTraj(:,i) = rTraj(:,i)./norm(rTraj(:,i));
    end
    
    for i=1:nRuns
        % Rotate p2
        R = dbkGetRM([1;0;0], (2*pi).*rand()) * ...
            dbkGetRM([0;1;0], (2*pi).*rand()) * ...
            dbkGetRM([0;0;1], (2*pi).*rand());
        p2R = R*p2;  % p2R ~ "P2 rotated"
        
        % Rotate p1
        R = dbkGetRM([1;0;0], (2*pi).*rand()) * ...
            dbkGetRM([0;1;0], (2*pi).*rand()) * ...
            dbkGetRM([0;0;1], (2*pi).*rand());
        p1R = R*p1;  % p2R ~ "P2 rotated"
        
        % Translate the rotated p2 (by its center of mass) up the z-axis a
        % distance startDistance
        p2RT = bsxfun(@plus, p2R, [0;0;startDistance]);  % COM initially at [0;0;0]
        
        % Calculate which atoms are going to be the ones to collide
        sigma = zeros(N1,N2);  % multiplier applied to rTraj
        collide = 0;  % Is there at least one collision?
        for j=1:N2
            for k=1:N1
                % sigma is the multiplier of the rTraj vector that solves
                % the eqn below for atoms j,k
                % mdSq = a1*sigma^2 + a2*sigma + a3
                a1 = (rTraj(:,i))'*(rTraj(:,i));
                a2 = 2*(rTraj(:,i)')*(p2RT(:,j) - p1R(:,k));
                a3 = (p2RT(:,j))'*(p2RT(:,j)) + (p1R(:,k))'*(p1R(:,k)) - 2*(p2RT(:,j)')*p1R(:,k) - mdSqm(k,j);

                sPlus = (1/(2*a1))*(-a2 + sqrt(a2^2 - 4*a1*a3));
                sMinus = (1/(2*a1))*(-a2 - sqrt(a2^2 - 4*a1*a3));
                sPlusMinus = [sPlus; sMinus];
                
                % Get rid of imaginary and negative solutions
                sPlusMinus = sPlusMinus(~imag(sPlusMinus));
                sPlusMinus = sPlusMinus(sPlusMinus > 0);
                sPlusMinus = min(sPlusMinus);
                
                % Detect and store collisions
                % sigma(j,k) contains the solution to the quadratic that
                % causes a collision between atoms j,k
                % Small # => early collision
                % sigma(k,j) = 0 => no collision
                if isempty(sPlusMinus)
                    sigma(k,j) = 0;
                else
                    collide = 1;  % At least one collision sets collide
                    sigma(k,j) = sPlusMinus;
                end
            end
        end
        %sigma
        % Store atomic collisions in the "c" matrix
        if collide == 1  % There was at least one collision
            firstToCollide = min(min(sigma(sigma ~= 0)));
            % Want to record not just the very first two atoms to collide, but
            % also atoms on a collision path that are quite close to also
            % colliding
            c = sigma.*(sigma <= (firstToCollide + stepSize)).*(sigma >= (firstToCollide - stepSize)).*(sigma ~= 0);
            % molecule 1 atoms form the rows, molecule 2 atoms form the
            % cols of the c matrix
            c = ones(size(c,1), size(c,2)).*(c>0);
        end
        
        % Classify and store detailed information about the collision
        if collide == 1  % There was at least one collision
            if sum(sum(c)) == 1  % 1 - 1
                c1 = c1 + c;
                numc1 = numc1 + 1;
                ctot = ctot + c;
                numctot = numctot + 1;
            elseif sum(sum(c)) == 2 && max(sum(c,1)) == 2  % type 2
                % one atom on molecule 2 is colliding with 2 atoms on
                % molecule 1
                %c
                ctot = ctot + c;  % overall collision matrix
                numc2 = numc2 + 1;
                numctot = numctot + 2;  % Note!!  Counts as two!

                % find row and col indices of nonzero elements of c
                [rIDX, cIDX] = ind2sub([size(c,1), size(c,2)], find(c));
                % 2 unique rows but only one unique col
                % use the c2 index to determine the proper row number of the c2 matrix
                % in which to store the collision.  The proper col of the
                % c2 matrix is the same as the col of the c matrix cIDX
                c2RowIDX = c2Index(:,1).*(c2Index(:,2) == rIDX(1)).*(c2Index(:,3) == rIDX(2));
                c2RowIDX = c2RowIDX(c2RowIDX ~= 0);  % Should be only one entry
                c2ColIDX = cIDX(1,1);

                c2(c2RowIDX,c2ColIDX) = c2(c2RowIDX,c2ColIDX) + 1;
            elseif sum(sum(c)) == 2 && max(sum(c,2)) == 2  % type 3
                % one atom on molecule 1 is colliding with 2 atoms on
                % molecule 2
                %sigma
                %c
                ctot = ctot + c;  % overall collision matrix
                numc3 = numc3 + 1;
                numctot = numctot + 2;  % Note!!  Counts as two!

                % find row and col indices of nonzero elements of c
                [rIDX, cIDX] = ind2sub([size(c,1), size(c,2)], find(c));
                % 2 unique cols but only one unique row
                % use the c3 index to determine the proper col number of the c3 matrix
                % in which to store the collision.  The proper row of the
                % c3 matrix is the same as the row of the c matrix rIDX
                c3ColIDX = c3Index(1,:).*(c3Index(2,:) == cIDX(1)).*(c3Index(3,:) == cIDX(2));
                c3ColIDX = c3ColIDX(c3ColIDX ~= 0);  % Should be only one entry
                c3RowIDX = rIDX(1,1);
                
                c3(c3RowIDX,c3ColIDX) = c3(c3RowIDX,c3ColIDX) + 1;
            elseif sum(sum(c)) == 400 % Type 4 collision... how to define???
                c
                c4 = c4 + c;
                numc4 = numc4 + 1;  % Note!!  Counts as four!
                ctot = ctot + c;
                numctot = numctot + 4;
            else  % some sort of many-atom collision
                noncounted = noncounted + 1;
            end
        end  % Finished classifying & recording the collision
        
        % Provide a time estimate... this does not slow things down at all
        if i == 1000
            estRTime = toc(runTime)*(nRuns/1000);
            sTime = now; estFTime = (sTime*86400+estRTime)/86400;  % 86400 sec/day
            fprintf('Estimated run time: %f;  Start time: %s;  Est end time: %s \n',...
                    estRTime, datestr(sTime, 'HH:MM:SS'),...
                    datestr(estFTime, 'HH:MM:SS'))
        end
    end  % Next random orientation and run (nRuns total)
    
    % End of simulation
    % Compile all the data (sum up over different atom types)
    % (atom types are close enough to surface types it's not too much of a
    % pain in the ass to do the final bit of compilation in Excel)
    
    % c1 and ctot collisions
    uAtom1 = unique(atomName1); uAtom2 = unique(atomName2);  % Cell arrays
    % Add surface types on molecule 1 (rows of c1)
    c1ByAtomType_Rows = [];  ctotByAtomType_Rows = [];
    for i=1:length(uAtom1)
        currRows = strcmp(atomName1, uAtom1(i,1));
        
        c1ByAtomType_Rows = [c1ByAtomType_Rows; sum(c1(currRows,:),1)];
        ctotByAtomType_Rows = [ctotByAtomType_Rows; sum(ctot(currRows,:),1)];
    end
    % Now add the cols
    c1ByAtomType_RowsCols = [];  ctotByAtomType_RowsCols = [];
    for i=1:length(uAtom2)
        currCols = strcmp(atomName2, uAtom2(i,1));  %currCols = currCols';
        
        c1ByAtomType_RowsCols = [c1ByAtomType_RowsCols, sum(c1ByAtomType_Rows(:,currCols),2)];
        ctotByAtomType_RowsCols = [ctotByAtomType_RowsCols, sum(ctotByAtomType_Rows(:,currCols),2)];
    end
    %c1ByAtomType_RowsCols
    %ctotByAtomType_RowsCols
    
    beta1 = c1ByAtomType_RowsCols.*(1/nRuns).*(pi*sampleRadius^2);
    betatot = (1/numctot).*ctotByAtomType_RowsCols;
    
    % Function parameters
    debug.input_fileName1 = fileName1; debug.input_fileName2 = fileName2;
    debug.input_nRuns = nRuns; debug.input_stepSize = stepSize;
    
    % Molecule 1, molecule 2 info
    debug.N1 = N1; debug.N2 = N2;
    debug.p1 = p1; debug.p2 = p2;
    debug.r1 = r1; debug.r2 = r2;
    debug.resName1 = o1.resName; debug.resName2 = o2.resName;
    debug.atomName1 = atomName1; debug.atomName2 = atomName2;
    debug.p1MaxWidth = p1MaxWidth; debug.p2MaxWidth = p2MaxWidth;
    %debug.c2Index = c2Index; debug.c3Index = c3Index;
    
    % Info about the simulation
    debug.mdSqm = mdSqm; debug.startDistance = startDistance;
    debug.sampleRadius = sampleRadius;
    %debug.trajectories = rTraj;  % Could eat up a _lot_ of space
    debug.lastp1R = p1R; debug.lastP2RT = p2RT;
    debug.runTime = toc(runTime);
    debug.lastSigma = sigma;

    % Collision data - uncompiled
    debug.ctot = ctot; debug.c1 = c1; debug.c2 = c2; debug.c3 = c3; debug.c4 = c4;
    debug.c2Index = c2Index;
    debug.numctot = numctot;
    debug.numc1 = numc1; debug.numc2 = numc2; debug.numc3 = numc3; debug.numc4 = numc4;
    debug.noncounted = noncounted;
    
    % Collision data - Compiled by atom type
    debug.beta1RowNames = uAtom1;  debug.beta1ColNames = uAtom2;
    debug.beta1 = beta1;
    debug.betatot = betatot;  % Row and col names same as for beta1

end
