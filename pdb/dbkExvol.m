function [Vex, debug] = dbkExvol(fileName1, fileName2, nP, varargin)

    if ~ischar(fileName1) || ~ischar(fileName2)
        error('PDB filenames must be specified as character arrays.  ');
    end
    if ~isnumeric(nP)
        error('Number of iterations nP must be an integer.  ');
    end

    o1 = dbkPDBToObject(fileName1); o2 = dbkPDBToObject(fileName2);
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
    
    % Create a box around p1 large enough to also fit p2 in any orientation
    % wrt p1 (actually a little larger than necessary:  2*(max(r1) + max(r2));)
    % This is pretty arbirtary :(
    p1Max = max(p1,[],2); p1Min = min(p1,[],2);  % returns: [maxX; maxY; maxZ]
    p2MaxWidth = sqrt(max(max(dSq(p2, p2))));  % scalar
    d = p2MaxWidth + 2*(max(r1) + max(r2));
    box = [(p1Min - d), (p1Max + d)];
    % box = [[-65;-65;-65], [65;65;65]];  % For testing
    boxV = prod(abs(box(:,2) - box(:,1)));
    
    % Generate a set of nP random points ("rP") [[x1;y1;z1], [x2;y2;z2], ...] inside
    % the box
    rP = bsxfun(@plus, bsxfun(@times, rand(3,nP), (box(:,2) - box(:,1))), box(:,1));
    
    if ~isempty(varargin) && strcmp(varargin{1}, 'single')
        r1 = single(r1); r2 = single(r2);
        p1 = single(p1); p2 = single(p2);
        rP = single(rP);
        mdSqm = single(mdSqm);
    end
    
    runTime = tic;
    nF = 0; nS = 0;
    nF1 = 0;  % Number of times only one atom overlaps
    sa = zeros(N1,N2);
    for i=1:nP
        % Rotate p2
        R = dbkGetRM([1;0;0], (2*pi).*rand()) * ...
            dbkGetRM([0;1;0], (2*pi).*rand()) * ...
            dbkGetRM([0;0;1], (2*pi).*rand());
        p2R = R*p2;  % p2R ~ "P2 rotated"
        
        % Translate the rotated p2 into the box
        p2RT = dbkTransPoints(p2R, rP(:,i));  % COM initially at [0;0;0]
        
        % THIS IS DANGEROUS
        %dbkObjectToPDB(strcat(num2str(i), '.pdb'),...
        %    dbkPointsToObject(p1, 'atomRadius', r1, 'atomName', o1.atomName, 'resName', o1.resName),...
        %    dbkPointsToObject(p2RT, 'atomRadius', r2, 'atomName', o2.atomName));
        
        % Test for collision
        d = bsxfun(@gt, dSq(p1, p2RT),  mdSqm); 
        if d % No collision:  all elements of d == 1 (?)
            nS = nS + 1;
        else  % Collision (one or more elements of d == 0 (?))
            nF = nF + 1;
            % Surface area calculation
%             if sum(sum(d)) == ((N1*N2) - 1) % (only one 0 element in d => only one overlapping pair)
%                 sa = sa + (-1.*d + 1);  % Turns 1 -> 0 and 0 -> 1
%                 nF1 = nF1 + 1;   % nF1 => "number of failures w/ only one collision"
%             end
        end
        
        % Provide a time estimate... this does not slow things down at all
        if i == 500
            estRTime = toc(runTime)*(nP/500);
            sTime = now; estFTime = (sTime*86400+estRTime)/86400;  % 86400 sec/day
            fprintf('Estimated run time: %f;  Start time: %s;  Est end time: %s \n',...
                    estRTime, datestr(sTime, 'HH:MM:SS'),...
                    datestr(estFTime, 'HH:MM:SS'))
        end
    end  % Next random orientation and placement (nP total)
    
    % Visualize box
    %o1copy = o1; o1copy.x = p1(1,:)'; o1copy.y = p1(2,:)'; o1copy.z = p1(3,:)';
    %dbkVisualizeSim(rP, o1copy);
    
    % Compute the excluded volume
    Vex = prod(abs(box(:,2) - box(:,1)))*(nF/nP); % boxVol*(frac ! accessible vol)
    
    debug.input_fileName1 = fileName1; debug.input_fileName2 = fileName2;
    debug.input_nP = nP;
    
    debug.Vex = Vex;
%     debug.saNorm1 = (sa./nF1);  % sum(sum(saNorm1)) = 1
%     debug.saNorm2 = sa.*(debug.boxV/nP);
    
    debug.N1 = N1; debug.N2 = N2;
    debug.p1 = p1; debug.p2 = p2;
    debug.r1 = r1; debug.r2 = r2;
    debug.resName1 = o1.resName; debug.resName2 = o2.resName;
    debug.atomName1 = o1.atomName; debug.atomName2 = o2.atomName;
    debug.mdSqm = mdSqm;
    debug.box = box; debug.boxV = boxV;
    % debug.rP = rP;  % Could use a lot of memory if we are calling
    % dbkExvol() multiple times and storing the results
    debug.nF = nF; debug.nS = nS;
    debug.lastDSq = dSq(p1, p2RT);
    debug.lastP2RT = p2RT;
    debug.runTime = toc(runTime);
    debug.sa = sa;
    
end
