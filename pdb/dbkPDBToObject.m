function [o] = dbkPDBToObject(fileName)
    % No-mercy PDB parser
    % There may be only one structure per file (multiple disconnected chains are OK)
    
    % Read in only the ATOM and HETATOM lines into cell array L
    f = fopen(fileName, 'r');
    cL = fgetl(f); i = 1;  numTER = 0;% cL ~ "current line"
    while ischar(cL)
        
        if strcmp(cL(1,1:3), 'END')
            break;
        end
        
        if strcmp(cL(1,1:4), 'ATOM') || strcmp(cL(1,1:6), 'HETATM')
            L{i,1} = cL;
            i = i + 1; % Only increment i if we want to read the line
        end
        
        % Count the number of TER lines
        if strcmp(cL(1,1:3), 'TER')
            % The TER line is often given an atom number (for reasons i
            % don't understand), which messes with one of my sanity checks
            % below.  
            numTER = numTER + 1;
        end
        cL = fgetl(f);
    end
    fclose(f);

    o = struct();
    
    % Read in using the col spacing definitions from the PDB standard.
    % This is the strictest method of pdb reading.  
    pProb = [];
    for i=1:size(L,1)
        % Atom serial number (cols 7-11)
        o.atomSerial(i,1) = str2double(L{i}(1,7:11));
        if isnan(o.atomSerial(i,1))
            error('Atom serial number %d is empty', i);
        end
        
        % Atom name (13-16) and residue name (18-20) (used for radius assignment)
        o.atomName{i,1} = strtrim(L{i}(13:16));
        if isempty(o.atomName{i,1})
            error('Atom name %d is empty', i);
        end
        o.resName{i,1} = strtrim(L{i}(18:20));
        if isempty(o.resName{i,1})
            error('Atom resName %d is empty', i);
        end
        
        % X (31-38), Y (39-46), Z (47-54) coordinates
        o.x(i,1) = str2double(strtrim(L{i}(1,31:38)));
        o.y(i,1) = str2double(strtrim(L{i}(1,39:46)));
        o.z(i,1) = str2double(strtrim(L{i}(1,47:54)));
        if isnan(o.x(i,:)) || isnan(o.y(i,:)) || isnan(o.z(i,:))
            error('Unable to read X or Y or Z coordinate for atom %d', i);
        end
    end
    
    % Assign radii to atoms in o
    % Want to match BOTH the atom name and the residue names to the set of
    % defined atom types from dbkReadRadiiDfns().  The most straightforward
    % way to do this is to concatenate the atomName and resName strings and
    % compare these mashed-together values rather than attempting two
    % separate comparisons.  
    definedAtoms = dbkReadRadiiDfns('radii.csv');
    defRNameAName = strcat(definedAtoms.resName, definedAtoms.atomName);  % concat resNames
                                                                          % and atomNames with
                                                                          % defined radii
    fileRNameAName = strcat(o.resName, o.atomName);  % concat resNames and atomNames from PDB
    for i=1:size(fileRNameAName,1)
        matchLoc = find(strcmp(strtrim(fileRNameAName(i)), defRNameAName));
        if numel(matchLoc) == 1
            o.atomRadius(i,1) = definedAtoms.atomRadius(matchLoc,1);
        elseif numel(matchLoc) > 1
            error('More than one match found for atom %d (%s)', i, fileRNameAName{i,1});
        elseif numel(matchLoc) < 1
            error('No match found for atom %d (%s)', i, fileRNameAName{i,1});
        end
    end
    
    % Some sanity checks
    % Make sure the atomSerial field goes from 1 -> max w/o skipping anything
    % Note the + numTER... for some reason the TER line is given an atom
    % serial number
    if (length(o.atomSerial) + numTER) ~= max(o.atomSerial)
        error('Highest serial number %i != to number of serial numbers %i',...
              (length(o.atomSerial) + numTER), max(o.atomSerial));
    end
    
    % Make sure there are not multiple atoms occupying the same point
    if length(unique([o.x, o.y, o.z], 'rows')) ~= length([o.x, o.y, o.z])
        error('Multiple atoms occupying the same point');
    end
    
end
