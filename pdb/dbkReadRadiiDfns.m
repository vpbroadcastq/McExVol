function [definedAtoms] = dbkReadRadiiDfns(fileName)

    % Read in the radius definition file omitting the comments and blank
    % lines (???)
    f = fopen(fileName, 'r');
    cL = fgetl(f); i = 1;  % cL ~ "current line"
    while ischar(cL)
        if ~strcmp(cL(1,1), '%') && ~isempty(cL)
            L{i,1} = cL;
            i = i + 1; % Only increment i if we want to read the line
        end
        cL = fgetl(f);
    end
    fclose(f);

    % Parse the atom type definitions (specifies which radii to assign to
    % which atom types, and defines all the allowable atom types)
    if strcmp(strtrim(L{1,1}(1,:)), '[SECTION:  RADIUS DEFINITION]');
        i = 2;  % Don't want the [SECTION... line
        j = 1;
        while ~strcmp(strtrim(L{i,1}(1,:)), '[SECTION:  ATOM TYPE]')
            [definedAtomType{j,1}(1,:), currRadius] = strtok(L{i,1}(1,:), ',');
            
            definedTypeRadii(j,1) = str2double(currRadius(1,2:end));  % Hack off the leading comma
            
            i = i + 1;
            j = j + 1;
        end
    else
        error('Could not find radius definition section');
    end

    % Parse the atom type table
    if strcmp(strtrim(L{i,1}(1,:)), '[SECTION:  ATOM TYPE]');  % i is assigned in the previous step
        i = i + 1;
        j = 1;
        while i <= size(L, 1)
            [definedResName{j,1}(1,:), remain] = strtok(L{i,1}(1,:), ',');
            
            remain = remain(1,2:end);  % Hack off the leading comma
            [definedAtomName{j,1}(1,:), currAtomType] = strtok(remain, ',');
            
            atomType{j,1}(1,:) = strtrim(currAtomType(1,2:end));  % Hack off the leading comma
            
            i = i + 1;
            j = j + 1;
        end
    else
        error('Could not find atom type section');
    end

    % Copy the proper radius to each atom type in the rDfn struct
    for i=1:size(definedAtomName, 1)  % CA, N, O, CG, O4', etc
        for j=1:size(definedAtomType, 1)  % sp3c, sp2c, etc
            if strcmp(strtrim(atomType{i,1}(1,:)), strtrim(definedAtomType{j,1}(1,:)))
                r(i,1) = definedTypeRadii(j,1);
                break;
            end

            if j == size(definedAtomType,1)
                % No match was found
                warning('No match was found for atom name %d, %s %s (%s)', ...
                        i, definedAtomName{i,1}(1,:), strtrim(atomType{i,1}(1,:)), ...
                        definedResName{i,1}(1,:));
            end
        end
    end
    definedAtoms = struct();
    definedAtoms.resName = definedResName;
    definedAtoms.atomName = definedAtomName;
    definedAtoms.atomType = atomType;
    definedAtoms.atomRadius = r;

end
