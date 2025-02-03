function [status] = dbkObjectToPDB(fileName, varargin)
    % Writes one or more objects to a (very minimal) PDB, where each object
    % is a different chain

    % TODO::  Don't write one like at a time!!
    
    f = fopen(fileName, 'w');

    for c=1:(nargin - 1)  % c for chain identifier
        o = varargin{c};
        x = o.x; y = o.y; z = o.z;

        L = {};
        for i=1:length(x)
            L{i} = 'ATOM  '; % next empty col is 7

            % Write atom number
            spacer = repmat(' ', 1, (5 - length(num2str(i))));
            L{i} = [L{i}, spacer, num2str(i), ' ']; % next empty col is 13

            % Write atom type
            spacer = repmat(' ', 1, (4-length(o.atomName{i})));
            L{i} = [L{i}, o.atomName{i}, spacer]; % start in col 13, Next empty col is 17
            
            % Write residue name (cols 18-20)
            if isfield(o, 'resName')
                spacer = repmat(' ', 1, (3-length(o.resName{i})));
                L{i} = [L{i}, ' ', o.resName{i}, spacer];
            else
                L{i} = [L{i}, '    '];
            end % Next empty col is 21

            % Write chain number (assuming only 1 digit)
            L{i} = [L{i}, ' ', num2str(c)]; % Next empty col is 23

            % Write x, y, z coords
            L{i} = [L{i}, '        ']; % Next empty col is 31
            strX = dbkNum2StrFixDecPlaces(o.x(i), 3);
            spacer = repmat(' ', 1, (8 - length(strX)));
            L{i} = [L{i}, spacer, strX]; % Next empty col is 39
            strY = dbkNum2StrFixDecPlaces(o.y(i), 3);
            spacer = repmat(' ', 1, (8 - length(strY)));
            L{i} = [L{i}, spacer, strY];
            strZ = dbkNum2StrFixDecPlaces(o.z(i), 3);
            spacer = repmat(' ', 1, (8 - length(strZ)));
            L{i} = [L{i}, spacer, strZ];
        end
        
        % Write L
        for i=1:length(L)
            fprintf(f, '%s \n', L{i});
        end
        
        fprintf(f, '%s \n', 'TER');
    end  % Go to next object

    fclose(f);

end
