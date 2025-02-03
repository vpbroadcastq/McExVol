function [idx, numOccur] = dbkStrFind(C, s)
    % Returns a col vector containing the indicies of the elements of C that match
    % s (idx) and and the number of occurances of s in C (numOccur).  
    % s is the string to find
    % C is a cell array
    
    x = strfind(C, s);  % x is a cell array
    numOccur = 0;
    idx = [];
    for i=1:numel(x)
        if ~isempty(x{i})
            numOccur = numOccur + 1;
            idx = i;
        end
    end

end
