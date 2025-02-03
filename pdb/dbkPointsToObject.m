function [o] = dbkPointsToObject(p, varargin)

    if mod((nargin - 1),2) ~= 0  % number of varargin is not even
        error('Number of variable arguments has to be even');
    end
    
    o = struct();
    
    o.x = p(1,:)'; o.y = p(2,:)'; o.z = p(3,:)';

    specifier = {varargin{1:2:(nargin-1)}};
    payload = {varargin{2:2:(nargin-1)}};
    % specifier(i) corresponds to payload{i}
    % Note both are cell arrays

    [idx, numOccur] = dbkStrFind(specifier, 'atomRadius');
    if numOccur == 1
        o.atomRadius = payload{idx};
    elseif numOccur > 1
        error('Can''t specify atomRadius field more than once');
    end
    
    [idx, numOccur] = dbkStrFind(specifier, 'atomName');
    if numOccur == 1
        o.atomName = payload{idx};
    elseif numOccur > 1
        error('Can''t specify atomName field more than once');
    end
    
    [idx, numOccur] = dbkStrFind(specifier, 'atomSerial');
    if numOccur == 1
        o.atomSerial = payload{idx};
    elseif numOccur > 1
        error('Can''t specify atomSerial field more than once');
    end
    
    [idx, numOccur] = dbkStrFind(specifier, 'resName');
    if numOccur == 1
        o.resName = payload{idx};
    elseif numOccur > 1
        error('Can''t specify resName field more than once');
    end
    
    %o

 % OLD %
%     o.x = p(1,:)'; o.y = p(2,:)'; o.z = p(3,:)';
%     
%     if ~isempty(varargin{1})
%         o.r = varargin{1};
%     end
%     if ~isempty(varargin{2})
%         o.atom = varargin{2};
%     end
% 

end
