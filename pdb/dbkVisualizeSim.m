function dbkVisualizeSim(rP, o1)
    
    % Put a water molecule at every random point within the box
    aSerial = 1:size(rP,2);
    aName = cellstr(repmat('O', [size(rP,2),1]));
    resName = cellstr(repmat('HOH', [size(rP,2),1]));
    
    oRP = dbkPointsToObject(rP, 'atomName', aName, 'atomSerial', aSerial, 'resName', resName);
    
    dbkObjectToPDB('visualPoints.pdb', oRP, o1);
end