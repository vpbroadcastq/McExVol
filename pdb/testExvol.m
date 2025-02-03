% dbkExvol simulation test script

% Reset the RNG... note the syntax for doing this has changed in the
% most recent version
s = RandStream.getDefaultStream;
reset(s);

%
% Test set 1:  PDB-PDB
%
% structureSet = {'2L63.pdb', '13-BD.pdb';
%               '13-BD.pdb', '2L63.pdb'};
%
% Test set 2:  Simple shape - Simple shape
%
structureSet = {'sphere1.pdb', 'sphere2.pdb';
                       'cyl1.pdb', 'sphere1.pdb';
                       'sphere1.pdb', 'cyl1.pdb';
                       '2L63.pdb', '13-BD.pdb';
                       '13-BD.pdb', '2L63.pdb'};

N = [500; 5000; 50000; 500000; 5000000];
repeat = 10;
for i=1:size(structureSet,1)
    result(i).name = strcat(structureSet{i,1}, ' - ', structureSet{i,2});
    
    stats = [];  % Clear before starting next pair
    for j=1:length(N)
        for k=1:repeat
            [v(k), d] = dbkExvol(structureSet{i,1}, structureSet{i,2}, N(j), 'single');
        end
        stats = [stats; N(j), mean(v), std(v)];
    end
    result(i).stats = stats;
    % Get the box volume before continuing to the next set
    result(i).boxV = d.boxV;
end

clear i j k N repeat v d avg stdev