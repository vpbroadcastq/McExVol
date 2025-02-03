% Get list of filenames in ./beta_calc/bases and ./beta_calc/solutes

clear all

bases = dir('./beta_calc/bases/');
solutes = dir('./beta_calc/solutes/');

numBases = numel(bases) - 2;  % -2 b/c one is '.' and another is '..'
numSolutes = numel(solutes) - 2;
for i=1:numBases
    baseName{i,1} = bases(i+2,1).name;
end
for i=1:numSolutes
    soluteName{i,1} = solutes(i+2,1).name;
end

for i=1:numBases
    o = dbkPDBToObject(['./beta_calc/bases/', baseName{i,1}]);
    baseAName{i,1} = o.atomName;
end

for i=1:numSolutes
    o = dbkPDBToObject(['./beta_calc/solutes/', soluteName{i,1}]);
    soluteAName{i,1} = o.atomName;
end

clear o i bases solutes numBases numSolutes