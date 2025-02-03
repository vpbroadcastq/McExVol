% Get list of filenames in ./beta_calc/bases and ./beta_calc/solutes

clear all

bases = dir('./beta_calc/bases/');
solutes = dir('./beta_calc/solutes/');
N = 1000000;

numBases = numel(bases) - 2;  % -2 b/c one is '.' and another is '..'
numSolutes = numel(solutes) - 2;
for i=1:numBases
    baseName{i,1} = bases(i+2,1).name;
end
for i=1:numSolutes
    soluteName{i,1} = solutes(i+2,1).name;
end

% % Pairwise excluded volume
% k=1;
% for i=1:numBases
%     for j=1:numSolutes
%         fprintf('%s:  %s - %s \n', 'Now calculating', baseName{i,1}, soluteName{j,1})
%         pair{k,1} = [baseName{i,1} ', ' soluteName{j,1}];
%         [v, v23Ex_d{k,1}] = dbkExvol(['./beta_calc/bases/', baseName{i,1}], ...
%                                        ['./beta_calc/solutes/', soluteName{j,1}], ...
%                                         N, 'single');
%         v23Ex_volumes(k,1) = v23Ex_d{k,1}.Vex;
%         k = k + 1;
%     end
% end

% Absolute volumes
for i=1:numSolutes
    fprintf('%s:  %s %s \n', 'Now calculating', soluteName{i,1}, 'absolute volume')
    [vSolute(i,1), vSolutes_d{i,1}] = dbkExvol(['./beta_calc/solutes/', soluteName{i,1}], ...
                                       'sphere3.pdb', ...
                                       N, 'single');
end
for i=1:numBases
    fprintf('%s:  %s %s \n', 'Now calculating', baseName{i,1}, 'absolute volume')
    [vBases(i,1), vBases_d{i,1}] = dbkExvol(['./beta_calc/bases/', baseName{i,1}], ...
                                'sphere3.pdb', ...
                                N, 'single');
end



clear bases solutes N numBases numSolutes baseName soluteName i j k v
