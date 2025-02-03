% Get list of filenames in ./beta_calc/bases and ./beta_calc/solutes

%clear all

% bases = dir('./beta_calc/bases/');
% solutes = dir('./beta_calc/solutes/');
N = 1000000;

% numBases = numel(bases) - 2;  % -2 b/c one is '.' and another is '..'
% numSolutes = numel(solutes) - 2;
% for i=1:numBases
%     baseName{i,1} = bases(i+2,1).name;
% end
% for i=1:numSolutes
%     soluteName{i,1} = solutes(i+2,1).name;
% end

baseName = {'ade.pdb'; 'thy.pdb'; 'cyt.pdb'; 'hyp.pdb'; 'tbr.pdb'; ...
            'tpy.pdb'; '5au.pdb'};
soluteName = {'eg.pdb'; 'dieg_linear.pdb'; 'trieg_linear.pdb'; ...
              'tetraeg_linear.pdb'; 'glycerol.pdb'; 'etoh.pdb'; 'meoh.pdb'; ...
              '1proh.pdb'; '2proh.pdb'; '13bd.pdb'; '14bd.pdb'; 'dioxane.pdb'};

numBases = numel(baseName); numSolutes = numel(soluteName);

k=1;
for i=1:numSolutes
    for j=1:numBases
        fprintf('%s:  %s - %s \n', 'Now calculating', baseName{j,1}, soluteName{i,1})
        pair{k,1} = [baseName{j,1} ', ' soluteName{i,1}];
        result{k,1} = dbkBeta(['./beta_calc/bases/', baseName{j,1}], ...
                             ['./beta_calc/solutes/', soluteName{i,1}], ...
                              N, 0.1, 'single');
        k = k + 1;
    end
end

% Produce a list of all the unique atom types on the bases
% Atom names for bases are the ROW names of the beta matrix
% Atom names for the solutes are the COL names of the beta matrix
% this is due to the order in which the filenames are pssed to dbkBeta(),
% above
allBeta1RowNames = [];  allBeta1ColNames = [];
for i=1:size(result,1)
    currResult = result{i,1};
    allBeta1RowNames = [allBeta1RowNames; currResult.beta1RowNames];
    allBeta1ColNames = [allBeta1ColNames, currResult.beta1ColNames'];
end
allBeta1RowNames = unique(allBeta1RowNames); allBeta1ColNames = unique(allBeta1ColNames);
allBeta1ColNames = {'OH', 'OG', 'CA'};  % Hardcoded HACK to make them come out in the right order
clear currResult

% Reorder the rows and cols in all of the beta1 matricies to correspond to
% the order in allBeta1RowNames and allBeta1ColNames
allBeta1Sorted = {};
for i=1:size(result,1)
    currResult = result{i,1};
    currBeta1 = currResult.beta1;
    currBeta1RowNames = currResult.beta1RowNames;
    currBeta1ColNames = currResult.beta1ColNames';  % Note transpose
    
    currBeta1SortedRows = zeros(size(allBeta1RowNames,1), size(currBeta1,2));
    currBeta1SortedRowsCols = zeros(size(allBeta1RowNames,1), size(allBeta1ColNames,2));
    
    % First sort the rows
    %currBeta1
    %allBeta1RowNames
    %currBeta1RowNames
    for j=1:size(allBeta1RowNames,1)
        rowIDX = find(strcmp(allBeta1RowNames{j,1}, currBeta1RowNames));  % Should be only 1 result
        
        if isempty(rowIDX)
            currBeta1SortedRows(j,:) = zeros(1,size(currBeta1,2));
        else
            currBeta1SortedRows(j,:) = currBeta1(rowIDX,:);
        end
        % currBeta1SortedRows ~ nRowsTot*nColsCurrBeta1
    end
    % Sort the cols
    %currBeta1SortedRowsCols
    for j=1:size(allBeta1ColNames,2)
        colIDX = find(strcmp(allBeta1ColNames(1,j), currBeta1ColNames));  % Should be only 1 result
        
        if isempty(colIDX)
            currBeta1SortedRowsCols(:,j) = zeros(size(allBeta1RowNames,1),1);
        else
            currBeta1SortedRowsCols(:,j) = currBeta1SortedRows(:,colIDX);
        end
        % currBeta1SortedRowsCols ~ nRowsTot*nColsTot
    end
    %currBeta1SortedRowsCols
    allBeta1Sorted{i,1} = currBeta1SortedRowsCols;
    
    clear rowIDX colIDX j currResult currBeta1 currBeta1RowNames currBeta1ColNames currBeta1SortedRows
    clear currBeta1SortedRowsCols
end

% Hardcoded hack:  Add rows 2 and 4 of all the sorted beta1 matricies to
% combine CR and NR into "Ring" ... sort into proper order (hack)
currBeta1Sorted_Ring = [];
allBeta1Sorted_Ring = {};
for i=1:size(allBeta1Sorted,1)
    currBeta1Sorted = allBeta1Sorted{i,1};
    currBeta1Sorted_Ring = [];
    
    currBeta1Sorted_Ring(1,:) = currBeta1Sorted(5,:); % 'O'
    currBeta1Sorted_Ring(2,:) = currBeta1Sorted(1,:); % 'CM'
    currBeta1Sorted_Ring(3,:) = currBeta1Sorted(3,:); % 'N'
    currBeta1Sorted_Ring(4,:) = currBeta1Sorted(2,:) + currBeta1Sorted(4,:);
    
    %currBeta1Sorted_Ring
    allBeta1Sorted_Ring{i,1} = currBeta1Sorted_Ring;
end
clear currBeta1Sorted currBeta1Sorted_Ring

% Flatten all the beta1 matricies into row vectors
beta1Prod = [];
for i=1:size(allBeta1Sorted_Ring,1)
    currAllBeta1Sorted_Ring = allBeta1Sorted_Ring{i,1};
    beta1Prod(i,:) = reshape(currAllBeta1Sorted_Ring, 1, 12);
end
clear currAllBeta1Sorted_Ring

% Get Beta1Prod rownames
for i=1:size(pair,1)
    [a, b] = strtok(pair{i,1}, '.');
    [c, d] = strtok(b, '.pdb, ');
    Beta1ProdRowNames{i,1} = strcat(a, '*', c);
end
clear a b c d

% k=1;
% for i=1:numBases
%     for j=1:numSolutes
%             Beta1ProdRowNames{k,1} = strcat(strtok(soluteName{j,1}, '.'), '*', strtok(baseName{i,1}, '.'));
%         k = k + 1;
%     end
% end

clear N bases solutes numSolutes numBases k i j baseName soluteName ans