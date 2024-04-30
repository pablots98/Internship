function [geneFreq, geneFreqTable, geneFreqCell] = countGeneFrequency(hkDetails, hkGenes)
    % Initialize the frequency map
    geneFreq = containers.Map(hkGenes, num2cell(zeros(length(hkGenes), 1)));
    
    % Iterate over each entry in hkDetails
    for i = 1:length(hkDetails)
        for j = 1:length(hkGenes)
            % Count occurrences of each gene in the geneAssociations
            geneCount = count(hkDetails(i).geneAssociations, hkGenes{j});
            geneFreq(hkGenes{j}) = geneFreq(hkGenes{j}) + geneCount;
        end
    end

    % Extract gene IDs and their frequencies
    geneNames = keys(geneFreq);
    frequencies = values(geneFreq);
    frequencies = cell2mat(frequencies);  % Convert from cell array to numeric array

    % Create a cell array with gene IDs and their frequencies
    % Ensure both are column vectors
    geneFreqCell = [geneNames(:), num2cell(frequencies(:))];

    % Create a table from the cell array for better readability
    geneFreqTable = cell2table(geneFreqCell, 'VariableNames', {'GeneID', 'Frequency'});
end
