function plotGeneFrequency(geneFreq)
    % Extract gene names and their frequencies
    geneNames = keys(geneFreq);
    frequencies = values(geneFreq);
    frequencies = cell2mat(frequencies);  % Convert from cell array to numeric array

    % Create a bar chart
    figure;
    bar(frequencies, 'FaceColor', [0.2, 0.2, 0.5]);
    set(gca, 'XTickLabel', geneNames, 'XTick', 1:numel(geneNames));
    xlabel('Gene ID');
    ylabel('Frequency');
    title('Frequency of Housekeeping Genes in Gene Associations');
    grid on;
end
