function [results, ListResults] = findRxnsFromGenes_Alicia(model, genes, ListResultsFlag)
% Print every reaction associated with a gene of interest
%
% USAGE:
%
%    [results, ListResults] = findRxnsFromGenes(model, genes, ListResultsFlag)
%
% INPUTS:
%    model:              COBRA model structure
%    genes:              string of single gene or cell array of multiple
%                        genes for which `rxns` are desired.
%
% OPTIONAL INPUTS:
%    ListResultsFlag:    1 if you want to output `ListResults` (Default = 0)
%
% OUTPUTS:
%    results:            structure containing cell arrays for each gene.
%                        Each cell array has one column of rxn abbreviations
%                        and one column containing the reaction formula
%    ListResults:        same as above, but in a cell array
%
% .. Author:
%       - Nathan Lewis 02/16/08
%       - edited 04/05/09 (MUCH faster now -- NL)
%       - edited 06/11/10 (yet even faster now -- NL)
%       - Ronan, Ines - edited interface for backward compatibility

if nargin< 3
    ListResultsFlag = 0;
end

disp("Input: string of single gene or cell array of multiple genes")
if ~iscell(genes)
    genes = {genes}; %convert genes to a cell array 
else %check if any nested cells in array
    addGenes = {}; %fill if nested cells in array
    for i = 1:length(genes)
        if iscell(genes{i})
            if length(genes{i}) == 1
                genes(i) = genes{i};
            else %more than one gene listed in nested cell
                addGenes = union(addGenes, genes{i});
                delGenes = i;
            end
        end
    end
    if ~isempty(addGenes)
        genes(delGenes) = []; %delete nested cell
        genes = union(genes, addGenes); %add genes from nested cell to list
    end
end

if isfield(model, 'geneNames') &&  isfield(model, 'genes')
    model.geneNames = regexprep(model.geneNames,'-','_DASH_');
    model.geneNames = regexprep(model.geneNames,'\.','_POINT_');
    model.genes = upper(regexprep(model.genes,'-','_DASH_')); 
    model.genes = upper(regexprep(model.genes,'\.','_POINT_'));
elseif isfield (model, 'genes')%to stay compatible with old style models where only 'genes' is present
    model.genes = regexprep(model.genes,'-','_DASH_');
    model.genes = regexprep(model.genes,'\.','_POINT_');
else 
    disp("The model needs to contain at least a 'genes' field")
end

genes = regexprep(genes,'-','_DASH_');
genes = regexprep(genes,'\.','_POINT_');

for i = 1:numel(genes)
    if ~startsWith(genes{i}, 'ENSG') % Check if the input is an Ensembl gene ID
        if isfield(model, 'geneNames') 
            % Find common elements between model.genes and genes{i}
            indModelGenes = find(ismember(model.geneNames, genes{i})); 
            if ~isempty(indModelGenes)
                % Index into model.geneNames using indices from model.genes
                genes{i} = model.genes{indModelGenes};
            end
        else 
            error('The input must be an Ensembl stable ID for this model')
        end
    end 
end

genes = unique(genes);

%find where the genes are located in the genes
GeneID = zeros(size(genes));
[~, geneIndModel, inModel] = intersect(model.genes, genes);
GeneID(inModel) = geneIndModel;%set location of genes in model

results = struct();
ListResults = {};
if any(GeneID == 0)
    warning('A gene was not found in the model!')
    if any(GeneID > 0)
        Ind = find(GeneID == 0);
        numGenesNotFound = numel(Ind); % Contar el número de genes no encontrados
        
        % Mostrar los IDs de Ensembl o símbolos de genes que no se encontraron en el modelo
        disp(['Number of genes not found in the model: ', num2str(numGenesNotFound)]); % Mostrar el número exacto de genes no encontrados
        %Display the Ensembl IDs or gene symbols of genes that were not found in the model for informational purposes (Alicia)
        disp(['The following genes were not found in the model: ', strjoin(genes(Ind), ', ')]);
        GeneID(Ind) = [];
        genes(Ind) = [];
    else
        return
    end
end

for i = 1:length(GeneID)
    %Ensures that geneids can become field names for structures
    tempGene = regexprep(genes{i}, '[^a-zA-Z0-9_]', '_');
    
    %If gene starts with a digit it cannot be a field name, prepend gene_ to correct
    %tempGene = cat(2, 'gene_', tempGene);
    
    %Reaction locations in model using rules field
    Ind_rxns = find(~cellfun(@isempty, strfind(model.rules, ...
        ['x(', num2str(GeneID(i)), ')'])));
    
    %Create gene field in results structure
    results.(tempGene) = cell(length(Ind_rxns), 5);
    
    %Fill in results
    results.(tempGene)(:, 1) = model.rxns(Ind_rxns);
    results.(tempGene)(:, 2) = printRxnFormula(model, model.rxns(Ind_rxns), 0);
    if isfield(model,'subSystems')
        results.(tempGene)(:, 3) = model.subSystems(Ind_rxns);
    end
    if isfield(model,'rxnNames')
        results.(tempGene)(:, 4) = model.rxnNames(Ind_rxns);
    end
    
    if isfield(model, 'geneNames')
        %Create an geneNames field in result structure (Alicia) 
        [~, geneIndModel, ~] = intersect(model.genes, tempGene);
        results.(tempGene)(:, 5) = model.geneNames(geneIndModel);
    end 

    if ListResultsFlag
        LR_RowCnt = size(ListResults, 1);
        ListResults(LR_RowCnt + 1 : LR_RowCnt + size(results.(tempGene), 1), 1:5) = results.(tempGene);
        ListResults(LR_RowCnt + 1 : LR_RowCnt + size(results.(tempGene), 1), 6) = {tempGene};
    end
end

%Convert the result cell array to a table (Alicia) 
%variableNames = {'Reaction', 'Equation', 'Subsystem', 'Reaction Name', 'Ensembl ID', 'Gene name'};
%ListResults = cell2table(ListResults, 'VariableNames', variableNames);

if isempty(results)
    warning('Your gene(s) was/were not associated with any reactions!')
end
