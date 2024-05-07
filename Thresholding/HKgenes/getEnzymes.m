function [enzymesInfo] = getEnzymesFromGeneList(model, geneList)

% USAGE:
%   [enzymesInfo] = getEnzymesFromGeneList(model, geneList)
%
% INPUTS:
%   model:      a COBRA model for which a list of enzymes are desired
%   geneList:   a cell array of strings containing the gene identifiers to be considered
%
% OUTPUTS:
%   enzymesInfo:  a structure containing information about enzymes associated with the specified genes
%
% AUTHORS:
%   Chintan Joshi:  for StanDep paper (May 2018)

% Parse and arrange GPRs
parsedGPR = GPRparser(model);
[parsedGPR, ix] = linearization_index(parsedGPR, 'rows');
corrRxns = model.rxns(ix);
corrSys = model.subSystems(ix);

% Filter out empty parsed GPRs and unrelated genes
ix = find(cellfun(@isempty, parsedGPR) | ~cellfun(@(x) any(ismember(geneList, x)), parsedGPR));
parsedGPR(ix) = []; 
corrRxns(ix) = [];
corrSys(ix) = [];

% Join GPR strings
for i = 1:length(parsedGPR)
    if length(parsedGPR) ~= 1
        parsedGPR{i,1} = strjoin(parsedGPR{i}, ' & ');
    end
end

% Find unique parsed GPRs and their associated reactions and systems
ugprs = unique(parsedGPR);
rxns = cell(length(ugprs), 1);
subSystems = cell(length(ugprs), 1);

for i = 1:length(ugprs)
    rxns{i, 1} = corrRxns(ismember(parsedGPR, ugprs{i}));
    subSystems{i, 1} = corrSys(ismember(parsedGPR, ugprs{i}));
end

% Store results
enzymesInfo.enzymes = ugprs;
enzymesInfo.rxns = rxns;
enzymesInfo.subSystems = subSystems;
nrxns = cellfun(@length, rxns);

% Keep only enzymes associated with at least one reaction
validEnzymes = nrxns >= 1;
enzymesInfo.enzymes = enzymesInfo.enzymes(validEnzymes);
enzymesInfo.rxns = enzymesInfo.rxns(validEnzymes);
enzymesInfo.subSystems = enzymesInfo.subSystems(validEnzymes);
enzymesInfo.nrxns = nrxns(validEnzymes);

end
