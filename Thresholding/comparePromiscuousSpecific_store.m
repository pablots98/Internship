function [enzymeData] = comparePromiscuousSpecific_store(spec,prom,modelData)

% USAGE:
% % [enzymeData] = comparePromiscuousSpecific(spec,prom,modelData)
% % calculates enzyme expression data

% INPUTS:
% % spec:       a specialist enzyme structure (see getSpecialistEnzymes)
% % prom:       a promiscuous enzyme structure (see getPromEnzymes)
% % modelData:  data for all the genes in the model

% OUTPUTS:
% % enzymeData:enzyme expression data for all the genes in the model

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

% EDITED BY LD:
    % 1. the function stores the info on enzymes and reactions for whose genes
    % we have no data in the dataset in the output structure
    % 2. the function considers enzymes for which we have enzymes made up of both genes
    % for which we have no data and genes for which we have data just
    % assigning the value of the gene for which we actually have data

% comparing distributions of specialist & promiscuous enzymes
specSubunits = regexp(spec.enzymes,' & ','split');
promSubunits = regexp(prom.enzymes,' & ','split');

missingspec=[];
missingspecrxns=[];
specExprMatrix = zeros(size(specSubunits,1),size(modelData.value,2));
for j=1:size(specSubunits,1)
        if sum(ismember(modelData.ID_geneMissing,specSubunits{j}))==length(specSubunits{j})
           missingspec=[missingspec; spec.enzymes(j)];
           missingspecrxns=[missingspecrxns; spec.rxns(j)];
        else
           for i=1:size(modelData.value,2)
               if sum(ismember(modelData.ID_geneMissing,specSubunits{j}))+sum(ismember(modelData.gene,specSubunits{j}))==length(specSubunits{j})
                specExprMatrix(j,i) = min(modelData.value(ismember(modelData.gene,specSubunits{j}),i));
               end
           end
        end
end
missingprom=[];
missingpromrxns=[];
promExprMatrix = zeros(size(promSubunits,1),size(modelData.value,2));
for j=1:size(promSubunits)
    if sum(ismember(modelData.ID_geneMissing,promSubunits{j}))==length(promSubunits{j})
       missingprom=[missingprom; prom.enzymes(j)];
       missingpromrxns=[missingpromrxns; prom.rxns(j)];
    else
        for i=1:size(modelData.value,2)
            if sum(ismember(modelData.ID_geneMissing,promSubunits{j}))+sum(ismember(modelData.gene,promSubunits{j}))==length(promSubunits{j})
            promExprMatrix(j,i) = min(modelData.value(ismember(modelData.gene,promSubunits{j}),i));
            end
        end
    end
end
enzAbund = [specExprMatrix; promExprMatrix];
specExprVec = reshape(specExprMatrix,numel(specExprMatrix),1);
promExprVec = reshape(promExprMatrix,numel(promExprMatrix),1);

specExprVec(specExprVec==-1) = [];
promExprVec(promExprVec==-1) = [];

figure;
hold on
histogram(log10(specExprVec),'BinEdges',[-8:0.25:5],'DisplayStyle','stairs');
histogram(log10(promExprVec),'BinEdges',[-8:0.25:5],'DisplayStyle','stairs');
hold off
legend('Specialist Enzymes','Promiscuous Enzymes','Location','NorthWest');

enzymeData.enzyme = [spec.enzymes;prom.enzymes];
enzymeData.value = enzAbund;
enzymeData.rxns = [spec.rxns;prom.rxns];
enzymeData.Tissue = modelData.Tissue;
enzymeData.missingEnz = [missingspec; missingprom];
enzymeData.missingRxns = [missingspecrxns; missingpromrxns];