function [ubiScore,uScore] = getUbiquityScore_2022(clustObj,edgeX,model)

% USAGE:
% % [ubiScore,uScore] = getUbiquityScore(clustObj,edgeX,model)
% % code needed to calculate inputs for mCADRE

% INPUTS:
% % clustObj:   cluster object calculated in geneExprDist_hierarchy
% % edgeX:      bins used in clustObj
% % model:      a COBRA model to be used

% OUTPUTS:
% % ubiScore:   a matrix describing ubiquity score of each reaction
% % uScore:     a matrix describing ubiquity score of enzymes

% AUTHORS:
% % Chintan Joshi:  for StanDep paper (May 2018)

objDist = zeros(size(clustObj.Data));
uci = 1:1:size(clustObj.C,1);
cidx = clustObj.cindex;
[~,thrVal] = clusterVariability1(clustObj,edgeX,false,0,[1 1]);
for j=1:size(clustObj.Data,2)
    for i=1:length(uci)
        ic = find(cidx==uci(i));
        objDist(ic,j) = clustObj.Data(ic,j) - (thrVal(i));
    end
end
uScore = objDist; uScore(uScore > 0) = 1;

for i=1:size(uScore,2)
    indx = find(uScore(:,i)~=1);
    m = uScore(:,i); m(m==-inf) = []; m = min(m);
    uScore(indx,i) = 1 - uScore(indx,i)/m;
end

[uScore2,rxns] = openUbiquityMatrix(clustObj,uScore);
ubiScore = repmat(-1,length(model.rxns),size(clustObj.Data,2));
for i=1:length(model.rxns)
    if sum(ismember(rxns,model.rxns{i}))~=0
        ubiScore(i,:) = max(uScore2(ismember(rxns,model.rxns{i}),:),[],1);
    end
end
ubiScore(ubiScore==-inf) = -1;

function [openUbi,rxns] = openUbiquityMatrix(clustObj,uMat)

    k = 0;
    for ii=1:size(uMat,1)
        if ischar(clustObj.objectMaps{ii})
            k = k+1;
            openUbi(k,:) = uMat(ii,:);
            rxns{k,1} = clustObj.objectMaps{ii};
        else
            n = length(clustObj.objectMaps{ii});
            for jj=1:n
                k = k+1;
                openUbi(k,:) = uMat(ii,:);
                rxns{k,1} = clustObj.objectMaps{ii}{jj,1};
            end
        end
    end
end

end