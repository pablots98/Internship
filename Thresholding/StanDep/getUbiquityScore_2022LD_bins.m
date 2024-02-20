function [ubiScore,uScore] = getUbiquityScore_2022LD_bins(clustObj,model)

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

%this updated (by LD) version of the function
   %1. assigns to missing-data rxns & to the gpr-less rxns the ubiValue of
   %the median of the intermediate ubiScore in a sample specific way (exluding 1
   %and -1e-6)
   %2. changes the value of the no-data reactions from -1 to -1e-6 as
   %requested by rankReactions in mCADRE extraction
   %3. puts back to 1e-6 the uScore for partially expressed enzymes 

 if isfield(clustObj,'missingobjectMaps')
    missingrxns_idx=[];
    kk=0;
    for f=1:size(clustObj.missingobjectMaps,1)
        if ischar(clustObj.missingobjectMaps{f})
            kk=kk+1;
            missingrxns_idx(kk,1)=find(ismember(model.rxns,clustObj.missingobjectMaps{f}));
        else
            a=size(clustObj.missingobjectMaps{f},1);
            for ff=1:a
                kk=kk+1;
                missingrxns_idx(kk,1)=find(ismember(model.rxns,clustObj.missingobjectMaps{f,1}{ff,1}));
            end
        end
    end
 end

objDist = zeros(size(clustObj.Data));
uci = 1:1:size(clustObj.C,1);
cidx = clustObj.cindex;
edgeX=clustObj.edgeX;

if isfield(clustObj,'partialzeroIdx')
    zeroidx=clustObj.partialzeroIdx;
end

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

if isfield(clustObj,'partialzeroIdx')
    for g=1:size(zeroidx,1)
        uScore(zeroidx(g,1), zeroidx(g,2)) = -1e-6;
    end
end

[uScore2,rxns] = openUbiquityMatrix(clustObj,uScore);
ubiScore = repmat(-1e-6,length(model.rxns),size(clustObj.Data,2));
for i=1:length(model.rxns)
    if sum(ismember(rxns,model.rxns{i}))~=0
        ubiScore(i,:) = max(uScore2(ismember(rxns,model.rxns{i}),:),[],1);            
    end
end
ubiScore(ubiScore==-inf) = -1e-6;
for i = 1:size(ubiScore, 2)
    ubi = ubiScore(:, i); ubi(ubi == -1e-6 | ubi == 1) = [];
    for j=1:size(ubiScore,1)
      if isfield(clustObj,'missingobjectMaps')
        if ismember(j,missingrxns_idx) & all(ubiScore(j,:) == -1e-6, 2) 
        ubiScore(j,i)=median(ubi);
        end
      end
      if isempty(model.grRules{j})
        ubiScore(j,i)=median(ubi);
      end
    end
end

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