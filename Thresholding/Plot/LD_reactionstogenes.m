%run the code StanDep_mCadre
% cluster.objects and uScore are the two outputs you will need

%sample OfInterest, insert the coloumn for which you want to make the
%comparison
uScoreOI=uScore(:,4);
ixrow=find(uScoreOI==1);

%enzymes that have a uScore of 1 are
enzymesOI=clustObj.objects(ixrow);

%open the enzyme matrix
for i = 1:length(enzymesOI)
    if contains(enzymesOI{i}, ' & ')
        splitParts = strsplit(enzymesOI{i}, ' & ');
        for ii=1:length(splitParts)
            genes{i,ii}=splitParts{ii};
        end
    else
        genes{i,1}=enzymesOI{i};
    end
end

linearGenes = reshape(genes.', [], 1);
linearGenes = linearGenes(~cellfun('isempty', linearGenes));
linearGenes = unique(linearGenes);

%what I get are the genes that are part of enzymes that have1 as uScore