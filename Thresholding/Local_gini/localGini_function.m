function [adjusted_matrix, expression_scoreMatrix] = localGini_function(mappedDat, unpruneddata, lowerThres, upperThres)
    % Esta función ajusta los datos de expresión génica en función de los umbrales de percentil proporcionados.
    % Convierte los datos de entrada en un formato binario donde los valores iguales o mayores a 1 se establecen en 1,
    % y los valores menores a 1 se establecen en 0.

    if exist('lowerThres', 'var') && exist('upperThres', 'var') && lowerThres < 1 && upperThres < 1
        warning('The thresholds should be percentiles in the 0 to 100 range, you may be using the 0 to 1 range here.')
    end
    if ~exist('lowerThres', 'var')
        lowerThres = 25;
    end
    if ~exist('upperThres', 'var')
        upperThres = 75;
    end
    
    Contexts = mappedDat.Properties.VariableNames;
    mappedDat = table2array(mappedDat);
    unpruneddata = table2array(unpruneddata);
    unpruneddata(unpruneddata==0) = [];
    

    linData = reshape(unpruneddata(:, 1:end-1), [], 1);
    lowerThres = quantile(linData, lowerThres / 100);
    upperThres = quantile(linData, upperThres / 100);

    coreMat = false(size(mappedDat));
    expression_ths_local_Gini_25_75 = zeros(size(mappedDat, 1), 1);

    for i = 1:size(mappedDat, 1)
        expressionValue = mappedDat(i, :);

        % Calcular el coeficiente de Gini usando la fórmula del paper
        giniCoefficient = ginicoeff(expressionValue);

        % Calcular el umbral Localgini
        %localGiniThreshold = prctile(expressionValue, (giniCoefficient));

        % % Asignar el umbral adecuado basado en lowerThres y upperThres
        % if giniCoefficient >= upperThres
        %     expression_ths_local_Gini_25_75(i) = upperThres;
        % elseif giniCoefficient <= lowerThres
        %     expression_ths_local_Gini_25_75(i) = lowerThres;
        % else
        %     expression_ths_local_Gini_25_75(i) = giniCoefficient;
        end
        giniCoefficient(giniCoefficient>=upperThres)=upperThres;giniCoefficient(giniCoefficient<=lowerThres)=lowerThres;

    end
    
    gene_exp = expressionValue-repmat(giniCoefficient,1,numel(Contexts));
    gene_exp = gene_exp./std(gene_exp,[],2); % modified gene expression values 

    % expression_scoreLocal_Gini_25_75 = zeros(size(mappedDat));
    % for i = 1:size(mappedDat, 2)
    %    expression_scoreLocal_Gini_25_75(:, i) = mappedDat(:, i) ./ expression_ths_local_Gini_25_75;
    % end

    expression_scoreMatrix = gene_exp; % expression_scoreLocal_Gini_25_75
    adjusted_matrix = gene_exp >= 1; % expression_scoreLocal_Gini_25_75 >= 1
end

function [gcp]=ginicoeff(data)
    data_sort=sort(data,2);
    n=size(data,2);
    
    for i=1:size(data,1)
        x=data_sort(i,:);
        G_num=sum(((2*[1:n])-n-1).*x);
        G_den=sum(x)*n;
        gc(i)=(G_num/G_den)*100;
    end
    
    for i=1:numel(gc)
        gcp(i)=prctile(data(i,:),gc(i));
    end
    gcp=gcp';
end 

function [enz_exp,enzymes,newparsedGPR] = geneToEnzymeExp(geneExpression, geneNames, parsedGPR)
    enzymes=cell(0,0); % to store enzyme names
    enz_exp=[]; % to store the expression values of enzymes
    newparsedGPR=cell(size(parsedGPR));
    for i=1:numel(parsedGPR)
        current_rxn = parsedGPR{i};
        if ~isempty(current_rxn{1})
            newparsedGPR{i}=cell(size(current_rxn));
            for j=1:numel(current_rxn)
                currentArr = current_rxn{j};
                if length(currentArr)==1
                    newparsedGPR{i}{j}=currentArr;
                    if sum(ismember(enzymes,currentArr{1}))<1
                        ids = find(ismember(geneNames,currentArr{1}));
                        if ~isempty(ids)
                            currentExp = geneExpression(ids);
                            enzymes=[enzymes;currentArr{1}];
                            enz_exp=[enz_exp;currentExp];
                        end
                    end
                else
                    enz_name = strjoin(sort(currentArr),' and ');
                    newparsedGPR{i}{j} = cellstr(enz_name);
                    if sum(ismember(enzymes,enz_name))<1
                        ids=find(ismember(geneNames,currentArr));
                        if numel(ids)==numel(currentArr)
                            enzymes=[enzymes;enz_name];
                            currentExp = min(geneExpression(ids));
                            enz_exp=[enz_exp;currentExp];
                        end
                    end
                end
            end  
        else
            newparsedGPR{i}=current_rxn;
        end
        
    end   
end

function rxn_exp = enzymeToRxnExp(enz_exp,enzymes,parsedGPR)
    rxn_exp = -1*ones(length(parsedGPR),1);
    for i = 1:length(parsedGPR)
        current_rxn=parsedGPR{i};
        curArr=[];
        for j=1:length(current_rxn)
            if length(current_rxn{j})>=1
                ids = find(ismember(enzymes,current_rxn{j}));
                if ~isempty(ids) 
                     curArr= [curArr, enz_exp(ids)]; 
                end
            end
        end
        if ~isempty(curArr)
            rxn_exp(i)=max(curArr);
        end
    end
end
