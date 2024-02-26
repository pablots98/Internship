%%%%%%%%%%%%%%%          RESULTS                        %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Load data                                            %%
sampleNames = readtable("sampleNames.xlsx", 'ReadVariableNames', false);
HK_G_acc_G = readtable("HK_G_acc_G.xlsx");
HK_R_acc_G = readtable("HK_R_acc_G.xlsx");
HK_G_acc_LT = readtable("HK_G_acc_LT.xlsx");
HK_R_acc_LT = readtable("HK_R_acc_LT.xlsx");
HK_G_acc_LG = readtable("HK_G_acc_LG.xlsx");
HK_R_acc_LG = readtable("HK_R_acc_LG.xlsx");
HK_G_acc_SD = readtable("MatchComparisonResults_G.xlsx")
HK_R_acc_SD = readtable("MatchComparisonResults.xlsx");

%%                     Create results table                              %%
numRows = height(sampleNames);
results_table = table('Size', [numRows 8], 'VariableTypes', repmat({'double'}, 1, 8));
RowNames = table2cell(sampleNames);
results_table.Properties.RowNames = RowNames;

%%                     Add the results                                   %%

results_table.(1) = HK_G_acc_G.(1);
results_table.Properties.VariableNames{1} = HK_G_acc_G.Properties.VariableNames{1};

results_table.(2) = HK_G_acc_LT.(1);
results_table.Properties.VariableNames{2} = HK_G_acc_LT.Properties.VariableNames{1};

results_table.(3) = HK_G_acc_LG.(1);
results_table.Properties.VariableNames{3} = HK_G_acc_LG.Properties.VariableNames{1};

results_table.(5) = HK_G_acc_SD.(1);
results_table.Properties.VariableNames{4} = HK_G_acc_SD.Properties.VariableNames{1};

results_table.(5) = HK_R_acc_G.(1);
results_table.Properties.VariableNames{5} = HK_R_acc_G.Properties.VariableNames{1};

results_table.(6) = HK_R_acc_LT.(1);
results_table.Properties.VariableNames{6} = HK_R_acc_LT.Properties.VariableNames{1};

results_table.(7) = HK_R_acc_LG.(1);
results_table.Properties.VariableNames{7} = HK_R_acc_LG.Properties.VariableNames{1};

results_table.(8) = HK_R_acc_SD.(1);
results_table.Properties.VariableNames{8} = HK_R_acc_SD.Properties.VariableNames{1};

%%                      transpose the matrix                             %%
transposed_table = array2table(table2array(results_table)');
transposed_table.Properties.RowNames = results_table.Properties.VariableNames;
transposed_table.Properties.VariableNames = results_table.Properties.RowNames;

%%                          Save the table                               %%
writetable(transposed_table, 'transposed_results_table.xlsx', 'WriteRowNames',true, 'WriteVariableNames',true);

%% TEST
%clear all;
%%                          load the data                                %%
% Sample Names
sampleNames = readtable("sampleNames.xlsx", 'ReadVariableNames', false);

% Global Thresholding results
HK_G_acc_G = readtable("HK_G_acc_G.xlsx");
HK_R_acc_G = readtable("HK_R_acc_G.xlsx");

% localT2 results
HK_G_acc_LT = readtable("HK_G_acc_LT.xlsx");
HK_R_acc_LT = readtable("HK_R_acc_LT.xlsx");

% LocalGini results
HK_G_acc_LG = readtable("HK_G_acc_LG.xlsx");
HK_R_acc_LG = readtable("HK_R_acc_LG.xlsx");

% StanDep results
HK_G_acc_SD = readtable("MatchComparisonResults_G.xlsx")
HK_R_acc_SD = readtable("MatchComparisonResults.xlsx");

%%                         Create the table                              %%
numRows = height(sampleNames);
results_table = table('Size', [numRows 8], 'VariableTypes', repmat({'double'}, 1, 8));

% Set the row names as the sample names
RowNames = table2cell(sampleNames);
results_table.Properties.RowNames = RowNames;

%%                          Add the results                              %%
% HK genes for global thresholding
columnName = HK_G_acc_G.Properties.VariableNames{1};
results_table.("Var1") = HK_G_acc_G.(columnName);
results_table.Properties.VariableNames{1} = columnName;

% HK genes for LocalT2 thresholding
columnName = HK_G_acc_LT.Properties.VariableNames{1};
results_table.("Var2") = HK_G_acc_LT.(columnName);
results_table.Properties.VariableNames{2} = columnName;

% HK genes for LocalGini thresholding
columnName = HK_G_acc_LG.Properties.VariableNames{1};
results_table.("Var3") = HK_G_acc_LG.(columnName);
results_table.Properties.VariableNames{3} = columnName;

% HK genes for StanDep thresholding
columnName = HK_G_acc_SD.Properties.VariableNames{1};
results_table.("Var4") = HK_G_acc_SD.(columnName);
results_table.Properties.VariableNames{4} = columnName;

% HK Reactions for Global Thresholding
columnName = HK_R_acc_G.Properties.VariableNames{1};
results_table.("Var5") = HK_R_acc_G.(columnName);
results_table.Properties.VariableNames{5} = columnName;

% HK Reactions for LocalT2 Thresholding
columnName = HK_R_acc_LT.Properties.VariableNames{1};
results_table.("Var6") = HK_R_acc_LT.(columnName);
results_table.Properties.VariableNames{6} = columnName;

% HK Reactions for LocalGini Thresholding
columnName = HK_R_acc_LG.Properties.VariableNames{1};
results_table.("Var7") = HK_R_acc_LG.(columnName);
results_table.Properties.VariableNames{7} = columnName;

% HK Reactions for StanDep Thresholding
columnName = HK_R_acc_SD.Properties.VariableNames{1};
results_table.("Var8") = HK_R_acc_SD.(columnName);
results_table.Properties.VariableNames{8} = columnName;

%%                       SAVE THE RESULTS TABLE                          %%
writetable(results_table, 'results_table.xlsx', 'WriteRowNames',true, 'WriteVariableNames',true);                       