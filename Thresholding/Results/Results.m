%%%%%%%%%%%%%%%          RESULTS                        %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                          load the data                                %%
% Global Thresholding results
HK_G_acc_G = readtable("HK_G_acc_G.xlsx");
HK_R_acc_G = readtable("HK_R_acc_G.xlsx");

% localT2 results
HK_G_acc_LT = readtable("HK_G_acc_LT.xlsx");
HK_R_acc_LT = readtable("HK_R_acc_LT.xlsx");

% LocalGini results
HK_G_acc_LG = readtable("HK_G_acc_LG.xlsx");
HK_R_acc_LG = readtable("HK_R_acc_LG.xlsx");

%%                         Create the table                              %%

results_table = table;


