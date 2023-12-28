Main scripts with the 4 algorithms we wanted to evaluate:
- **Global_Threshold_reactions.m**: Script that use just a percentile of the expression in a sample to establish the threshold, above which, the genes or reactions are considered core.
- **LocalT2_Threshold_reactions.m**: Script that use a function called localT2_new to establish different thresholds for each gene or reactions.
- **Localgini_threshold.m**: Script based on ... paper to create thresholds for each genes based on the gini coefficient.
- **Standep...**: Script that use the SanDep algorithm to set the thresholds and establish the core reactions.
