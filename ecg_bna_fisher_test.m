function p = ecg_bna_fisher_test(unit_counts_rest, unit_counts_task)
% ecg_bna_fisher_test
% This function computes the exact Fisher's test for counts of 
% heart-modulated (pos./neg.)units. Fisher's p is computed separately for 
% positively and negatively modulated units.
%
% USAGE:
%	p = ecg_bna_fisher_test(unit_counts_task, unit_counts_rest)
%
% INPUTS:
%       unit_counts_task - a 3-element vector with unit counts for
%       positively, negatively, and non-modulated units for task
%
%       unit_counts_rest - a 3-element vector with unit counts for
%       positively, negatively, and non-modulated units for rest
%
%       session_info - the info on the current session that contains data
%       paths (path to waveclus files is required here)
%       
% OUTPUTS:
%		p - a 2-element vector with Fisher's p-values for positively and 
%       negatively modulated units (without any correction for multiple 
%       comparisons)
%
% Author(s):	L.N. Vasileva, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2023-12-04:	Created function
% 
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

% % for debugging - unit counts
% count_task_pos = 20;
% count_task_neg = 21;
% count_task_non = 97;
% 
% count_rest_pos = 8;
% count_rest_neg = 13;
% count_rest_non = 117;

count_rest_pos = unit_counts_rest(1);
count_rest_neg = unit_counts_rest(2);
count_rest_non = unit_counts_rest(3);

count_task_pos = unit_counts_task(1);
count_task_neg = unit_counts_task(2);
count_task_non = unit_counts_task(3);

% Fisher's exact test for positively modulated neurons
contingency_pos = [count_rest_pos, count_task_pos; ...
                   sum([count_rest_neg, count_rest_non]), sum([count_task_neg, count_task_non])];

[~, pval_fisher_pos] = fishertest(contingency_pos);
disp(['Fisher''s exact test p-value for positively modulated: ', num2str(pval_fisher_pos)]);

% Fisher's exact test for negatively modulated neurons
contingency_neg = [count_task_neg, count_rest_neg; ...
                   sum([count_task_pos, count_task_non]), sum([count_rest_pos, count_rest_non])];

[~, pval_fisher_neg] = fishertest(contingency_neg);
disp(['Fisher''s exact test p-value for negatively modulated: ', num2str(pval_fisher_neg)]);

p = [pval_fisher_pos pval_fisher_neg];

end


