% Author: Hashem Sadeghiyeh
% March 16, 2020
% Paper: On the psychology of the psychology subject pool: an exploratory test of the good student effect

clear all;
[A, title, ~] = xlsread ('time_of_semester.xlsx');%105 item
l = size(title);
l = l(1,2);
title = strrep(title (1,6:l),'_',' ');
title = strrep(title,'-11','');
title(18:51) = lower(title(18:51));
title = regexprep(title,'(\<[a-z])','${upper($1)}');
title = strrep(title,'Epistemic Curiosity','');
title = strrep(title,'Dispositional Positive Emotions','DPE');
title = strrep(title,'Arnett','');
title = strrep(title,'Revised','');
title = strrep(title,'Santa Barbara','');
title = strrep(title,'Today','Discounting');
%% date
for i=5:l-1
    [R,P] = corrcoef (A(:,2), A(:,i), 'rows', 'complete');
    c(1,i-4).title = (title(i-4));
    c(1,i-4).r = R(1,2);
    c(1,i-4).p = P(1,2);
    t=0;
    for j=1:485
        if ~isnan(A(j,i))
            t=t+1;
        end
    end
    c(1,i-4).n = t;
end
%% time
for i=5:l-1
    [R,P] = corrcoef (A(:,3), A(:,i), 'rows', 'complete');
    c(2,i-4).title = (title(i-4));
    c(2,i-4).r = R(1,2);
    c(2,i-4).p = P(1,2);
    t=0;
    for j=1:485
        if ~isnan(A(j,i))
            t=t+1;
        end
    end
    c(2,i-4).n = t;
end 
%%
t=0;
for i=1:l-5
    if ~isempty (c(1,i).title) | ~isempty (c(2,i).title)
        if c(1,i).p < .05
            t= t+1;
            if ~isempty (c(1,i).title)
                C(t).trait = c(1,i).title;
            else
                C(t).trait = c(2,i).title;
            end
            C(t).R_date = c(1,i).r;
            C(t).P_date = c(1,i).p;
            C(t).R_time = c(2,i).r;
            C(t).P_time = c(2,i).p;
            C(t).N = c(2,i).n;
        end
    end
end

%% Multiple Comparisons Correction: Bonferroni-Holm Method

% Date
for i=5:l-1
    [R,P] = corrcoef (A(:,2), A(:,i), 'rows', 'complete');
         c(i-4).title = (title(i-4));
         c(i-4).r = R(1,2);
         c(i-4).p = P(1,2);
end
% sort
for i=1:l-5
    S(i,1) = i;
    S(i,2) = c(i).p;
end
S2 = sortrows(S, 2);
%
[corrected_p, h]=bonf_holm(S2(:,2),0.05);
%
for i=1:l-5
    if round(corrected_p(i,1),2) >= 0.05
        break
    end
    d = S2(i,1);
    Final_date(i).title = c(d).title;
    Final_date(i).r = c(d).r;
    Final_date(i).p = S2(i,2);
    Final_date(i).p_adjusted = corrected_p(i,1);
end
% time
for i=5:l-1
    [R,P] = corrcoef (A(:,3), A(:,i), 'rows', 'complete');
         c(i-4).title = (title(i-4));
         c(i-4).r = R(1,2);
         c(i-4).p = P(1,2);
end
% sort
for i=1:l-5
    S(i,1) = i;
    S(i,2) = c(i).p;
end
S2 = sortrows(S, 2)
%
[corrected_p, h]=bonf_holm(S2(:,2),0.05)
%
for i=1:l-5
    if corrected_p(i,1) >= 0.05
        break
    end
    d = S2(i,1);
    Final_time(i).title = c(d).title;
    Final_time(i).p_adjusted = corrected_p(i,1)
    Final_time(i).r = c(d).r;
end

%% Multiple Comparisons Correction: Bonferroni protected alpha method 
t=0;
for i=1:l-5
    if ~isempty (c(1,i).title) | ~isempty (c(2,i).title)
        if c(1,i).p < .05/210
            t= t+1;
            if ~isempty (c(1,i).title)
                C(t).trait = c(1,i).title;
            else
                C(t).trait = c(2,i).title;
            end
            C(t).R_date = c(1,i).r;
            C(t).P_date = c(1,i).p;
            C(t).R_time = c(2,i).r;
            C(t).P_time = c(2,i).p;
            C(t).N = c(2,i).n;
        end
    end
end
%% correlation between Grit and Consciousness
[r,p]=corrcoef(A(:,56),A(:,105), 'row', 'complete')

%% descriptive stats
for i=1:length(title)
    table(i).trait = title{1,i};
    table(i).Min = nanmin(A(:,i+4));
    table(i).Max = nanmax(A(:,i+4));
    table(i).Mean = nanmean(A(:,i+4));
    table(i).SD = nanstd(A(:,i+4));
    table(i).N = sum(~isnan(A(:,i+4)));
end
%% Arizona color palette
AZred = [171,5,32]/256;
AZblue = [12,35,75]/256;
AZcactus = [92, 135, 39]/256;
AZsky = [132, 210, 226]/256;
AZriver = [7, 104, 115]/256;
AZsand = [241, 158, 31]/256;
AZmesa = [183, 85, 39]/256;
AZbrick = [74, 48, 39]/256;
%%
for i=1:20
    figure(1);
    subplot(5, 4, i)
    histogram(A(:,i+4), 'Facecolor', AZsand);
    set(gca,'fontsize',15);
    xlabel(title{1,i}, 'FontSize', 15);
end
subplot(5, 4, 1)
xticklabels({'M', 'F', 'Trans'});

subplot(5, 4, 3)
xticks([1 2 3 4 5 6 7 8]);
xticklabels({'A', 'B', 'L', 'N', 'W', 'M','I', 'O'});

subplot(5, 4, 4)
xticks([0 1]);
xticklabels({'No', 'Yes'});
%%
for i=21:40
    figure(2);
    subplot(5, 4, i-20)
    histogram(A(:,i+4), 'Facecolor', AZsand);
    set(gca,'fontsize',15);
    xlabel(title{1,i}, 'FontSize', 15);
end
%%
for i=41:60
    figure(3);
    subplot(5, 4, i-40)
    histogram(A(:,i+4), 'Facecolor', AZsand);
    set(gca,'fontsize',15);
    if (i-40 == 13) | (i-40 == 14) | (i-40 == 19) | (i-40 == 20)
        xlabel(title{1,i}, 'FontSize', 10);
    else
        xlabel(title{1,i}, 'FontSize', 15);
    end
end
%%
for i=61:80
    figure(4);
    subplot(5, 4, i-60)
    histogram(A(:,i+4), 'Facecolor', AZsand);
    set(gca,'fontsize',15);
    if 2<(i-60) & (i-60)<8
        xlabel(title{1,i}, 'FontSize', 11);
    else
        xlabel(title{1,i}, 'FontSize', 15);
    end
end
%%
for i=81:105
    figure(5);
    subplot(5, 5, i-80)
    histogram(A(:,i+4), 'Facecolor', AZsand);
    set(gca,'fontsize',15);
    if (i-80) == 8 | (i-80) == 7
        xlabel(title{1,i}, 'FontSize', 9);
    else
        xlabel(title{1,i}, 'FontSize', 13);
    end
end
%%
% Bonferroni-Holm (1979) correction for multiple comparisons.  This is a
% sequentially rejective version of the simple Bonferroni correction for multiple
% comparisons and strongly controls the family-wise error rate at level alpha.
%
% It works as follows:
% 1) All p-values are sorted in order of smallest to largest. m is the
%    number p-values.
% 2) If the 1st p-value is greater than or equal to alpha/m, the procedure
%    is stopped and no p-values are significant.  Otherwise, go on.
% 3) The 1st p-value is declared significant and now the second p-value is
%    compared to alpha/(m-1). If the 2nd p-value is greater than or equal 
%    to alpha/(m-1), the procedure is stopped and no further p-values are 
%    significant.  Otherwise, go on. 
% 4) Et cetera.
%
% As stated by Holm (1979) "Except in trivial non-interesting cases the 
% sequentially rejective Bonferroni test has strictly larger probability of
% rejecting false hypotheses and thus it ought to replace the classical 
% Bonferroni test at all instants where the latter usually is applied."
%
%
% function [corrected_p, h]=bonf_holm(pvalues,alpha)
%
% Required Inputs:
%  pvalues - A vector or matrix of p-values. If pvalues is a matrix, it can
%            be of any dimensionality (e.g., 2D, 3D, etc...).
%
% Optional Input:
%  alpha   - The desired family-wise alpha level (i.e., the probability of
%            rejecting one of more null hypotheses when all null hypotheses are
%            really true). {default: 0.05}
%
% Output:
%  corrected_p  - Bonferroni-Holm adjusted p-values. Any adjusted p-values
%                 less than alpha are significant (i.e., that null hypothesis 
%                 is rejected).  The adjusted value of the smallest p-value
%                 is p*m.  The ith smallest adjusted p-value is the max of
%                 p(i)*(m-i+1) or adj_p(i-1).  Note, corrected p-values can
%                 be greater than 1.
%  h            - A binary vector or matrix of the same dimensionality as
%                 pvalues.  If the ith element of h is 1, then the ith p-value
%                 of pvalues is significant.  If the ith element of h is 0, then
%                 the ith p-value of pvalues is NOT significant.
%
% Example:
% >>p=[.56 .22 .001 .04 .01]; %five hypothetical p-values
% >>[cor_p, h]=bonf_holm(p,.05)
% cor_p =
%    0.5600    0.4400    0.0050    0.1200    0.0400
% h =
%     0     0     1     0     1
% 
% Conclusion: the third and fifth p-values are significant, but not the
% remaining three.
%
% Reference:
% Holm, S. (1979) A simple sequentially rejective multiple test procedure.
% Scandinavian Journal of Statistics. 6, 65-70.
%
%
% For a review on contemporary techniques for correcting for multiple
% comparisons that are often more powerful than Bonferroni-Holm see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis 
% of event-related brain potentials/fields I: A critical tutorial review. 
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x 
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010
% 
% 
function [corrected_p, h]=bonf_holm(pvalues,alpha)

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvalues<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvalues>1,1)),
        fprintf('WARNING: Some uncorrected p-values are greater than 1.\n');
    end
end

if nargin<2,
    alpha=.05;
elseif alpha<=0,
    error('Alpha must be greater than 0.');
elseif alpha>=1,
    error('Alpha must be less than 1.');
end

s=size(pvalues);
if isvector(pvalues),
    if size(pvalues,1)>1,
       pvalues=pvalues'; 
    end
    [sorted_p sort_ids]=sort(pvalues);    
else
    [sorted_p sort_ids]=sort(reshape(pvalues,1,prod(s)));
end
[dummy, unsort_ids]=sort(sort_ids); %indices to return sorted_p to pvalues order

m=length(sorted_p); %number of tests
mult_fac=m:-1:1;
cor_p_sorted=sorted_p.*mult_fac;
cor_p_sorted(2:m)=max([cor_p_sorted(1:m-1); cor_p_sorted(2:m)]); %Bonferroni-Holm adjusted p-value
corrected_p=cor_p_sorted(unsort_ids);
corrected_p=reshape(corrected_p,s);
h=corrected_p<alpha;
end
