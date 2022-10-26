% Hazel Stewart SARS-CoV-2 3C Interactions
clear
clc

dat = struct();
dat.input = readtable('D:\Collab\Hazel_2021_221026\combined\txt\proteinGroups.txt');

% Experimental design:
% Exp 1: L = 3c, M = control
% Exp 2: M = 3c, H = control
% Exp 3: H = 3c, L = control

% Note: intensity determination setting in MaxQuant was selected as 'Total Sum' rather than the default
% which is 'Maximum'.

%% Data cleanup

% Remove reverse
dat.input.Reverse = categorical(dat.input.Reverse);
dat.input = dat.input(dat.input.Reverse ~= '+',:);

% Remove contaminants
dat.input.PotentialContaminant = categorical(dat.input.PotentialContaminant);
dat.input = dat.input(dat.input.PotentialContaminant ~= '+',:);

% Create table
% Table in format Mock 1-3, IP 1-3
dat.dat = [dat.input.IntensityM1, dat.input.IntensityH2, dat.input.IntensityL3...
    , dat.input.IntensityL1, dat.input.IntensityM2, dat.input.IntensityH3];
% %%
% temp = [dat.input.IntensityL1, dat.input.IntensityM1, dat.input.IntensityM2, dat.input.IntensityH2, dat.input.IntensityL3, dat.input.IntensityH3];
% temp(temp == 0) = NaN;
% 
% figure
% subplot(1,2,1)
% bar(sum(isnan(temp)))
% 
% subplot(1,2,2)
% bar(nanmedian(temp))

% Convert 0 to NaN;
dat.dat(dat.dat == 0) = NaN;

% Normalising for loading
% Subset those which are common to all runs
% Take the median of those.

% Identify rows with no missing data
numNaN = sum(isnan(dat.dat),2); 
notNaN = numNaN == 0;
%%
% Take the median
channelMed = nanmedian(dat.dat(notNaN,:));

% Divide all by median to control for sample loading
dat.dat ./ channelMed;

% Remove rows with more than 3 missing datapoints.
logNaN = numNaN > 3;
dat.dat = dat.dat(~logNaN,:);
dat.input = dat.input(~logNaN,:);

% Imputing missing data.
% Record missing data
dat.NaN = isnan(dat.dat);

% Knn imputation
dat.dat = knnimpute(dat.dat);

% Per sample normalisation: 3C/Mock
dat.ratio = dat.dat(:,[4:6]) ./ dat.dat(:,[1:3]);
%dat.ratio = log2(dat.dat(:,[4:6])) ./ log2(dat.dat(:,[1:3]));

% Re-normalise again, this time on complete dataset
dat.ratio = dat.ratio ./ median(dat.ratio);

%% Further normalisation?

[h,p] = ttest(dat.ratio');

% Check using Storeys approach for multiple hypothesis testing
[fdr,q] = mafdr(p);

% Note - when you investigate the above, all have a FDR < 0.05 and Q-value < 0.05.

%%
% Work out p-value and cutoff hits
phit = p < 0.05;
rhit = mean(log2(dat.ratio')) > 0.5;
hit = phit + rhit == 2;




figure
fa = scatter(mean(log2(dat.ratio)'),-log10(p),'filled','k','MarkerFaceAlpha',0.15)
hold on
% Add lines to indicate cutoffs
% P-value 0.05
line([-6,6],[-log10(0.05),-log10(0.05)],'Color','k','LineStyle','--')
% +/- log2 0.5
line([-0.5,-0.5],[0,3],'Color','k','LineStyle','--')
line([0.5,0.5],[0,3],'Color','k','LineStyle','--')

% PGAM5
fb = scatter(mean(log2(dat.ratio(hit,:))'),-log10(p(hit)),'filled','b')
fc = scatter(mean(log2(dat.ratio(346,:))'),-log10(p(346)),'filled','r')

subs = [fa(1),fb(1),fc(1)]'

legend(subs,'Not Significant','Significant Proteins','PGAM5')
xlabel('Mean Log_2 SILAC Ratio (HA-3c/HA-only)')
ylabel('-Log_1_0 p-value')

set(gca,'FontSize',15)

%% Reformat tables for export (all, and hits)
dat.ratioT = array2table(dat.ratio,'VariableNames',{'Exp1_3CoverMock','Exp2_3CoverMock','Exp3_3CoverMock'})

dat.fullExport = [dat.input , dat.ratioT];
% Export full - supplemental table 1
writetable(dat.fullExport,'D:\Collab\Hazel_2021_221026\fullExport.csv')

% Export hits - supplemental table 2
dat.hits = dat.fullExport(hit,:);
writetable(dat.hits,'D:\Collab\Hazel_2021_221026\hits.csv')

%writetable(dat.input(hit),'D:\Collab\Hazel_2021\combined\txt\hits.csv')