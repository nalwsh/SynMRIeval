% Evaluation of synthetic data with statistical similarity
% Evaluation methods : 1) Kolmogorov-Smirnov test, 2) Kullback–Leibler
% Divergence 3) Jensen-Shannon divergence
% Example data set : simulated synthetic data. v.s real data
% Dimension of data set : 1 X 118 vector
clearvars;
%%% User specifiecd variables here %%%
% load example data set
x1 = load('SynMMSE.mat'); x2 = load('RealMMSE.mat');
%%% data convert
d1 = cell2mat(struct2cell(x1)); d2 = cell2mat(struct2cell(x2));
%%% 2 sample KS-test
[h,p,ks2stat] = kstest2(d1, d2);
if h == 0
    sprintf('The null hypothesis that the data of vectors x1 and x2 are extracted from the same continuous distribution was not rejected, and p value is %f', p)
elseif h == 1
    sprintf('The null hypothesis that the data of vectors x1 and x2 are extracted from the same continuous distribution was rejected, and p value is %f', p)
end
%%% Kullback-Leibler Divergence
h1 = histogram(d1, 'NumBins', 10, 'BinWidth', 2, 'BinEdges', [min(d1):max(d1)],'FaceAlpha',0.3, 'EdgeAlpha',0.1, 'FaceColor','r');
h1.Normalization = 'probability';
h1Data = h1.Data; h1Count = h1.Values; h1Pr = h1Count./sum(h1Count); h1Bin = h1.BinEdges;
hold on;
h2 = histogram(d2, 'NumBins', 10, 'BinWidth', 2, 'BinEdges', [min(d1):max(d1)],'FaceAlpha',0.3, 'EdgeAlpha',0.1, 'FaceColor','b');
h2.Normalization = 'probability';
h2Data = h2.Data; h2Count = h2.Values; h2Pr = h2Count./sum(h2Count); h2Bin = h2.BinEdges;
hold off;
h1Pr = h1Pr + eps; h2Pr = h2Pr + eps;
h1Bin(1) = [];
KL = kldiv(h1Bin, h1Pr, h2Pr);
sprintf('Kullback–Leibler Divergence is %f\n', KL)
KLjs = kldiv(h1Bin, h1Pr, h2Pr, 'js');
sprintf('Jensen-Shannon divergence is %f\n', KLjs)
%%% The function below is from David Fass (2023). KLDIV (https://www.mathworks.com/matlabcentral/fileexchange/13089-kldiv)
function KL = kldiv(varValue,pVect1,pVect2,varargin)
if ~isequal(unique(varValue),sort(varValue)),
    warning('KLDIV:duplicates','X contains duplicate values. Treated as distinct values.')
end
if ~isequal(size(varValue),size(pVect1)) || ~isequal(size(varValue),size(pVect2)),
    error('All inputs must have same dimension.')
end
% Check probabilities sum to 1:
if (abs(sum(pVect1) - 1) > .00001) || (abs(sum(pVect2) - 1) > .00001),
    error('Probablities don''t sum to 1.')
end

if ~isempty(varargin),
    switch varargin{1},
        case 'js',
            logQvect = log2((pVect2+pVect1)/2);
            KL = .5 * (sum(pVect1.*(log2(pVect1)-logQvect)) + ...
                sum(pVect2.*(log2(pVect2)-logQvect)));

        case 'sym',
            KL1 = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
            KL2 = sum(pVect2 .* (log2(pVect2)-log2(pVect1)));
            KL = (KL1+KL2)/2;

        otherwise
            error(['Last argument' ' "' varargin{1} '" ' 'not recognized.'])
    end
else
    KL = sum(pVect1 .* (log2(pVect1)-log2(pVect2)));
end
end