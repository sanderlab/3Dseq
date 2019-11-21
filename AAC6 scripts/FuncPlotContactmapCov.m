function [CorrectMonomerContacts, CorrectDimerContacts] = FuncPlotContactmapCov(CovPairs,MonomerContactList,L,varargin)
%FuncPlotContactmapCov plots the contact map and ECs of pairs of residues in a protein.
%   FuncPlotCintactmapCov(CovPairs,MonomerContactList,L,'Lcutoff',100,...) 
%
%INPUTS
%   CovPairs - ordered m-by-3 list of residues and their covariance values (or ECs)
%   MonomerContactList - m-by-2 list of residues contacting in the monomer
%   L - length of protein (number of residues)
%
%OPTIONAL
%
%       dimercontacts - m-by-2 list of residues contacting in the dimer
%       PDistCutoff - primary sequence distance truncation
%       LCutoff - truncation to first Lcutoff CovPairs (after applying a potential PDistCutoff)
%       MSize - plot marker size, default is 20
%
%
%DEPENDENCIES: 
%
%Frank J. Poelwijk (c) 02/14/2018
%% Initialization
MSize=20;


%% Optional Parameters
nvarargin = numel(varargin);
if nvarargin > 0
    if rem(nvarargin,2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'dimercontacts' 'PDistCutoff' 'LCutoff' 'MSize'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), lower(okargs)); 
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % contacts in dimer
                    DimerContactList=pval;
                case 2  % primary sequence distance truncation
                    PDistCutoff = pval;
                case 3  % truncation to first Lcutoff CovPairs
                    LCutoff = pval;
                case 4  % marker size
                    MSize = pval;
            end
        end
    end
end

if exist('PDistCutoff','var')
    resprimdistCovPairs=abs(CovPairs(:,1)-CovPairs(:,2));
    CovPairsPTrunc=CovPairs(resprimdistCovPairs>PDistCutoff,:);
    resprimdistMono=abs(MonomerContactList(:,1)-MonomerContactList(:,2));
    MonomerContactListTrunc=MonomerContactList(resprimdistMono>PDistCutoff,:);
else
    CovPairsPTrunc=CovPairs;
    MonomerContactListTrunc=MonomerContactList;
end
        
if exist('LCutoff','var')
    LCutoff=ceil(LCutoff);
    CovPairsPLTrunc=CovPairsPTrunc(1:LCutoff,:);
else
    CovPairsPLTrunc=CovPairsPTrunc;
end

CorrectMonomerContacts=intersect(CovPairsPLTrunc(:,1:2),MonomerContactListTrunc(:,1:2),'rows');
if exist('DimerContactList','var')
    CorrectDimerContacts=intersect(CovPairsPLTrunc(:,1:2),DimerContactList(:,1:2),'rows');
end

plot([1,L],[1,1],'k') % QnD way to get the axes on top and right...
hold all
scatter(MonomerContactList(:,1),MonomerContactList(:,2),MSize,[0.9 0.9 0.9],'filled');
scatter(MonomerContactList(:,2),MonomerContactList(:,1),MSize,[0.9 0.9 0.9],'filled');
scatter(MonomerContactListTrunc(:,1),MonomerContactListTrunc(:,2),2*MSize,[0.7 0.7 0.7],'filled');
scatter(MonomerContactListTrunc(:,2),MonomerContactListTrunc(:,1),2*MSize,[0.7 0.7 0.7],'filled');
if exist('DimerContactList','var')
    scatter(DimerContactList(:,1),DimerContactList(:,2),2*MSize,[0.7 0.7 0.9],'filled');
    scatter(DimerContactList(:,2),DimerContactList(:,1),2*MSize,[0.7 0.7 0.9],'filled');
end
scatter(CovPairsPLTrunc(:,1),CovPairsPLTrunc(:,2),MSize,[0 0 0],'filled');
scatter(CovPairsPLTrunc(:,2),CovPairsPLTrunc(:,1),MSize,[0 0 0],'filled');
scatter(CorrectMonomerContacts(:,1),CorrectMonomerContacts(:,2),2*MSize,[1 0 0],'filled');
scatter(CorrectMonomerContacts(:,2),CorrectMonomerContacts(:,1),2*MSize,[1 0 0],'filled');
if exist('DimerContactList','var')
    scatter(CorrectDimerContacts(:,1),CorrectDimerContacts(:,2),2*MSize,[0.3 0.3 1],'filled');
    scatter(CorrectDimerContacts(:,2),CorrectDimerContacts(:,1),2*MSize,[0.3 0.3 1],'filled');
end
axis([1 L 1 L])
axis square
axis ij
hold off

