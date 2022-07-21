function [isoEnzymes,Promiscuous,Complexes,RxnWithKcat] =  rxnsCounter(model,originalModel,barplots)
%
% Function that counts the number of reactions with isoenzymes, enzyme 
% complexes, Promiscuous enzymes and rxns with  Kcat values in an ecModel.
%
%   model
%   originalModel
%   barPlots
%
%   isoEnzymes
%   Promiscuous
%   Complexes
%   RxnWithKcat
%
% Usage: [isoEnzymes,Promiscuous,Complexes,RxnWithKcat] =  rxnsCounter(model,originalModel,barplots)
%
% Ivan Domenzain.       Last edited: 2019-06-21

if nargin<3
    barplots = false;
end
% Get enzymes indxs
origRxns    = originalModel.rxns;
%Find enzyme pseudometabolites indexes
enzMetIndxs = find(contains(model.metNames,'prot_'));
%Find enzyme usage pseudoreaction indexes
enzRxnIndxs = find(contains(model.rxnNames,'prot_'));
enzRxnIndxs = enzRxnIndxs(3:end);
isoEnzymes  = find(contains(model.rxnNames,'(arm)'));
isoEnzymes  = length(isoEnzymes);
S_enzymes   = model.S(enzMetIndxs,1:enzRxnIndxs(1)-1);
[m, n]      = size(S_enzymes);
promDist    = [];
compDist    = [];
isoEDist    = [];
ylab        = 'Frequency';
% Isoenzymes
for i=1:length(origRxns)
    ID = origRxns{i};
    %Find isoenzyme reactions related to each original reaction
    ID = [ID 'No'];
    nIso = sum(contains(model.rxns,ID));
    if nIso>1
        isoEDist = [isoEDist; nIso];
    end
end
% Promiscuous
Promiscuous = 0;
for i=1:m
    nProm = length(find(S_enzymes(i,:)~=0));
    if nProm>1
        Promiscuous = Promiscuous+1;
        promDist    = [promDist; nProm];
    end
end
% Complexes
Complexes = 0; RxnWithKcat = 0;
for i=1:n
    nComp = length(find(S_enzymes(:,i)~=0));
    if nComp>0
        RxnWithKcat = RxnWithKcat+1;
        if nComp>1
            Complexes = Complexes+1;
            compDist  = [compDist; nComp];
        end
    end
end
%Plot results 
if barplots
    titlestr = [num2str(isoEnzymes) ' Reactions with isoenzymes'];
    plotHist(isoEDist,'# of isoenzymes per reaction',ylab,titlestr,100);
    titlestr = [num2str(Promiscuous) ' Promiscuous enzymes'];
    plotHist(promDist,'# of reactions per enzyme',ylab,titlestr,100);
    titlestr = [num2str(Complexes) ' Enzyme complexes'];
    plotHist(compDist,'# of subunits per complex',ylab,titlestr,100);
end
end
%--------------------------------------------------------------------------
function plotHist(var,xlab,ylab,titlestr,nBins)
figure
hist(var,nBins)
title(titlestr)
xlabel(xlab)
ylabel(ylab)
end

