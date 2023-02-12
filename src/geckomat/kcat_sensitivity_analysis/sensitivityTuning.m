function [model, tunedKcats] = sensitivityTuning(model, desiredGrowthRate, modelAdapter, foldChange)
% sensitivityTuning
%    Function that relaxes the most limiting kcats until a certain growth rate 
%    is reached. The function will update kcats in model.ec.kcat.
%    
% Input:
%   model              an ecModel in GECKO 3 version
%   desiredGrowthRate  kcats will be relaxed until this growth rate is reached
%   modelAdapter       a loaded model adapter (Optional, will otherwise use the
%                      default model adapter).
%   foldChange         kcat values will be increased by this fold-change.
%                      (Opt, default 10)
% Output:
%   model              ecModel with updated model.ec.kcat
%
% Usage:
%   model = sensitivityTuning(model, 0.07, humanGEMAdapter)

if nargin < 4 || isempty(foldChange)
    foldChange = 10;
end
if nargin < 3 || isempty(modelAdapter)
    modelAdapter = ModelAdapterManager.getDefaultAdapter();
    if isempty(modelAdapter)
        error('Either send in a modelAdapter or set the default model adapter in the ModelAdapterManager.')
    end
end
params = modelAdapter.params;
if nargin < 2 || isempty(desiredGrowthRate)
    desiredGrowthRate = params.gR_exp;
end

kcatList = [];
m = model;
m.c = double(strcmp(m.rxns, params.bioRxn));% Make sure that growth is maximized

%To avoid numerical issues, make sure no kcat is below 0.1
%This is not desirable to do, but unfortunately necessary - otherwise the 
%solver says the solution is infeasible
%m.ec.kcat(m.ec.kcat < 0.1) = 0.1;
%m = applyKcatConstraints(m);

if ~m.ec.geckoLight
    %for the full model, we first find the draw reaction with the most flux
    drawRxns = startsWith(m.rxns, 'usage_prot_');
    iteration = 1;
    while true
        res = solveLP(m,0); %skip parsimonius, only takes time
        lastGrowth = -res.f;
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]) 
        iteration = iteration + 1;
        %find the highest draw_prot rxn flux
        drawFluxes = zeros(length(drawRxns),1);
        drawFluxes(drawRxns) = res.x(drawRxns);
        [~,sel] = max(drawFluxes);
        %Now get the metabolite
        metSel = m.S(:,sel) > 0;
        %now find the reaction with the largest consumption of this protein
        protFluxes = m.S(metSel,:).' .* res.x; %negative
        [~,rxnSel] = min(protFluxes);
        kcatList = [kcatList, rxnSel];
        rxn = m.rxns(rxnSel);
        targetSubRxn = strcmp(m.ec.rxns, rxn);
        m.ec.kcat(targetSubRxn) = m.ec.kcat(targetSubRxn) .* foldChange;
        m = applyKcatConstraints(m,targetSubRxn);
    end
    
else
    origRxns = extractAfter(m.ec.rxns,4);
    iteration = 1;
    while true
        res = solveLP(m,0); %skip parsimonius, only takes time
        lastGrowth = -res.f;
        if (lastGrowth >= desiredGrowthRate)
            break;
        end
        %If you get an error here, it is likely due to numerical issues in the solver
        %The trick where we don't allow low kcats is to fix that, but maybe
        %it is not enough.
        disp(['Iteration ' num2str(iteration) ': Growth: ' num2str(lastGrowth)]) 
        iteration = iteration + 1;
        %find the highest protein usage flux
        protPoolStoich = m.S(strcmp(m.mets, 'prot_pool'),:).';
        [~,sel] = min(res.x .* protPoolStoich); %max consumption
        kcatList = [kcatList, find(sel)];
        rxn = m.rxns(sel.');
        targetSubRxns = strcmp(origRxns, rxn);
        m.ec.kcat(targetSubRxns) = m.ec.kcat(targetSubRxns) .* foldChange;
        m = applyKcatConstraints(m,targetSubRxns);
    end
end

kcatList        = unique(kcatList);
tunedKcats.rxns = m.rxns(kcatList);
tunedKcats.rxnNames = m.rxnNames(kcatList);
[~, rxnIdx]     = ismember(tunedKcats.rxns,m.ec.rxns);
tunedKcats.enzymes = cell(numel(kcatList),1);
for i=1:numel(rxnIdx)
    [~, metIdx]     = find(m.ec.rxnEnzMat(rxnIdx(i),:));
    tunedKcats.enzymes{i}   = strjoin(m.ec.enzymes(metIdx),';');
end
tunedKcats.oldKcat = model.ec.kcat(rxnIdx);
tunedKcats.newKcat = m.ec.kcat(rxnIdx);

model = m;

end