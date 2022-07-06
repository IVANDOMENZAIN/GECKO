%% Find essential growth requirements for rhtoGEM
% 1.- Load |rhto-GEM |into the MATLAB working space

load('models/rhtoGEM.mat')
%% 
% 2.- Use the RAVEN function getExchangeRxns in order to find all reactions 
% for uptake of metabolites from the environment

[exch_rxns,indexes] = getExchangeRxns(model,'in');
disp([num2str(numel(indexes)) ' uptake reactions were found'])
%% 
% No uptake reactions (of the form {ø} -> metabolite[cell compartment]) 
% were found in the model. This suggests that exchange reactions are defined as 
% in the opposite direction, as secretion steps. 
% 
% 3.- Get all the exchange reactions in the model

[exch_rxns,indexes] = getExchangeRxns(model);
disp([num2str(numel(indexes)) ' exchange reactions were found'])
%% 
% 4.- In order to confirm the directionality of these reactions print the 
% reaction formulas of a random set of exchanges

constructEquations(model,indexes(randi(numel(indexes),3)))
%% 
% As supposed, it seems that all reactions are defined as secretion steps. 
% This means that a negative flux represents uptake of nutrients by the cell and, 
% positive fluxes indicate secretions.
% 
% 5.- Unconstrain all exchange reactions in both directions with an arbitrarily 
% high flux value

model.ub(indexes) = 1000;
model.lb(indexes) = -1000;
%% 
% 6.- Set production of biomass as an objective to maximize

growthRxn = find(strcmpi(model.rxnNames,'growth'));
model     = setParam(model,'obj',growthRxn,1);
model     = setParam(model,'lb',growthRxn,0);
%% 
% 7.- Verify biomass production by solving the FBA problem

sol = solveLP(model);
disp(sol.x(growthRxn))
%% 
% 8.- Find all essential components for growth by blocking uptake reactions 
% one by one and verifying growth

essential = zeros(numel(indexes),1);
for i = 1:numel(indexes)
    tempModel = setParam(model,'lb',indexes(i),0);
    sol       = solveLP(tempModel);
    if sol.x(growthRxn)<=1E-3
        essential(i) = 1;
    end
end
essential = indexes(logical(essential));
rxns      = model.rxnNames(essential);
rxnNames  = model.rxns(essential);
essential = table(rxns,rxnNames);
disp(essential)
%% 
% 9.- Identify other relevant exchange flux reactions

for reaction = parameters.exch_names
    reaction = strrep(reaction{1},' (reversible)','');
    rxn = find(strcmp(model.rxnNames,reaction));
    disp([model.rxns{rxn} ': ' reaction])
end