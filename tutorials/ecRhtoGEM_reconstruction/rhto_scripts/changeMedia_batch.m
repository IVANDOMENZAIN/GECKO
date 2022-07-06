function model = changeMedia_batch(model,c_source,flux)
%changeMedia_batch
%   function that modifies the ecModel and makes it suitable for batch growth
%   simulations on different carbon sources.
%
%   model       (struct) An enzyme constrained model
%   c_source    (string) Rxn name for the main carbon source uptake reaction
%   flux        (doule) Experimental flux value [mmol/gDw h] for the main 
%               carbon source uptake reaction.
%
%
%   Usage: model = changeMedia_batch(model,c_source,flux)
%
% Ivan Domenzain        2022-07-06

% Give the carbon source (c_source) input variable with the following
% format: c_source  = 'D-glucose exchange (reversible)'
if nargin<3
    flux = 1000;
end
%first block any uptake
[rxnIDs,exchange]  = getExchangeRxns(model);
exchange           = exchange(find(contains(rxnIDs,'_REV')));
model.ub(exchange) = 0;
%Allow main carbon source uptake
c_id  = model.rxns(strcmp(model.rxnNames,c_source));
model = setParam(model,'ub',c_id,flux);
model = setParam(model, 'ub', 'r_1861_REV', +1000); % 'iron(2+) exchange';
model = setParam(model, 'ub', 'r_2060_REV', +1000); % 'sulphate exchange';
model = setParam(model, 'ub', 'r_1714', 0);         % block D-glucose production';
model = setParam(model, 'ub', 'r_1992', 0);         % block oxygen production';
model = setParam(model, 'ub', 'r_2111', +1000);     % allow biomass production
end