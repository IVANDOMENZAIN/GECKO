function model = fillProtConcs(model, protData)
% fillProtConcs
%   Uses the protein concentrations from protData to fill model.ec.concs.
%   Protein levels should be provided in mg/gDCW. If no data is provided
%   a particular protein, its level is NaN. Existing entries in
%   model.ec.concs are overwritten.
%
% Input:
%   model       an ecModel in GECKO 3 format (with ecModel.ec structure)
%   protData    structure with proteome data, from loadProtFluxData
%               uniprotIDs      cell arrray with Uniprot IDs matching
%                               protData.abundances
%               abundances      matrix of proteomics data
%
% Output:
%   model       an ecModel where model.ec.concs is populated with protein
%               concentrations.
%
% Note: to also constrain the model with the content of model.ec.concs, you
% should run constrainProtConcs.
%
% Usage:
%   model = fillProtConcs(model, protData)

uniprotIDs = protData.uniprotIDs;
abundances = protData.abundances;

%Redefine an empty model.ec.concs vector
model.ec.concs=nan(numel(model.ec.enzymes),1);

[a,b] = ismember(uniprotIDs, model.ec.enzymes);
model.ec.concs(b(a)) = abundances(a);
end
