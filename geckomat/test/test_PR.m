%load yeastGEM model file
load('yeastGEM.mat')
%enhance model with enzyme-constraints
name = 'ecYeastGEM';
toolbox = 'COBRA';
cd ..
%[ecModel,ecModel_batch] = enhanceGEM(model,toolbox,name);
cd ../change_model
[ecModel,modifications] = manualModifications(ecModel);
cd ../limit_proteins
[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,modifications,'ecYeastGEM');
cd ../test
save('first_iteration.mat','ecModel','ecModel_batch')
%find differences between curated Kcats and those in the ecModel
kcats_table = find_kcat_discrepancies(ecModel_batch);
writetable(kcats_table,'kcat_data_comparison.txt','delimiter','\t','QuoteStrings',false);
%Check the status of this repo and branch in git, no changes should be
%present in the ecModel output files (probably just white spaces)

%The ecModel is already calibrated for growth maximization FBA simulations,
%let's obtain a parsimonious flux distribution and save it as a text file
sol = solveLP(ecModel_batch,1);
%print exchange fluxes in command window
printFluxes(ecModel_batch,sol.x)
%Generate a .txt file with the whole flux distribution
subSystems = cellfun(@strjoin,ecModel_batch.subSystems,repelem({' // '},length(ecModel_batch.subSystems))','UniformOutput',false);
fluxDist   = table(ecModel_batch.rxns,ecModel_batch.rxnNames,sol.x,ecModel_batch.grRules,subSystems);
fluxDist.Properties.VariableNames = {'rxns' 'rxnNames' 'flux' 'grRules' 'subSystems'};
writetable(fluxDist,'maxGrowth_fluxDist.txt','delimiter','\t','QuoteStrings',false);
%Obtain a distribution of "net" fluxes
netFluxDist = getNetFluxes(fluxDist.rxns,fluxDist.rxnNames,fluxDist.flux,fluxDist.grRules,fluxDist.subSystems);
formulas = constructEquations(model,netFluxDist.rxns);
netFluxDist.formulas = formulas;
writetable(netFluxDist,'maxGrowth_net_fluxDist.txt','delimiter','\t','QuoteStrings',false);
%get metabolites turnover
metTurnovers = met_turnovers(netFluxDist,model);
writetable(metTurnovers,'maxGrowth_metsTurnover.txt','delimiter','\t','QuoteStrings',false);
%get flux control coefficients
target = find(model.c);
ECCs = getECCs(ecModel_batch,target,1.01);
writetable(ECCs,'maxGrowth_FCCs.txt','delimiter','\t','QuoteStrings',false);
%find intersect between enzymes with non-zero FCCs and curated info
nonZeroFCCs = ECCs(ECCs.CC>0,:);
[presence,iB] = ismember(nonZeroFCCs.genes,kcats_table.genes);
%

