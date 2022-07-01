%% Adaptation of sumBiomass.m function
% This LIVE script guides the user through all the necessary steps for adapting 
% the |sumBiomass.m| GECKO function to a _Rhodotorula toruloides_ metabolic model 
% (<https://github.com/SysBioChalmers/rhto-GEM rhto-GEM>).
%% 1.- Open the original GECKO |sumBiomass.m| function

open ../../geckomat/limit_proteins/sumBioMass.m
%% 
% This script was initially created with the purpose of computing the mass 
% fraction occupied by the main macromolecular groups (proteins, carbohydrates, 
% lipids, RNA and DNA) per gram of biomass produced by the biomass pseudoreaction 
% in <https://github.com/SysBioChalmers/yeast-GEM YeastGEM>.
% 
% This script consists on a main function body, in which all relevant biomass 
% precursors are listed, including their metabolite ID in the model, molecular 
% weight and macromolecule class. Additionally, an auxiliary nested function is 
% called, |getFraction|, and used for the mass fraction computations. 
%% 2.- Explore biomass and macromolecular synthesis pseudoreactions in rhtoGEM
% The biomass pseudoreaction, has been identified and explored previously by 
% the |explore_rhtoGEM.xlm| live script in this tutorial. In order to analyse 
% its structure, the name of the biomass reaction and the ones for production 
% of its main precursors can be retrieved from the previously created |getParameters.m| 
% script by executing it

cd rhto_scripts
parameters = getModelParameters;
cd ..
%% 
% Load |rhtoGEM| and display the biomass pseudoreaction formula

load('../models/rhtoGEM.mat')
constructEquations(model,parameters.bioRxn)
%% 
% Display the main macromolecular synthesis pseudoreactions and get a list 
% of their constituent metabolites

for macroM = [{'protein'} {'carbohydrate'} {'lipid'} {'RNA'} {'DNA'}]
    disp([macroM{1} ' precursors:'])
    macroM_pos = find(strcmp(model.metNames,macroM));
    reaction   = find(model.S(macroM_pos,:)>0);
    constructEquations(model,reaction)
    precursors = find(model.S(:,reaction)<0);
    IDs        = model.mets(precursors);
    metName    = model.metNames(precursors);
    precTable  = table(IDs,metName)
    disp(' ')
end