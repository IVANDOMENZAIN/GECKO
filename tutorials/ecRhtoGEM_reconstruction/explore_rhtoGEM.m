%% *Explore rhtoGEM*
%% 1.- Download GECKO and rhtoGEM repositories
% In order to enhance the *rhtoGEM* with enzyme constraints, the original source 
% repository and a copy of are downloaded from GitHub.
% 
% clone GECKO repository

git('clone https://github.com/SysBioChalmers/GECKO')

% clone model repository

git('clone https://github.com/SysBioChalmers/rhto-GEM')
clc
%% 
% Load rhtoGEM model file:

model = load('rhto-GEM/ModelFiles/MAT/rhto.mat');
model = model.model;
mkdir('../models')
%% 
% Save original GEM file in our repository
% 
% %Probably unnecessary

save('../models/rhtoGEM.mat','model')
%% 2.-Explore rhtoGEM - adapt getModelParameters.m
% Open the |getModelParameters.m script available in the GECKO toolbox|

open GECKO/geckomat/getModelParameters.m
%% 
% This script contains model-specific information that is retrieved by the 
% different steps of the model enhancement pipeline in |GECKO/geckomat, |including:
% 
% * physiological parameters for _S. cerevisiae_
% * identifiers and names for reactions related to biomass production, glucose 
% * uptake, oxygen uptake, CO2 production, energy maitenance, etc.
% * Polimerization cost for the main macromolecules (mol of ATP per unit of 
% product).
% * Reaction IDs for the reactions involved in the oxidative phosphorylation 
% pathway
% 
% As all information in this script refers to  <https://github.com/SysBioChalmers/yeast-GEM 
% yeastGEM> and _S. cerevisiae_,  the *rhtoGEM* structure needs to explored in 
% order to create a new version of the |getModelParameters.m| script for _Rhodotorula 
% toruloides_.
%% *2.1 Add organism-specific information*
% Substitute the following variables with the indicated values:
% 
% * |parameters.org_name = 'rhodotorula toruloides';|
% * |parameters.Ptot = 0.2669; (Tiukova, et al. 2019)|
% * |parameters.gR_exp = 0.3; (|<https://link.springer.com/content/pdf/10.1007/s11274-011-0920-2.pdf 
% Andrade, et al. 2012>|)|
% 
% As the yeast _R. toruloides_ is not present in the <https://www.genome.jp/kegg/genome/ 
% KEGG genome> database, provide an empty string to the KEGG ID variable to optimize 
% GECKO running:
% 
% * |parameters.keggID = '';|
%% 2.2 Extract biomass reaction data
% Find the biomass pseudoreaction in the model and obtain its full name and 
% unique reaction identifier:

bioRxns = find(contains(model.rxnNames,'biomass'));
if ~isempty(bioRxns)
    disp(model.rxnNames(bioRxns))
    disp(model.rxns(bioRxns))
end
%% 
% Display biomass reaction formula, using the <https://github.com/SysBioChalmers/RAVEN 
% RAVEN> |constructEquations built-in function:|

constructEquations(model,bioRxns)
%% 
% Retrieve biomass precursors

prec_Ids   = model.mets(find(model.S(:,bioRxns)<0));
prec_names = model.metNames(find(model.S(:,bioRxns)<0));
%% 
% Display biomass precursors

disp(prec_names)
%% 
% Display biomass reaction products

prodIDs   = model.mets(find(model.S(:,bioRxns)>0));
prodNames = model.metNames(find(model.S(:,bioRxns)>0));
disp(prodNames)
%% 
% The biomass reaction is a modular reaction, therefore individual synthesis 
% reactions for each of the main macromolecular groups are also present in the 
% model (proteins, carbohydrates, lipids, RNA and DNA). The GECKO toolbox rescales 
% the coefficients in the biomass pseudoreaction for the most abundant macromolecular 
% groups, according to the provided cellular protein content. It is then needed 
% to find the names for these biosynthetic pseudoreactions.

macromolecules = {'protein' 'carbohydrate' 'lipid'};
for macro=macromolecules
    index = find(strcmp(model.metNames,macro));
    rxns  = find(model.S(index,:)>0);
    disp(['The synthesis reaction for ' macro{1} ' is called: ' model.rxnNames{rxns}])
end
%% 
% The name for the lipid pseudoreaction (|lipid pseudoreaction - merge|) 
% suggests that this is a <https://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-018-0673-8 
% SLIMEr>-type reaction, where the synthesis of lipid backbones and fatty acid 
% chains are synthesized separately. Find the corresponding reaction names for 
% each of these, by retrieving the reactions that synthesize the substrates of 
% the |lipid |pseudometabolite

lipid     = find(strcmpi(model.metNames,'lipid'));
lipRxn    = find(model.S(lipid,:)>0);
lip_parts = find(model.S(:,lipRxn)<0);
for i=1:numel(lip_parts)
    lipRxn = find(model.S(lip_parts(i),:)>0);
    disp(['Reaction name found: ' model.rxnNames{lipRxn}])
end
%% 
% Now that all the biomass pseudoreaction macromolecular components have 
% been identified the correspoding variables should be updated in the |getModelParameters.m| 
% script as follows:
% 
% * |parameters.bio_comp{1} = 'protein';|
% * |parameters.bio_comp{2} = 'carbohydrate';|
% * |parameters.bio_comp{3} = 'lipid backbone';|
% * |parameters.bio_comp{4} = 'lipid chain';|
%% 2.3 Assign enzymes to a single cellular compartment
% For simplicity, GECKO incorporates all enzymes as pseudometabolites located 
% in the same cellular compartment, it is recommended to locate them in the cytoplasmic 
% compartment. In order to provide a compartment name that complies with the compartment 
% names in the model, display the content of the |model.compNames |field and modify 
% the variable |parameters.enzyme_comp |accordingly.

disp(model.compNames)
%% 
% * |parameters.enzyme_comp = 'cytoplasm';|
%% 2.4 Search main exchange reactions
% In order to reconstruct a functional ecModel suitable for simulations, the 
% growth-associated maintenance energy expenditure (GAM) is fitted in the raw 
% ecModel using experimental data on exchange fluxes for glucose, oxygen and carbon 
% dioxide at two, or more, growth rates, within the purely respiratory regime 
% of the cell. As the GECKO pipeline performs this parameter fitting automatically, 
% it is necessary to provide the name of the relevant exchange flux reactions.
% 
% To identify all exchange reactions in the model use the built-in <https://github.com/SysBioChalmers/RAVEN 
% RAVEN> function  |getExchangeRxns|:

[~, positions] = getExchangeRxns(model);
%% 
% Each exchange reaction comprises a single compound and ussually take the 
% following form:
% 
% * |compound[compartment] ->|
% 
% In order to find the desired exchange reactions, retrieve the metabolite 
% that is exchanged in each of the identified reactions in the previous step and 
% identify those that are specific for D-glucose, oxygen, carbon dioxide and biomass 
% (cellular growth)

exchCompounds = {'biomass' 'D-glucose' 'oxygen' 'carbon dioxide'};
for i=1:numel(positions)
    metPos  = find(model.S(:,positions(i)));
    metName = model.metNames(metPos);
    if ismember(metName,exchCompounds)
        disp(model.rxnNames(positions(i)))
    end
end
%% 
% Update the correspoding variables in the parameters script:
% 
% * |parameters.exch_names{1} = 'growth';|
% * |parameters.exch_names{2} = 'D-glucose exchange (reversible)';*|
% * |parameters.exch_names{3} = 'oxygen exchange (reversible)';*|
% * |parameters.exch_names{4} = 'carbon dioxide exchange';|
% 
% * Note that the substring |' (reversible)'| has been added to the glucose 
% and oxygen exchange reactions. This is due to the fact that all unbounded exchange 
% reactions in the original model are split into an  uptake and a secretion reaction 
% for the same exchanged metabolite. The exchange reaction formula above illustrates 
% that, ussually, exchange reactions are formulated as secretions in the model 
% estructure, therefore the  uptake reactions created by GECKO are named as "|compound 
% exchange (reversible)"|. For the purpose of GAM fitting, secretion directionality 
% is required for biomass and carbon dioxide reactions, while  uptakes are required 
% for the oxygen and D-glucose cases.
%% Get total protein content from biomass reaction
% Get list of precursors for protein pseudoreaction

protRxn  = find(strcmpi(model.rxnNames,'protein pseudoreaction'))
protPrec = find(model.S(:,protRxn)<0);
disp(model.mets(protPrec))
%% 3. Provide parameters for the proteomics integration utility
% The GECKO pipeline returns an enzyme-constrained model version of the original 
% provided GEM, in two different files |ec_modelName_batch.mat |and |ec_modelName.mat|, 
% the former being an ecModel suitable for simulations, and constrained by a total 
% pool of enzyme that is proportional to the total protein content of the cell; 
% the latter being an ecModel with unconstrained enzyme usage reactions, suitable 
% for integration of proteomics data as upper bounds for these reactions. A <https://github.com/SysBioChalmers/GECKO/tree/master/geckomat/utilities 
% pipeline for proteomics data integration> is available in the GECKO utilities 
% folder.
% 
% Proteomics data are one of the several high-throughput omics data types 
% available for cell physiology. Proteomics measurements tend to display high 
% variability across, both technical and biological replicates (<https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/pmic.202000093 
% Sánchez et al. 2021>)
%% 3.1 .- Find NGAM pseudoreaction
% Cells require a free energy supply (ATP) not just for reproduction but also 
% for maintaining all the bioprocesses involved in homeostasis, this is represented 
% in a GEM as a "non-growth associated maintenance" pseudoreaction (NGAM) which 
% consumes intracellular ATP. When proteomics data are provided as constraints, 
% GECKO flexibilizes individual measurements that hamper simulation of biologically 
% meaningful phenotypes, for this step it is necessary that experimentally-derived 
% bounds are imposed on this pseudoreaction. Therefore, the |getModelParameters.m| 
% script requires the name of such reaction for easy accession in the <https://github.com/SysBioChalmers/GECKO/tree/master/geckomat/utilities/integrate_proteomics 
% GECKO data integration pipeline>.
% 
% Try to find the NGAM reaction by searching different similar strings in 
% the model.rxnNames field:


for str = {'NGAM' 'non-growth' 'maintenance'}
    NGAM = find(contains(model.rxnNames,str));
    if ~isempty(NGAM)
        disp(model.rxnNames(NGAM))
        disp(model.rxns(NGAM))
    end
end
%% 
% If successful, then substitute the corresponding variable in the |getModelParameters.m| 
% script:
% 
% * |parameters.NGAM = 'r_4046';|
%% 2.4 Search oxidative phosphorylation reactions
% The associated metabolic pathways for each reaction are stored in the model 
% field "subSystems"

oxPhos = [];
for i =1:length(model.subSystems)
    rxnSubSyst = model.subSystems{i};
    if any(contains(rxnSubSyst,'Oxidative phosphorylation'))
        oxPhos = [oxPhos; i];
    end
end
disp(model.rxnNames(oxPhos))
disp(' ')
disp(model.rxns(oxPhos))