# Reconstruction of an enzyme-constrained model for *Rhodotorula toruloides*

This is an interactive tutorial for the enhancement of [rhto-GEM](https://onlinelibrary.wiley.com/doi/abs/10.1002/bit.27162https://onlinelibrary.wiley.com/doi/abs/10.1002/bit.27162), a genome-scale metabolic model for the lipid accumulating yeast *R. toruloides*, with enzyme constraints using the GECKO toolbox. Even though **rhto-GEM** was reconstructed using [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM) as a template, some of the GECKO functions still need to be adapted to this particular case, due to their dependence in the specific structure and nomenclature of some reactions in the original GEM.

The process for generating a new ecModel can be summarized in the following basic steps:

- Supplying of organism specific data files
- Modification of specific GECKO functions
- Enhance model with enzyme constraints by running `geckomat/get_enzyme_data/updateDatabases.m` and `geckomat/enhanceGEM.m`

## 1. Supplying of organism specific data files

The README file, available in the [main page](https://github.com/SysBioChalmers/GECKO) of the GECKO toolbox repository, lists the of data files that should be updated with information regarding the modelled organism. All the updated files should be substituted in the indicated GECKO subfolders. The list below describes each of such files and details the process to adapt them to the rhto-GEM case.

- `databases/prot_abundance.txt`: integrated protein abundance data from Pax-DB. This file is used to compute the fraction of the total cellular proteome that is actually captured by the created ecModel, in order to impose a constraint on the total enzyme pool that is proportional to a measured total protein content. GECKO offers two ways of dealing with organisms that are not present in this database, by either providing a custom relative proteomics dataset file or by deleting the prot_abundance.txt, which enables the assumption that 50% of the total proteome mass is accounted for as metabolic enzymes in the model. As R. toruloides is not present in Pax-DB, deletion of the prot_abundance.txt file is recommended for simplicity.

-	`databases/uniprot.tab`: Data connecting genes to proteins from the uniprot database1. This tab separated file should be manually downloaded from [www.uniprot.org](https://www.uniprot.org/). Searching for the specific scientific name of the desired organism (*R. toruloides*, in this case) will display all protein-gene relations available for it. Once the data is displayed in the uniprot website, the columns “protein names”, “Gene names”, “EC number” and “Sequence” must be selected,  in this specific order, from the “customize columns” tool in the website. The resulting data table should be downloaded as a tab separated file and stored in the local GECKO/databases directory under the name “uniprot.tab”. 

-	`databases/chemostatData.tsv`: File consisting of a table with the main exchange fluxes (D-glucose, O2 and CO2), obtained from chemostat cultures at different dilution rates. Fluxes must be reported in units of mmol/gDw h and dilution rates in h-1. These data are used by GECKO in order to fit the growth-associated maintenance requirements in the created ecModel. For R. toruloides, D-glucose specific uptake rates at varying dilution rates have been measured and reported by a previous study2, here we extracted these values and provided them as a tab-separated tsv file under the required file name. Introducing a chemostatData.tsv file with a different structure, or with a different number of reported exchange fluxes may require custom changes into the fitGAM.m script in the GECKO toolbox.

-	`databases/manual_data.txt`: A flat text file with two columns, consisting of protein identifiers and their corresponding, manually curated, kcat values. This is an optional requirement, as it is expected that in a first iteration of ecModel reconstruction the enzyme parameters that need curation, in order to simulate meaningful phenotypes, are not known a priori. For the enhancement of rhto-GEM no manual_data.txt file has been provided, as organism-specific reported kcat values are scarce in the literature and specialized databases.

## 2. Modification of specific GECKO functions

The GECKO pipeline is designed with the purpose of returning a functional ecModel, suited for numerical simulation, that can be used to compute biologically meaningful phenotypes. To do so, growth-associated maintenance requirements are fitted according to experimental data and cellular growth on a specified media is verified in silico. As GEMs vary drastically in their formulation of the biomass pseudoreaction (lack of standardized formats and identifiers), streamlining of fitting and verification steps cannot be completely ensured, therefore, the GECKO scripts related to these steps are needed to be adapted, considering the specific biomass pseudoreaction and the organism nutrient requirements, in order to prevent execution errors. A description of the model-specific functions in the GECKO pipeline, and the procedure followed to adapt them to the rhto-GEM case is listed below:

-	`geckomat/getModelParameters.m`: Function that saves organism specific data, physiological parameters, and model-specific identifiers for a handful of relevant reactions into a MATLAB structure for easy accession and code simplification along the pipeline execution. An extensive guide on adaptation of this script to the rhto-GEM case is explained in the MATLAB live script [rhtoGEM_exploration.mlx](https://github.com/SysBioChalmers/GECKO/blob/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction).

-	`geckomat/change_model/manualModifications.m`: Function that performs several modifications on the created ecModel, namely, addition of manually curated kcat values, removal of duplicated reactions and unused enzymes, etc. For the ecRhto-GEM case, this function has been simplified by removing all steps regarding the addition of curated *kcat* values (removal of lines **9-117** and **186-462**). 

-	`geckomat/limit_proteins/sumBioMass.m`: Function that computes the mass fraction for the most abundant macromolecular groups in the biomass composition of the model (lipids, carbohydrates, protein, nucleic acids). As each GEM has a particular biomass pseudoreaction structure, it is not possible to generalize this process. A step-by-step guide for update of this function in GECKO to the rhto-GEM case is provided as a MATLAB live script [rhtoGEM_adapt_sumBiomass.mlx](https://github.com/SysBioChalmers/GECKO/blob/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction).

-	`geckomat/limit_proteins/scaleBioMass.m`: Function that rescales the stoichiometric coefficients of the carbohydrates and lipids precursors in the biomass biosynthesis pseudoreactions according to the provided total protein content of the modeled strain. The rescaling procedure is based on the two following assumptions.

1. The total summed mass of protein, carbohydrates and lipids is constant.
2. The carbohydrates to lipids mass ratio is constant.

For the particular case of rhto-GEM, this function does not need to be updated, as the biomass pseudoreaction, and those for synthesis of its macromolecular components, conserve the same nomenclature and structure as those in the original script for enhancement of yeast-GEM.

-	`geckomat/kcat_sensitivity_analysis/changeMedia_batch.m`: Function that sets constraints for batch growth (unlimited availability of resources) according to a specific medium formulation suitable for the modelled organism. This function sets arbitrarily high flux values (e.g., 1000 mmol/ gDw h) as upper bounds for the exchange reactions, in the uptake direction (¬{ø} -> compound[extracellular]), for all compounds present in the desired culture medium, while blocking all other uptake reactions. All ecModels require uptake of the same compounds as their GEM counterpart for simulation of cellular growth, as the metabolites present in the biomass pseudoreaction are not changed by the GECKO pipeline.

In order to adapt this function to any other model, uptake of the main carbon source, oxygen and essential micronutrients should be allowed; production of biomass must be enabled, and production of the carbon source and oxygen must be blocked, for avoiding non-meaningful solutions. 

For the rhto-GEM case, information about its minimal nutrient requirements for growth was not found in its original publication, neither in its respective code repository. Therefore, an iterative approach for finding the minimal set of active uptake reactions is described in the live script [rhtoGEM.mlx](https://github.com/SysBioChalmers/GECKO/tree/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction). This script identified iron (2+) and sulphate uptake reactions, `r_1861` and `r_2060`, respectively, as essential for growth in rhto-GEM. Additionally, some other exchange reactions relevant for growth are identified by this script. In order to adapt the changeMedia_batch.m function to the rhto-GEM case, substitute lines 30 to 53 by the following:

	model = setParam(model, 'ub', 'r_1861_REV', +1000); % 'iron(2+) exchange';
	model = setParam(model, 'ub', 'r_2060_REV', +1000); % 'sulphate exchange';
	model = setParam(model, 'ub', 'r_1992_REV', +1000);  % oxygen uptake;
	model = setParam(model, 'ub', 'r_2005_REV', +1000);  % phosphate uptake;
	model = setParam(model, 'ub', 'r_1654_REV', +1000);  % ammonium uptake;
	model = setParam(model, 'ub', 'r_1714', 0);                    % block D-glucose production';
	model = setParam(model, 'ub', 'r_1992', 0);                    % block oxygen production';
	model = setParam(model, 'ub', 'r_2111', +1000);            % allow biomass production
	
-	`geckomat/limit_proteins/fitGAM.m`: Function that fits the stoichiometric coefficient for ATP in the biomass pseudoreaction (growth-associated maintenance, GAM), after rescaling of biomass macromolecular precursors by the scaleBioMass.m function. This function iterates, spanning GAM values from 20 – 100 mmol ATP per unit of biomass created. In each of these iterations, chemostat conditions are simulated by the simulateChemostat nested function, using a different GAM value at a time in the biomass pseudoreaction, by fixing the lower bound on the growth pseudoreaction equal to the desired dilution rate (indicated by the data), and solving an FBA problem with minimization of glucose uptake rate as an objective function. From the solution vector, values corresponding to exchange fluxes for glucose, O2 and CO2 are obtained, then the rooted-squared mean error (RMSE) is calculated between the data and the model predictions. After running a finite number of iterations, the GAM value which induced the lowest RMSE is chosen and set in the biomass pseudoreaction. 

Typically, if the experimental data in the chemostatData.tsv file complies with the format specified above, then this function does not need any further customization. However, as just glucose uptake rates at different dilution rates were found for the *R. toruloides* case, a modified chemostatData.tsv, with less columns than expected, is provided. Therefore, the fitGAM.m function needs to be adapted in order to be able to deal with such a modified data input (line 26 in the original script). For this case, the introduced modifications enable the function to read a two-columns data file instead of a four-columns one. Additionally, the lines corresponding to the calculation of the RMSE between the solution values and the experimental data were modified in order to base this calculation just on the glucose uptake rate (lines 106-109). Finally, the lines encoding for the plotting of the experimental data and model simulation outputs also need to be adjusted to the specific number of exchange fluxes available in this example (line 64). A modified version of the fitGAM function for the rhto-GEM case can be found in the folder: [rhtoGEM_scripts](https://github.com/SysBioChalmers/GECKO/tree/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction/rhtoGEM_scripts).

## 3. Generate ecRhtoGEM

A detailed step-by-step guide for enhancing rhto-GEM with enzyme constraints is presented in the live script [ecRhtoGEM.mlx](https://github.com/SysBioChalmers/GECKO/tree/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction)


## 4. Manual curation of kinetic parameters in ecRhtoGEM

As GECKO may produce overconstrained ecModels (not capable of reaching the provided experimental maximum growth rate in FBA simulations), specially for non-model organisms for which kinetic and protein data are not widely avaialable, manual curation of the kinetic parameters for the top-limiting enzymes in the model is required.

A detailed guide for curation of the limiting kinetic constraints in ecRhtoGEM, in order to generate a functional ecModel, is provided in [ecRhtoGEM_curation.mlx](https://github.com/SysBioChalmers/GECKO/tree/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction).

## 5. Validation of ecRhtoGEM predictions

A hands-on example on how to validate phenotype predictions of an ecModel is provided in [ecRhtoGEM_validation.mlx](https://github.com/SysBioChalmers/GECKO/tree/feat/add_tutorial/tutorials/ecRhtoGEM_reconstruction). In this live script prediction of maximum growth rates, subjected to constrained uptake rates of different carbon sources is presented and compared to experimental measurements.
