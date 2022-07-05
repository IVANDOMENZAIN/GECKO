% ecRhtoGEMUpdate
%
%   Ivan Domenzain, 2022-07-05
%

%Clone the necessary repos:
git clone https://github.com/SysBioChalmers/GECKO.git
%Load rhto model:
git clone https://github.com/SysBioChalmers/rhto-GEM.git
model    = load('rhto-GEM/ModelFiles/mat/rhto.mat');
model    = model.model;
modelVer = model.description(strfind(model.description,'_v')+1:end);
%Replace scripts in GECKO:
fileNames = dir('rhto_scripts');
for i = 1:length(fileNames)
    fileName = fileNames(i).name;
    if ~ismember(fileName,{'.' '..' '.DS_Store'})
        fullName   = ['rhto_scripts/' fileName];
        GECKO_path = dir(['GECKO/**/' fileName]);
        GECKO_path = GECKO_path.folder;
        copyfile(fullName,GECKO_path)
    end
end
%Run GECKO pipeline:
cd GECKO
delete databases/prot_abundance.txt
GECKOver = git('describe --tags');
cd geckomat/get_enzyme_data
updateDatabases;
cd ..
[ecModel,ecModel_batch] = enhanceGEM(model,'COBRA','ecRhtoGEM',modelVer);
cd ../..
%Move model files:
%rmdir('model', 's')
movefile GECKO/models/ecRhtoGEM models
save('models/ecRhtoGEM.mat','ecModel')
save('models/ecRhtoGEM_batch.mat','ecModel_batch')
%Save associated versions:
fid = fopen('dependencies.txt','wt');
fprintf(fid,['GECKO\t' GECKOver '\n']);
fprintf(fid,['rhto-GEM\t' modelVer '\n']);
fclose(fid);
%Remove the cloned repos:
rmdir('GECKO', 's')
rmdir('rhto-GEM', 's')
