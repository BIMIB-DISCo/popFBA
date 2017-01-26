function [popModel] = createPopModel(model, idxExRxns, idxCoopkRxn, CharExtComp, nPop, otherFeat, rxnsFeat, metsFeat)
%work with model loaded with CobraToolbox
%model must have at least the fields 'rxns', 'mets', 'S', 'c', 'lb', 'up'
%the metabolites must are in the form of mets[x] where x = char for the
%compartment

%le ExRnxs saranno trattate nella forma [x] <=> '' se x ~= CharExtComp la reazione
%viene moltiplicata per numPop. Si considera il fatto che l'utente può voler far
%entrare o uscire un componente direttamente dalla cellula senza farlo
%passare dalla matrice extracellulare (opzione sconsigliata)

%le UptkRxn saranno trattate nella forma [y] <=> [CharExtComp] (in <=> out)
%se la reazione di partenza è [y] <=> '' allora viene modificata in [y_n] <=> [CharExtComp]
%e verrà creata la corrispondente razione di exchange

%tutte le altre reazioni saranno trattate nella forma [?] <=> [?]
%trasformandole in [?_n] <=> [?_n]


%CharExtComp = char to discriminate between the extracellular and intercellular
%compartment e.g. 's'

%numPop is the number of time of model multiplication

%rxnsFeat, metsFeat, otherFeat are string array with the name of field
%linked rispectivelly to rxns, mets or niether of them. If theese argument
%are not passed the feature will automatically determinate

RxnsFeatures = table();
MetsFeatures = table();
OtherFeature = [];
NumRxns = size(model.rxns,1); %problema se ci sono stesso numero di metaboliti e reazioni
NumMets = size(model.mets,1); %prima soluzione approssimativa è far dare nome dei campi che si vogliono conservare direttamente dall'utente

if nargin < 5
    disp('Not enought argument');
    popModel = model;
    return
elseif nargin < 6
    if(NumRxns==NumMets)
        disp('Warning there is the same numbers of metabolite and reaction. It is not possible to distinguish between rnxs feature and mets feature.');
        disp('Please pass all the argument');
        return
    else
        %Parse all model field
        NameOfField = fieldnames(model);
        NameOfField(find(strcmp(NameOfField,'rxns')==1))=[]; %remove field rxns, mets, S from field to parse
        NameOfField(find(strcmp(NameOfField,'mets')==1))=[];
        NameOfField(find(strcmp(NameOfField,'S')==1))=[];
        
        for i=1:length(NameOfField)
            if(size(model.(NameOfField{i}),1)==NumRxns)
                tmpTable = table(model.(NameOfField{i}), 'VariableNames', NameOfField(i)); %feature table: name of columns is the name of model fields
                RxnsFeatures = [RxnsFeatures tmpTable]; %  is a table Rxns x (Feature belong to rxns)
            elseif (size(model.(NameOfField{i}),1)==NumMets)
                tmpTable = table(model.(NameOfField{i}), 'VariableNames', NameOfField(i)); %feature table: name of columns is the name of model fields
                MetsFeatures = [MetsFeatures tmpTable]; %  is a table Mets x (Feature belong to Mets)
            else
                
                OtherFeature = [OtherFeature NameOfField(i)]; % is a table within features not linked neither mets or rxns e.g. List of genes
                
            end
        end
    end
    
elseif nargin < 7 % only other features are explicitate so the Rxns and Mets feature will be determinate
    %Parse all model field
    NameOfField = fieldnames(model);
    NameOfField(find(strcmp(NameOfField,'rxns')==1))=[]; %remove field rxns, mets, S from field to parse
    NameOfField(find(strcmp(NameOfField,'mets')==1))=[];
    NameOfField(find(strcmp(NameOfField,'S')==1))=[];
    for i=1:length(otherFeat)
        NameOfField(find(strcmp(NameOfField,otherFeat(i))==1))=[]; %remove also the "other" feature passed by user       
    end
    OtherFeature = otherFeat;
    for i=1:length(NameOfField)
        tmpTable = table(model.(NameOfField{i}), 'VariableNames', NameOfField(i)); %feature table: name of columns is the name of model fields
        if(size(model.(NameOfField{i}),1)==NumRxns)
            RxnsFeatures = [RxnsFeatures tmpTable]; %  is a table Rxns x (Feature belong to rxns)
        else
            MetsFeatures = [MetsFeatures tmpTable]; %  is a table Mets x (Feature belong to Mets)
        end
    end
    
else %all argument was passed
    for i=1:length(rxnsFeat)
        tmpTable = table(model.(rxnsFeat{i}), 'VariableNames', rxnsFeat(i));
        RxnsFeatures = [RxnsFeatures tmpTable];
    end
    for i=1:length(metsFeat)
        tmpTable = table(model.(metsFeat{i}), 'VariableNames', metsFeat(i));
        MetsFeatures = [MetsFeatures tmpTable];
    end
    OtherFeature = otherFeat;
end

%Initialize var
S = []; %new matrix S: n x 3 | idxMets, idxRxn, value
nPop = nPop - 1;
MatrNewMets = table(); %tabella (numPop*Mets) x 3 + k con nuovo metabolita, nuovo indice e vecchio indice + k features
MatrNewRxns = table(); %tabella (numPop*Rxns) x 3 + p con nuovo metabolita, nuovo indice e vecchio indice + p features

%Exchage reaction add prefix "Ex_" for debug
for i=1:length(idxExRxns)
    
    
    idxMets = find(model.S(:,idxExRxns(i))<0); %deve essere solo uno a sinistra della reazione,
    %se è più di uno ipotizzo che sia lo stesso metabolita
    %con il compartimento diverso quindi prendo solo il primo dell'elenco, rimuovo
    %quello che c'è tra [] e aggiungo '[s]'
    tmpMet = model.mets{idxMets(1)};
    compChar = tmpMet(end-1);
    if(compChar == CharExtComp)
        tmpTable = table({strcat(model.rxns{idxExRxns(i)})}, idxExRxns(i), 'VariableNames', {'NewRxnName', 'RxnsOldIdx'});
        tmpTable = [tmpTable RxnsFeatures(idxExRxns(i),:)]; % add Rxn(i) feature
        MatrNewRxns = [MatrNewRxns; tmpTable];
        
        tmpTable = table({tmpMet}, idxMets(1), 'VariableNames', {'NewMetName', 'MetsOldIdx'});%rimossa parte [x]
        tmpTable = [tmpTable MetsFeatures(idxMets(1),:)];
        MatrNewMets = [MatrNewMets; tmpTable];
        S = [S; size(MatrNewMets,1) size(MatrNewRxns,1) model.S(idxMets(1),idxExRxns(i))]; %save only one rxn and one met for each i
    else
        for j=0:nPop
            tmpTable = table({strcat(model.rxns{idxExRxns(i)}, '_', num2str(j))}, idxExRxns(i), 'VariableNames', {'NewRxnName', 'RxnsOldIdx'});
            tmpTable = [tmpTable RxnsFeatures(idxExRxns(i),:)]; % add Rxn(i) feature
            MatrNewRxns = [MatrNewRxns; tmpTable];
            
            tmpTable = table({strcat(tmpMet(1:end-3), '_', num2str(j), '[',compChar,']')}, idxMets(1), 'VariableNames', {'NewMetName', 'MetsOldIdx'});%rimossa parte [x]
            tmpTable = [tmpTable MetsFeatures(idxMets(1),:)];
            MatrNewMets = [MatrNewMets; tmpTable];
            S = [S; size(MatrNewMets,1) size(MatrNewRxns,1) model.S(idxMets(1),idxExRxns(i))];
        end
    end
end

%Uptake reaction:"Uptk_" for debug
for i=1:length(idxCoopkRxn)
    idxMets = find(model.S(:,idxCoopkRxn(i)));
    if length(idxMets) == 1 % [y] <=> '' to [y_n] <=> [CharExtComp] and [CharExtComp] <=> ''
        idxRxnSt = size(MatrNewRxns,1)+1;
        for j=0:nPop %crea le n rxn
            tmpTable = table({strcat(model.rxns{idxCoopkRxn(i)}, '_', num2str(j))}, idxCoopkRxn(i), 'VariableNames', {'NewRxnName', 'RxnsOldIdx'}); %crea le n rxn
            tmpTable = [tmpTable RxnsFeatures(idxCoopkRxn(i),:)];
            MatrNewRxns = [MatrNewRxns; tmpTable];
        end
        tmpMet = model.mets{idxMets};
        compChar = tmpMet(end-1);
        for j=0:nPop %crea gli n mets[y]
            tmpTable = table({strcat(tmpMet(1:end-3), '_', num2str(j), '[',compChar,']')}, idxMets, 'VariableNames', {'NewMetName', 'MetsOldIdx'});
            tmpTable = [tmpTable MetsFeatures(idxMets,:)];
            MatrNewMets = [MatrNewMets; tmpTable];
            S = [S; size(MatrNewMets,1) (idxRxnSt + j)  model.S(idxMets,idxCoopkRxn(i))];
        end
        %creazione metabolita [CharExtComp] (nuovo metabolita)
        tmpTable = table({strcat(tmpMet(1:end-3), '[',CharExtComp,']')}, -1, 'VariableNames', {'NewMetName', 'MetsOldIdx'});
        tmpTable = [tmpTable MetsFeatures(idxMets,:)];
        if isfield(table2struct(tmpTable),'metNames') && isfield(table2struct(tmpTable),'metCompartment')
        tmpTable.('metNames') = {strcat(tmpMet(1:end-3),CharExtComp)};
        tmpTable.('metCompartment') = {CharExtComp};
        end
        MatrNewMets = [MatrNewMets; tmpTable];
        for j=0:nPop
            S = [S; size(MatrNewMets,1) (idxRxnSt + j)  -1*(model.S(idxMets,idxCoopkRxn(i)))]; %
        end
        tmpTable = table({strcat(model.rxns{idxCoopkRxn(i)}, '_new')}, -1, 'VariableNames', {'NewRxnName', 'RxnsOldIdx'});
        tmpTable = [tmpTable RxnsFeatures(idxCoopkRxn(i),:)];
        tmpTable.('NewRxnName') = {strcat(model.rxns{idxCoopkRxn(i)}, '_new')};
        MatrNewRxns = [MatrNewRxns; tmpTable];
        S = [S; size(MatrNewMets,1) size(MatrNewRxns,1)  -1];
    else        %[y] <=> [CharExtComp] to [y_n] <=> [CharExtComp]
        %Bisogna iniziare a controllare che il metabolità non sia già in elenco
        idxRxnSt = size(MatrNewRxns,1)+1;
        for j=0:nPop %crea le n rxn
            tmpTable = table({strcat(model.rxns{idxCoopkRxn(i)}, '_', num2str(j))}, idxCoopkRxn(i), 'VariableNames', {'NewRxnName', 'RxnsOldIdx'}); %crea le n rxn
            tmpTable = [tmpTable RxnsFeatures(idxCoopkRxn(i),:)];
            MatrNewRxns = [MatrNewRxns; tmpTable];
        end
        tmpMet = model.mets(idxMets); %nomi metaboliti coinvolti nella reazione
        for k=1:length(idxMets)
            if ~isempty(MatrNewMets)
                newIdx = find(ismember(MatrNewMets.('MetsOldIdx'), idxMets(k)));
            else
                newIdx = [];
            end
            if isempty(newIdx)
                if tmpMet{k}(end-1) == CharExtComp %nuovo metabolita in [CharExtComp]
                    tmpTable = table({tmpMet{k}}, idxMets(k), 'VariableNames', {'NewMetName', 'MetsOldIdx'});
                    tmpTable = [tmpTable MetsFeatures(idxMets(k),:)];
                    MatrNewMets = [MatrNewMets; tmpTable];
                    for j=0:nPop
                        S = [S; size(MatrNewMets,1) (idxRxnSt + j)  (model.S(idxMets(k),idxCoopkRxn(i)))]; %
                    end
                else
                    compChar = tmpMet{k}(end-1);
                    for j=0:nPop %crea gli n mets[y]
                        tmpTable = table({strcat(tmpMet{k}(1:end-3), '_', num2str(j), '[',compChar,']')}, idxMets(k), 'VariableNames', {'NewMetName', 'MetsOldIdx'});
                        tmpTable = [tmpTable MetsFeatures(idxMets(k),:)];
                        MatrNewMets = [MatrNewMets; tmpTable];
                        S = [S; size(MatrNewMets,1) (idxRxnSt + j)  model.S(idxMets(k),idxCoopkRxn(i))];
                    end
                end
            else %metabolita già in elenco
                if tmpMet{k}(end-1) == CharExtComp %nuovo metabolita in [CharExtComp]
                    for j=0:nPop
                        S = [S; newIdx(1) (idxRxnSt + j)  (model.S(idxMets(k),idxCoopkRxn(i)))]; %
                    end
                else
                    for j=0:nPop %crea gli n mets[y]
                        S = [S; (newIdx(1)+j) (idxRxnSt + j)  model.S(idxMets(k),idxCoopkRxn(i))];
                    end
                end
            end
        end
    end
end
%Moltiplicazione Inter reaction: from [?] <=> [?] to [?_n] <=> [?_n]
%saranno n = numPop e chiamate rxn_n
idxRxn = 1:NumRxns;
idxRxn = setdiff(idxRxn, idxExRxns);
idxRxn = setdiff(idxRxn, idxCoopkRxn); %rimosse reazioni già trattate
for i=1:length(idxRxn) % [?]<=>[?] to [?_n]<=>[?_n]
    idxRxnSt = size(MatrNewRxns,1)+1;
    for j=0:nPop %crea le n rxn
        tmpTable = table({strcat(model.rxns{idxRxn(i)}, '_', num2str(j))}, idxRxn(i), 'VariableNames', {'NewRxnName', 'RxnsOldIdx'}); %crea le n rxn
        tmpTable = [tmpTable RxnsFeatures(idxRxn(i),:)];
        MatrNewRxns = [MatrNewRxns; tmpTable];
    end
    idxMets = find(model.S(:,idxRxn(i)));
    tmpMet = model.mets(idxMets);
    for k=1:length(idxMets)
        compChar = tmpMet{k}(end-1);
        if ~isempty(MatrNewMets)
            newIdx = find(ismember(MatrNewMets.('MetsOldIdx'), idxMets(k)));
        else
            newIdx = [];
        end
        if isempty(newIdx)
            for j=0:nPop
                tmpTable = table({strcat(tmpMet{k}(1:end-3), '_', num2str(j), '[',compChar,']')}, idxMets(k), 'VariableNames', {'NewMetName', 'MetsOldIdx'});
                tmpTable = [tmpTable MetsFeatures(idxMets(k),:)];
                MatrNewMets = [MatrNewMets; tmpTable];
                S = [S; size(MatrNewMets,1) (idxRxnSt + j)  model.S(idxMets(k),idxRxn(i))];
            end
        else
            for j=0:nPop
                S = [S; (newIdx(1)+j) (idxRxnSt + j)  model.S(idxMets(k),idxRxn(i))];
            end
        end
    end
end
popModel.S = sparse(S(:,1),S(:,2),S(:,3));
popModel.rxns = MatrNewRxns.('NewRxnName');
popModel.mets = MatrNewMets.('NewMetName');
for i=2:size(MatrNewRxns,2)
    popModel.(MatrNewRxns.Properties.VariableNames{i}) = MatrNewRxns{:,i};
end
for i=2:size(MatrNewMets,2)
    popModel.(MatrNewMets.Properties.VariableNames{i}) = MatrNewMets{:,i};
end
for i=1:size(OtherFeature,2)
    popModel.(OtherFeature{i}) = model.(OtherFeature{i});
end
end