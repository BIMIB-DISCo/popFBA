function [popModel, singleModel, optFlux] = popFBA(nameSBML, nameExRxns, nameCoopRxn, CharExtComp, nPop, otherFeat, rxnsFeat, metsFeat)

%nameExRxns and nameUpRxn are cell array with name of reaction will be
%consider exchange reaction and uptake reaction rispectivelly.

%CharExtComp =  char to define the tumor micro envirorment. default 's'

%numPop = number of subpopulation in the popModel

% for specify the struct field to be incorporated in the new model pass the
% following arguments

% otherFeat = name of field with no link to rxns or mets e.g. list of genes
% rxnsFeat = name of field with rxns feature
% metsFeat = name of field with mets feature
% if not pass the function will try to determine from the size of rxns and
% mets. if this two value are the same must pass all the argument (or add a fake metabolite)


singleModel = readCbModel(nameSBML); %load from sbml file the single model

idxExRxns = [];
idxUptkRxn = [];

for i=1:length(nameExRxns)
    idxtmp = (find(strcmp(singleModel.rxns, nameExRxns(i))==1));
    if isempty(idxtmp)
        disp(nameExRxns(i)); disp('Not found');
    else
        idxExRxns = [idxExRxns idxtmp];
    end
end
for i=1:length(nameCoopRxn)
    idxtmp = (find(strcmp(singleModel.rxns, nameCoopRxn(i))==1));
    if isempty(idxtmp)
        disp(nameCoopRxn(i)); disp('Not found');
    else
        idxUptkRxn = [idxUptkRxn idxtmp];
    end
end
if nargin < 5
    CharExtComp = 's';
elseif nargin < 6
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, CharExtComp, nPop);
elseif nargin < 7
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, CharExtComp, nPop, otherFeat);
else
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, CharExtComp, nPop, otherFeat, rxnsFeat, metsFeat);
end

optFlux = optimizeCbModel(popModel);

end