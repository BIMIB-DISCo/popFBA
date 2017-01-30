function [popModel, singleModel, optFlux] = popFBA(nameSBML, nameExRxns, nameCoopRxn, nPop, CharExtComp, otherFeat, rxnsFeat, metsFeat)
%[popModel, singleModel, optFlux] = popFBA(nameSBML, nameExRxns, nameCoopRxn, nPop, CharExtComp, otherFeat, rxnsFeat, metsFeat)
%
%
%nameExRxns and nameCoopRxn are cell arrays with names of reactions that will be
%considered as exchange reactions and uptake reactions respectively.
%
%CharExtComp =  char to define the tumor micro envirorment. Default: 's'.
%
%numPop = number of subpopulations in the popModel
%
%OPTIONAL:
% In order to specify specific names fot the model structure fields to be incorporated in the new model the
% following arguments can be used:
%
% otherFeat = name of any field with no correspondence with rxns or mets e.g. list of genes
% rxnsFeat = name of any field associated with model.rxns
% metsFeat = name of any field associated with model.mets
%
% if not specified the function will try to deduce from the size of model.rxns and
% model.mets. If these two value are the same all arguments must be specified

if nargin < 4
    disp('Not enough input arguments');
    return
end
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
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, nPop, CharExtComp);
elseif nargin < 6
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, nPop, CharExtComp);
elseif nargin < 7
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, nPop, CharExtComp, otherFeat);
else
    popModel = createPopModel(singleModel, idxExRxns, idxUptkRxn, nPop, CharExtComp, otherFeat, rxnsFeat, metsFeat);
end

optFlux = optimizeCbModel(popModel);

end