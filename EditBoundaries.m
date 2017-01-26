function [modelOut, TableRes] = EditBoundaries(model,rxn, newLb, newUb, NotSure)
%User friendly function to manage the boundaries of a FBA model or a
%population model

%Find all the reactions whose name contain the substring rxn in model.rxns
 
%With first and second arguments only show the reaction boundaries

%With also third and fourth arguments the function show the reaction
%boundaries and ask user confirmation before change them.

%With the fifth argument (NotSure) = false change the boundaries without ask 
%user confirmation, for use the function inside a script
%(Default true)


idxS =  strfind(lower(model.rxns),lower(rxn));
idx = find(not(cellfun('isempty', idxS)));
TableRes = table();

if nargin < 3
    if(isempty(idx))
        disp(rxn);
        disp('No reaction found');
        modelOut = model;
    else
        TableRes = table(idx, model.rxns(idx), model.lb(idx),  model.ub(idx), 'VariableNames', {'ID' 'Reaction', 'LowerBound', 'UpperBound'})
        disp('No change was made');
        modelOut = model;
        return
    end
elseif nargin < 5
    NotSure = true;
end
if(isempty(idx))
    disp(rxn);
    disp('No reaction found');
    modelOut = model;
elseif(NotSure)
    TableRes = table(idx, model.rxns(idx), model.lb(idx),  model.ub(idx), 'VariableNames', {'ID' 'Reaction', 'LowerBound', 'UpperBound'})
    choice = questdlg('Would you like to edit these boundaries?', 'Edit Boundaries', 'Yes','No','No');   
    switch choice
        case 'Yes'
            model.lb(idx) = newLb;
            model.ub(idx) = newUb;
            modelOut = model;
        case 'No'
            modelOut = model;
    end
else
    TableRes = table(idx, model.rxns(idx), model.lb(idx),  model.ub(idx), 'VariableNames', {'ID' 'Reaction', 'LowerBound', 'UpperBound'});
    model.lb(idx) = newLb;
    model.ub(idx) = newUb;
    modelOut = model;
end
end
