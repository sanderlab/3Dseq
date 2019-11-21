function [results] = n_choose_2_rmsd(tblFile,n_resis)

% Reads in .csv file with column 'filename' that points to model files
list_pdbs = readtable(tblFile,'Delimiter',',');

expPDBfile = list_pdbs.filename{1};
expPDB = pdbread(expPDBfile);
if nargin < 2
    n_resis = length(unique([expPDB.Model.Atom(:).resSeq]));
end

% Initialize variables
N = height(list_pdbs);
l = nchoosek(N,2);

results.dist_matrix = zeros(l,n_resis);
results.model_i = cell(l,1);
results.model_j = cell(l,1);

% Pre-read pdbs into memory once to make faster
pdbs = cell(N,1);
for i=1:N
    pdbs{i} = pdbread(list_pdbs.filename{i});
    pdbs{i}.Sequence.ChainID = pdbs{i}.Model.Atom(1).chainID;
end

% Count how many of the N choose 2 we've gone through so far
x = 1;

% Iterate over all models
for i=1:N
    disp(i);
    
    % Set current model as the base to compare all others against
    expPDBfile = list_pdbs.filename{i};
    expPDB = pdbs{i};
    
    % Get locations of C-alpha atoms of base model
    chain_of_interest = expPDB.Model.Atom(1).chainID;
    z_i = find(ismember({expPDB.Model.Atom.AtomName},'CA') & ismember({expPDB.Model.Atom.chainID},chain_of_interest) & (cellfun('isempty',{expPDB.Model.Atom.altLoc}) | ismember({expPDB.Model.Atom.altLoc},'A')));
    
    % Need to map into actual numbering
    res_i = [expPDB.Model.Atom(z_i).resSeq];
    
    X_i = nan(1,n_resis);
    X_i(res_i) = [expPDB.Model.Atom(z_i).X];
    Y_i = nan(1,n_resis);
    Y_i(res_i) = [expPDB.Model.Atom(z_i).Y];
    Z_i = nan(1,n_resis);
    Z_i(res_i) = [expPDB.Model.Atom(z_i).Z];
    
    % Go through all possible pairs of models that include the base model
    % that we haven't already calculated
    for j=i+1:N
        
        
        tempPDBfile = list_pdbs.filename{j};
        
        results.model_i{x} = expPDBfile;
        results.model_j{x} = tempPDBfile;
        
        tempPDB = pdbs{j};
        
        % Superpose the two models and return the non-base coordinates that
        % have been structurally aligned to the base
        [~,~,~,pdbm] = pdbsuperpose(expPDB,tempPDB,'Display',false);
        
        % Find locations of alpha carbons for the second model
        chain_of_interest = pdbm.Model.Atom(1).chainID;
        z_j = find(ismember({pdbm.Model.Atom.AtomName},'CA') & ismember({pdbm.Model.Atom.chainID},chain_of_interest) & (cellfun('isempty',{pdbm.Model.Atom.altLoc}) | ismember({pdbm.Model.Atom.altLoc},'A')));
        
        res_j = [pdbm.Model.Atom(z_j).resSeq];
        
        X_j = nan(1,n_resis);
        X_j(res_j) = [pdbm.Model.Atom(z_j).X];
        Y_j = nan(1,n_resis);
        Y_j(res_j) = [pdbm.Model.Atom(z_j).Y];
        Z_j = nan(1,n_resis);
        Z_j(res_j) = [pdbm.Model.Atom(z_j).Z];
        
        results.dist_matrix(x,:) = sqrt((X_i - X_j).^2 + (Y_i - Y_j).^2 + (Z_i - Z_j).^2);
        
        x = x+1;
    end
    
end

end