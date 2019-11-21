function [data] = pdbdists3(pdbname,chain,varargin)
%PDBDISTS calculates the distances between all pairs of residues in a protein.
%   [hncodistances,sc_residues,hncoprobes] = pdbdists(pdbname,chain,...) 
%
%INPUTS
%   pdbname - 4 character PDB name
%   chain - 'A', etc.
%   OPTIONAL
%       ires - restrict to given residues in 1st or "i'th" dimension
%       jres - restrict to given residues in 2nd or "j'th" dimension
%       eclude_prolines - exclude prolines in 1st or "i'th" dimension - useful to compare with NMR data
%
%COMMENTS: Multiple types of residue to residue distances are calculated
%
%DEPENDENCIES: getpdb.m (Bioinformatics toolbox)
%
%Alan Poole, 05/06/2009
%Bug fixes and adaptations by Frank J. Poelwijk, 02/13/2018

%% Initialization
% gets pdb from NCBI, needs internet
if ischar(pdbname)
    pdb = getpdb(pdbname);
else
    pdb = pdbname;
end
% get PDBID and sequence
try data.pdbid = pdb.Header.idCode; end

backboneatoms = char({'C';'N';'O';'OXT'});  % not side-chain atoms

% PDB numbering for each residue in the model
data.res = unique([pdb.Model.Atom(upper([pdb.Model.Atom.chainID]) == upper(chain)).resSeq]);
data.ires = data.res;
data.jres = data.res;
ixires = 1:numel(data.res);
ixjres = 1:numel(data.res);
exclude_prolines = false;

%% Optional Parameters
nvarargin = numel(varargin);
if nvarargin > 0
    if rem(nvarargin,2) == 1
        error('Bioinfo:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'ires' 'jres' 'exclude_prolines'};
    for j=1:2:nvarargin
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs); 
        if isempty(k)
            error('Bioinfo:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % ires (res)
                    [data.ires,ixires] = intersect(data.res,pval); 
                case 2  % jres (muts)
                    [data.jres,ixjres] = intersect(data.res,pval);
                case 3  % exclude_prolines
                    exclude_prolines = pval;
            end
        end
    end
end

%% Loop over all residues to calculate the center of side-chains % FJP: better would be to only calculate for ires and jres (instead of now cutting to those only at the end of the calculation).
for ii = 1:numel(data.res)
    indices = find([pdb.Model.Atom.resSeq] == data.res(ii) & upper([pdb.Model.Atom.chainID]) == upper(chain));  % Find atoms corresponding to each residue %FJP: only for the specoified chain!
    data.sequence(ii) = aminolookup(pdb.Model.Atom(indices(1)).resName);
    count = 1;
    clear sidechainatom*
    data.atomnames{ii} = {pdb.Model.Atom(indices).AtomName};    % get atom names
    data.atomnums{ii} = indices;
    for jj = 1:numel(indices)
        if ~numel(strmatch(data.atomnames{ii}{jj},backboneatoms,'exact'))    % for all sidechain atoms (not backbone)
            if ~strcmpi(pdb.Model.Atom(indices(jj)).resName,'GLY') && strcmp(data.atomnames{ii}{jj},'CA')
                % skip non-glycine C-alpha atoms
            else
                scatoms{count} = data.atomnames{ii}{jj};
                scatomnums(count) = indices(jj);
                count = count + 1;
            end
        end
    end
    data.atom_coords{ii} = [[pdb.Model.Atom(indices).X]' [pdb.Model.Atom(indices).Y]' [pdb.Model.Atom(indices).Z]'];
    data.res_center(ii,:) = mean(data.atom_coords{ii},1);
    data.sidechainatomnames{ii} = scatoms;
    data.sidechainatomnums{ii} = scatomnums;
    data.sc_coords{ii}=[[pdb.Model.Atom(data.sidechainatomnums{ii}).X]' [pdb.Model.Atom(data.sidechainatomnums{ii}).Y]' [pdb.Model.Atom(data.sidechainatomnums{ii}).Z]'];
    data.sc_center(ii,:) = mean(data.sc_coords{ii},1);
    clear scatoms scatomnums
end
clear count ii jj indices
data.sc2sc_dists = squareform(pdist(data.res_center));

%% Proline handling - this is for comparing with NMR data
prolines = strfind(upper(data.sequence(ixires)),'P');
if exclude_prolines
    [ixires] = setdiff(ixires,prolines);
    data.ires = data.res(ixires);
end

%% Loop over all residues to find the center of H-N-C(O) clusters
% This is specifically for comparing with NMR data, so prolines (no amide proton) and the 1st residue (no i-1 CO) are
% excluded.
clear hnco_centers nitrogen carbonyl
% since HNCO clusters the res (i) N and res (i-1) CO, we can't calculate the cluster for residue #1
ixhncores = setdiff(ixires,[1 prolines]);
data.hncores = data.res(ixhncores);

for ii = 1:numel(ixhncores)
%    nitrogen = strmatch('N',data.atomnames{ixhncores(ii)},'exact');     % find N atom if it exists, this is probably unnecessary - all residues have N
%    carbonyl = strmatch('C',data.atomnames{ixhncores(ii)-1},'exact');
    nitrogen = find(strcmp('N',data.atomnames{ixhncores(ii)}),1,'first');     % find N atom if it exists, this is probably unnecessary - all residues have N
    carbonyl = find(strcmp('C',data.atomnames{ixhncores(ii)-1}),1,'first');   % FJP: for cases where there are alternative locations (altLoc) in the pdb, pick first
    
    if exist('nitrogen','var') && exist('carbonyl','var') % && ~strcmpi(data.sequence(ii),'P');
        data.hnco_coordinates(ii,:,:) = [data.atom_coords{ixhncores(ii)}(nitrogen,:); data.atom_coords{ixhncores(ii)-1}(carbonyl,:)]; 
        
        data.hnco_center(ii,:) = mean(squeeze(data.hnco_coordinates(ii,:,:)));
    else
        data.hnco_center(ii,:) = NaN;
    end
    clear nitrogen carbonyl
end
clear ii

% Calculate distances
data.sc2hnco_dists = pdist2(data.hnco_center,data.sc_center);

%% Calculate contacts
data.allatoms = [[pdb.Model.Atom.X]' [pdb.Model.Atom.Y]' [pdb.Model.Atom.Z]'];
data.atomdists = squareform(pdist(data.allatoms));

elem = [pdb.Model.Atom.element];
atom_type=['C', 'O', 'N', 'S', 'H'];
atom_radius=[1.9, 1.4, 1.5, 1.85, 0];
for ii = 1:numel(atom_type)
    elemr(strmatch(atom_type(ii),elem')) = atom_radius(ii);
end

elemrm = ones(numel(elemr),1) * elemr;
contactatomdists = 1.2*(elemrm + elemrm');
contactatoms = data.atomdists <= contactatomdists;

junk = zeros(size(contactatoms,1),numel(data.res));
for ii = 1:numel(data.res)
    junk(:,ii) = sum(contactatoms(:,data.atomnums{ii}),2);
end
for ii = 1:numel(data.res)
    data.contacts(ii,:) = sum(junk(data.atomnums{ii},:)) > 0;
end

%% Restrict to given residues
data.sc2sc_dists = data.sc2sc_dists(ixires,ixjres);
data.sc2hnco_dists = data.sc2hnco_dists(:,ixjres);
data.contacts = data.contacts(ixires,ixjres);

