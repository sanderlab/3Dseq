% Process Pacbio Sequel data from PSE-1 evolution experiments
% 
% INPUT
% 1. The fastq file containing Pacbio Sequel reads for round 20 of PSE-1 
%    evolution (P20.fastq) 
%
% OUTPUT
% Outputs processed sequences in fasta format:
% 1. Between forward and reverse anchor sequences that correspond to
%    primer sequences used for PSE-1 epPCR
% 2. All nucleotides with Q-score > 30
% 3. Full-length and truncated to 266 amino acids of pdb 1g68 sequence
% 4. Only unique amino acid sequences
% 5. Remove possible contaminating sequences having >=6 adjacent amino acid
%    mutations
%
% Copyright (c) 2019 Michael A. Stiffler    mstiffler@post.harvard.edu

%% Preliminary settings

Qcut=30; % Qscore cutoff for all nucleotides in gene

% forward read anchor sequences
fwd_anchor1='GGTTTTTGC';% 9 bp of 3' end of forward mutagenesis primer
fwd_anchor2='TAACTGTCAG';% 7 bp of 3' end of reverse mutagenesis primer + 3 bp of one of the stop codons
fwd_anchor3='TAGCTGTCAG';% 7 bp of 3' end of reverse mutagenesis primer + 3 bp of one of the stop codons
fwd_anchor4='TGACTGTCAG';% 7 bp of 3' end of reverse mutagenesis primer + 3 bp of one of the stop codons
fwd_anchor1_len=size(fwd_anchor1,2);
fwd_anchor2_len=size(fwd_anchor2,2);

% reverse complement read anchor sequences
rev_anchor1='CTGACAGTTA'; % 7 bp of 3' end of reverse mutagenesis primer + 3 bp of one of the stop codons
rev_anchor2='CTGACAGCTA'; % 7 bp of 3' end of reverse mutagenesis primer + 3 bp of one of the stop codons
rev_anchor3='CTGACAGTCA'; % 7 bp of 3' end of reverse mutagenesis primer + 3 bp of one of the stop codons
rev_anchor4='GCAAAAACC'; % 9 bp of 3' end of forward mutagenesis primer 
rev_anchor1_len=size(rev_anchor1,2);
rev_anchor4_len=size(rev_anchor4,2);

%% load fastq file - reads and quality scores
[~,reads,qual]=fastqread('P20.fastq'); % process round 20
num_reads_init=size(reads,2); % note initial number of reads

%% Process reads - anchors and quality score
% process for forward anchor sequences
parfor i=1:size(reads,2)
    f1=strfind(reads{i},fwd_anchor1);
    r1=strfind(reads{i},fwd_anchor2);
    r2=strfind(reads{i},fwd_anchor3);
    r3=strfind(reads{i},fwd_anchor4);
    if f1 % check if the 5' anchor is there
        f1=f1(1); 
        if r1 % check if the 3' anchor is there
            r1=r1(1);
            seqsf{i}=(reads{i} (f1+fwd_anchor1_len+1:r1-1));
            seqsf_qual{i}=(qual{i} (f1+fwd_anchor1_len+1:r1-1));
        end
        if r2 % check if the 3' anchor is there
            r2=r2(1);
            seqsf{i}=(reads{i} (f1+fwd_anchor1_len+1:r2-1));
            seqsf_qual{i}=(qual{i} (f1+fwd_anchor1_len+1:r2-1));
        end
        if r3 % check if the 3' anchor is there
            r3=r3(1);
            seqsf{i}=(reads{i} (f1+fwd_anchor1_len+1:r3-1));
            seqsf_qual{i}=(qual{i} (f1+fwd_anchor1_len+1:r3-1));
        end
    end
end
parfor i=1:size(seqsf,2)
    seqsf(i)=cellstr(char(seqsf{i}));
    seqsf_qual(i)=cellstr(char(seqsf_qual{i}));
end
seqsf_trunc=seqsf(~strcmp(seqsf,{''}));
seqsf_qual_trunc=seqsf_qual(~strcmp(seqsf,{''}));
num_reads_pass_fwdanchor=size(seqsf_trunc,2); % note number of reads with forward anchor sequences
clear seqsf seqsf_qual
seqsf_Q_cut=false(size(seqsf_trunc,2),1);
parfor i=1:size(seqsf_trunc,2)
    seqsf_Q_cut(i,1)=(min(double(char(seqsf_qual_trunc{i})))-33)>=Qcut; % make Q score cutoff
end
seqs_fwd=seqsf_trunc(seqsf_Q_cut);
clear seqsf_qual_trunc seqsf_Q_cut seqsf_trunc

% process for reverse complement read anchors
parfor i=1:size(reads,2)
    f1=strfind(reads{i},rev_anchor1);
    f2=strfind(reads{i},rev_anchor2);
    f3=strfind(reads{i},rev_anchor3);
    r1=strfind(reads{i},rev_anchor4);
    if f1 % check if the 3' anchor is there
        f1=f1(1); 
        if r1 % check if the 5' anchor is there
            r1=r1(1);
            seqsr{i}=(reads{i} (f1+rev_anchor1_len:r1-2));
            seqsr_qual{i}=(qual{i} (f1+rev_anchor1_len:r1-2));
        end
    end
    if f2 % check if the 3' anchor is there
        f2=f2(1); 
        if r1 % check if the 5' anchor is there
            r1=r1(1);
            seqsr{i}=(reads{i} (f2+rev_anchor1_len:r1-2));
            seqsr_qual{i}=(qual{i} (f2+rev_anchor1_len:r1-2));
        end
    end
    if f3 % check if the 3' anchor is there
        f3=f3(1); 
        if r1 % check if the 5' anchor is there
            r1=r1(1);
            seqsr{i}=(reads{i} (f3+rev_anchor1_len:r1-2));
            seqsr_qual{i}=(qual{i} (f3+rev_anchor1_len:r1-2));
        end
    end
end
parfor i=1:size(seqsr,2)
    seqsr(i)=cellstr(char(seqsr{i}));
    seqsr_qual(i)=cellstr(char(seqsr_qual{i}));
end
clear reads qual
seqsr_trunc=seqsr(~strcmp(seqsr,{''}));
seqsr_qual_trunc=seqsr_qual(~strcmp(seqsr,{''}));
num_reads_pass_revanchor=size(seqsr_trunc,2); % note number of reads with reverse complement anchor sequences
clear seqsr seqsr_qual
seqsr_Q_cut=false(size(seqsr_trunc,2),1);
parfor i=1:size(seqsr_trunc,2)
    seqsr_Q_cut(i,1)=(min(double(char(seqsr_qual_trunc{i})))-33)>=Qcut; % make Q score cutoff
end
seqs_rev=seqsr_trunc(seqsr_Q_cut);
clear seqsr_qual_trunc seqsr_Q_cut seqsr_trunc
parfor i=1:size(seqs_rev,2)
    seqs_rev{i}=seqrcomplement(char(seqs_rev(i)));
end

% combine sequences (after this point refer to as seqs instead of
% reads)
num_reads_pass_fwdanchor_Qcut=size(seqs_fwd,2);
num_reads_pass_revanchor_Qcut=size(seqs_rev,2);
seqs=seqs_fwd;
seqs(size(seqs_fwd,2)+1:size(seqs_fwd,2)+size(seqs_rev,2))=seqs_rev;
clear seqs_fwd seqs_rev
num_seqs_init=size(seqs,2); % note number of sequences passing anchors and Q score cutoff

%% Process sequences for full-length and no internal stop codons

% load wt amino acid sequence and calculate length
pse_wt_aa='SSSKFQQVEQDVKAIEVSLSARIGVSVLDTQNGEYWDYNGNQRFPLTSTFKTIACAKLLYDAEQGKVNPNSTVEIKKADLVTYSPVIEKQVGQAITLDDACFATMTTSDNTAANIILSAVGGPKGVTDFLRQIGDKETRLDRIEPDLNEGKLGDLRDTTTPKAIASTLNKFLFGSALSEMNQKKLESWMVNNQVTGNLLRSVLPAGWNIADRSGAGGFGARSITAVVWSEHQAPIIVSIYLAQTQASMAERNDAIVKIGHSIFDVYTSQSR';
pse_wt_aa_len=size(pse_wt_aa,2);

% eliminate sequences not full-length (271 amino acids)
parfor i=1:size(seqs,2)
    seqs_len(i)=size(char(seqs{i}),2); % compute sequence lengths
end
seqs=seqs(seqs_len==(3*pse_wt_aa_len)); % take only sequences with same length as wt 
clear seqs_len
num_seqs_pass_len=size(seqs,2); % note number of full-length sequences

% convert to amino acid sequences
seqs_aa=nt2aa(seqs,'AlternativeStartCodons',false); % takes a bit

% eliminate sequences with internal stop codons
parfor i=1:size(seqs_aa,2)
    s1=strfind(seqs_aa{i},'*');
    if s1 % check if stop is present
        seqs_aa{i}=({''});
    end
end
parfor i=1:size(seqs_aa,2)
    seqs_aa(i)=cellstr(char(seqs_aa{i}));
end
seqs_aa=seqs_aa(~strcmp(seqs_aa,{''}));
num_seqs_pass_stopcodon=size(seqs_aa,2); % note number of sequences passing without internal stop codons

% truncate to length of pdb 1g68 sequence
seqs_aa=char(seqs_aa);
seqs_aa=seqs_aa(:,3:268); % eliminate 2 N-terminal and 3 C-terminal amino acids
seqs_aa=cellstr(seqs_aa);
% same for wt pse sequence
pse_wt_aa_trunc=pse_wt_aa(1,3:268);

% unique amino acid sequences
seqs_aa_uniq=unique(seqs_aa);

%% Clean for possible contaminating sequences (sequences with >=6 adjacent amino acid mutations)

% Look at mutations per sequence and position - mutations not matching wt
% PSE1 truncated to 1g68 sequence
seqs_aa_uniq=char(seqs_aa_uniq);
parfor a=1:size(seqs_aa_uniq,2)
    seq_id_psewt(:,a)=ismember(seqs_aa_uniq(1:size(seqs_aa_uniq,1),a),pse_wt_aa_trunc(1,a)); % mutations = false (0)
end
% sum number of mutations in 6 amino acid windows
parfor i=1:(size(pse_wt_aa_trunc,2)-5)
    seq_nummut_6scan(:,i)=6-sum(seq_id_psewt(:,i:i+5),2);
end
seq_nummut_6scan_max=max(seq_nummut_6scan');

% eliminate contaminant sequences having 6 mutations in any window
seqs_aa_uniq_fin=seqs_aa_uniq(find(seq_nummut_6scan_max<6),:);

clear seq_id_psewt seq_nummut_6scan seq_nummut_6scan_max

%% save data
save Round20_ProcessPacBioSequel.mat

%% save final alignment
fastawrite('PSE-1_Round20.fas',num2str((1:size(seqs_aa_uniq_fin,1))'),seqs_aa_uniq_fin);