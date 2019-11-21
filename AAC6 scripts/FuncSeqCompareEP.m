%%% function assigning type of mutations based on sequence alignment. Mind
%%% that this procedure is not unique, but the script works well for typical
%%% mutation rates encountered in the mutagenesis study.
%%% Copyright (c) 2019 Frank J. Poelwijk, poelwijk@gmail.com.

function [MUTS,MUTpos,MUTcode,MUTcodeAA,MUTAAname,classMUT,score]=FuncSeqCompareEP(seq, refseq, startpos, stoppos)

[score, algn, seqstarts]=swalign(refseq(startpos:stoppos),seq,'Alphabet','NT','GapOpen', 15,'ExtendGap', 3);

refseqaa=nt2aa(refseq(startpos:stoppos),'AlternativeStartCodons', false); 
% MIND startpos and stoppos need to be in frame for this to work.
% A generalized implementation is not made for now, since the functionality is not currently needed.

startbase=find(algn(1,:)~='-',1,'first');
endbase=find(algn(1,:)~='-',1,'last');

refnum=zeros(1,numel(algn(1,:)));
p=seqstarts(1)+startpos-2;
for i=1:endbase
    if algn(1,i)~='-'
        p=p+1;
    end
    refnum(i)=p;
end

if refnum(end)<stoppos  %%% since swalign truncates the wt seq if the last bases of wt and query don't match, they will need to be put back, otherwise things go wrong later on 
    extraseqlength=stoppos-refnum(end);
    numgapsquery=sum(algn(3,1:end)=='-');
    queryendbase=size(algn,2)-numgapsquery+seqstarts(2)-1;
    algn=[algn char(32*ones(3,extraseqlength))];
    algn(1,end-extraseqlength+1:end)=refseq(refnum(end)+1:stoppos);
    algn(3,end-extraseqlength+1:end)=char(45*ones(1,extraseqlength));
    if queryendbase+1<=length(seq)
        addseqlength=min(extraseqlength,length(seq)-queryendbase);
        algn(3,end-extraseqlength+1:end-extraseqlength+addseqlength)=seq(queryendbase+1:queryendbase+addseqlength);
    end
    refnum=[refnum refnum(end)+1:stoppos];
    startbase=find(algn(1,:)~='-',1,'first');
    endbase=find(algn(1,:)~='-',1,'last');
end

if refnum(1)>startpos  %%% also swalign truncates the wt seq if the first bases of wt and query don't match 
    extraseqlength=refnum(1)-startpos;
    algn=[char(32*ones(3,extraseqlength)) algn];
    algn(1,1:extraseqlength)=refseq(startpos:startpos+extraseqlength-1);
    algn(3,1:extraseqlength)=char(45*ones(1,extraseqlength)); % adds dashes
    if seqstarts(2)~=1 %check for this how swalign numbers seqstarts
        addseqlength=min(extraseqlength,seqstarts(2)-1);
        algn(3,extraseqlength-addseqlength+1:extraseqlength)=seq(seqstarts(2)-addseqlength:seqstarts(2)-1);
    end
    refnum=[startpos:refnum(1)-1 refnum];
    startbase=find(algn(1,:)~='-',1,'first');
    endbase=find(algn(1,:)~='-',1,'last');
end

algnnum=1:size(algn,2);

INS=algn(1,startbase:endbase)=='-';
DEL=algn(3,startbase:endbase)=='-';
SUBS=(algn(2,startbase:endbase)~='|')&(algn(1,startbase:endbase)~='-')&(algn(3,startbase:endbase)~='-');

numINS=sum(INS);
numDEL=sum(DEL);
numSUBS=sum(SUBS);

MUTlocs=logical([zeros(1,startbase-1),INS+DEL+SUBS]); % logical array of size endbase having 1's at mutated positions
SUBSlocs=logical([zeros(1,startbase-1),SUBS]); % logical array of size endbase having 1's at substituted positions

MUTS=[algn(1,MUTlocs);algn(3,MUTlocs)]'; % subset of alignment containing the mutations
MUTpos=refnum(MUTlocs)'; % positions of the mutations in refseq numbering   
SUBSpos=refnum(SUBSlocs)'; % substitution position in refseq numbering
SUBSposalgn=algnnum(SUBSlocs); % substitution position in the alignment
aaSUBSpos=ceil(SUBSpos/3); % substitution position in refseq amino acid numbering
codonSUBSpos=mod(SUBSpos-1,3)+1; % substitution position within the codon

MUTcode=zeros(numel(MUTpos),1);
dummyctr=1:4;
dummyarr=['A','C','G','T'];

for i=1:numSUBS 
    MUTcode(SUBSpos(i)==MUTpos)=dummyctr(MUTS(SUBSpos(i)==MUTpos,2)==dummyarr);
end
    
MUTcodeAA=zeros(size(MUTpos));
classMUT=cell(size(MUTpos));
MUTAAname=zeros(size(MUTpos));
% Sorry, unelegant way to do this:
dummy=1:numel(MUTpos);
dummyins=dummy(MUTS(:,1)=='-');
dummyctr=1:8;
dummyarr=['x','x','x','x','A','C','G','T'];
for i=dummyins
    classMUT{i,1}='insertion';
    MUTcode(i)=dummyctr(MUTS(i,2)==dummyarr);
    MUTcodeAA(i)=3;
    MUTAAname(i)=22;
end
dummydel=dummy(MUTS(:,2)=='-');
for i=dummydel
    classMUT{i,1}='deletion';
    MUTcode(i)=9;
    MUTcodeAA(i)=3;
    MUTAAname(i)=22;
end

for i=1:numSUBS
    if codonSUBSpos(i)==1
        seqcodon=algn(3,SUBSposalgn(i):SUBSposalgn(i)+2);
    elseif codonSUBSpos(i)==2
        seqcodon=algn(3,SUBSposalgn(i)-1:SUBSposalgn(i)+1);
    else
        seqcodon=algn(3,SUBSposalgn(i)-2:SUBSposalgn(i)); 
    end
    if sum(seqcodon=='-')>0
        muttype='gapped'; % this option needs to be there in case there is a substitution near a deletion 
        codeAA=3;
        MUTAAname(SUBSpos(i)==MUTpos)=22;
    else
        mutatedcodon=nt2aa(seqcodon,'AlternativeStartCodons', false);
        if refseqaa(aaSUBSpos(i))==mutatedcodon
            muttype='synonymous';
            codeAA=0;
        else
            muttype=[num2str(refseqaa(aaSUBSpos(i))),num2str(aaSUBSpos(i)),num2str(mutatedcodon)];
            codeAA=1;
            if mutatedcodon=='*'
                codeAA=2;
            end
            MUTAAname(SUBSpos(i)==MUTpos)=find(mutatedcodon==['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','*']); % very QnD assignment of aa mutations
        end

    end
    classMUT{SUBSpos(i)==MUTpos}=muttype;
    MUTcodeAA(SUBSpos(i)==MUTpos)=codeAA;
end
                
        

