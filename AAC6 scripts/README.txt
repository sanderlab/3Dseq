Scripts for filtering and analyzing fastq sequences for AAC6 experimental evolution.

1 - extract Q15 filtered sequences directly from fastq (linux)
    a) ExtractQualAndSeqQ15newT.bat for libraries exclusively in fwd direction
    b) ExtractQualAndSeqFwdQ15newT.bat and ExtractQualAndSeqRevQ15newT.bat for libraries in both directions, improving sequencing quality.
2 - process the Q15 filtered libraries further, and perform mutation assessment, FuncIllumProc.m (Matlab)
3 - inferred contact filtering, scriptMaskContactMap.m (Matlab) 

All other files in the directory are functions called by the main scripts.


Mind: uses third-party matlab function 'fastsmooth.m':

fastsmooth(Y,w,type,ends) smooths vector Y with smooth 
of width w. Version 3.0, October 2016.
Copyright (c) 2012, Thomas C. O'Haver