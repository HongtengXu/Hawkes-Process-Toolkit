function Z=SoftThreshold_S( A, thres )


tmp=A;
S=sign(tmp);
tmp=(abs(tmp)-thres);
tmp(tmp<=0)=0;
Z=(S.*tmp);





