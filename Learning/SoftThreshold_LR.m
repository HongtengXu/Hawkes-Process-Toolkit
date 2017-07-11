function Z=SoftThreshold_LR( A, thres )


Z = zeros(size(A));

for t = 1:size(A,2)
    tmp = A(:,t,:);
    tmp = reshape(tmp, [size(A,1), size(A,3)]);
    [Ut, St, Vt]=svd(tmp);
    St=St-thres;
    St(St<0)=0;
    Z(:,t,:)=reshape(Ut*St*Vt', [size(A,1), 1, size(A,3)]);
end


