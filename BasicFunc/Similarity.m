function weight = Similarity(f1, f2, sigma)

if isempty(f1) || isempty(f2)
    weight = 1;
else
    weight = exp(-norm(f1(:)-f2(:))^2/(2*sigma^2));
end