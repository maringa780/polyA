function rho=corrpercolumn(a,b)

rho = nan(size(a,2),1);

for k=1:numel(rho)
rho(k) =corr(a(:,k),b(:,k));
end

