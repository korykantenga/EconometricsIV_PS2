
%=========================================================================
%           TABLE 2: MEANS & CREDIBLE SETS
%=========================================================================

%Means
meanPhi_0      = mean(Phip(:,17,:),1);
meanLargeEigen = mean(largeeig);

tempPhi_0        = zeros(size(Phip,1),size(Phip,3));
tempPhiMean      = zeros(1,4);
tempPhiMean(1,:) = meanPhi_0(1,1,:);

for q=1:size(Phip,3)
    tempPhi_0(:,q) = Phip(:,17,q);
end

%Credible Sets
CSetPhi_0      = hpdi(tempPhi_0,90);
CSetLargeEign  = hpdi(largeeig,90);

%Latex Code
CSets          = [CSetPhi_0'; CSetLargeEign'];
means          = [tempPhiMean';meanLargeEigen];
disp('Latex Code for Table 2: Means & Credible Sets');
for j=1:size(CSets,1)
    fprintf('%8.3f %8.3f & %8.3f\\\\ \n',...
        means(j,1), CSets(j,1), CSets(j,2))
end
