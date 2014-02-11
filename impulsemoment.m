function [y,y05,y95] = impulsemoment(IRF)

[nsim,n,H] = size(IRF);
y = squeeze(mean(IRF));
y05 = zeros(n,H);
y95 = zeros(n,H);
for i=1:H
 
  d= sort(squeeze(IRF(:,:,i))); 
  y05(:,i) = d(round(0.05*nsim),:);
  y95(:,i) = d(round(.95*nsim),:);
end
end
