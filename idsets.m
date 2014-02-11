
%=========================================================================
%         IDENTIFIED SETS: IRF MONETARY SHOCK ON R AT HORIZON 1
%=========================================================================

clc

% Sign restrictions say AA*q row2,row4 are NON-POS, row3 ie NON-NEG

qsim   = 100000;
Qgrid  = zeros(n,qsim);
Rgrid  = zeros(n,qsim);
QgridV = zeros(1,qsim);
AA     = chol(Sigma)' ;

for s=1:qsim
    e           = randn(n,1)   ;
    Qgrid(:,s)  = (e./norm(e)) ;
    Rgrid(:,s)  = AA*Qgrid(:,s);
    if Rgrid(2,s)>0 || Rgrid(4,s)>0 || Rgrid(3,s)<0
        QgridV(s)=1;
    end
end

for j = 1:qsim
    if QgridV(j)==0
        Qgrid(:,j)=Qgrid(:,j);
        Rgrid(:,j)=Rgrid(:,j);
        
    else
        Qgrid(:,j)=NaN;
        Rgrid(:,j)=NaN;
    end
end

Qgrid = Qgrid(isnan(Qgrid)==0);
Rgrid = Rgrid(isnan(Rgrid)==0);
Qgrid = reshape(Qgrid,n,length(Qgrid)/n);
Rgrid = reshape(Rgrid,n,length(Rgrid)/n);

ttheta = zeros(1,size(Qgrid,2));
ttheta = Rgrid(3,:);

CSet = [min(ttheta) max(ttheta)];


% Compute Impulse Responses using q

IRF_ID = zeros(n,H,size(Qgrid,2));

for j=1:size(Qgrid,2)
for i=1:H
    IRF_ID(:,i,j) = impulseresponce(Phi_tilde,Sigma,Qgrid(:,j),i);
end
end

[IRF_MIN, I] = min(IRF_ID,[],3);
[IRF_MAX, I] = max(IRF_ID,[],3);

% Plots

figure;
subplot(2,2,1)
hold on
p1 = plot(1:H,IRF_MIN(1,:),'-');
p2 = plot(1:H,IRF_MAX(1,:),'-');
title('ln(GDP per Capita)')
hline = refline(0,0);
xlim([0 H])
set(hline,'Color','k','LineStyle','--')
set(p1,'Color','b')
set(p2,'Color','b')

subplot(2,2,2)
hold on
p1 = plot(1:H,IRF_MIN(2,:),'-');
p2 = plot(1:H,IRF_MAX(2,:),'-');
title('Inflation')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')
set(p1,'Color','b')
set(p2,'Color','b')

subplot(2,2,3)
hold on
p1 = plot(1:H,IRF_MIN(3,:),'-');
p2 = plot(1:H,IRF_MAX(3,:),'-');
title('Interest Rate')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')
set(p1,'Color','b')
set(p2,'Color','b')

subplot(2,2,4)
hold on
p1 = plot(1:H,IRF_MIN(4,:),'-');
p2 = plot(1:H,IRF_MAX(4,:),'-');
title('Money Balances')
xlim([0 H])
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')
set(p1,'Color','b')
set(p2,'Color','b')

print -depsc2 VAR_IRF_Q2_6.eps


% Histogram Plots

figure;
subplot(2,2,1)
hold on
hist(IRF_ID(1,:),50)
title('ln(GDP per Capita) IR at Horizon 0')
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')

subplot(2,2,2)
hold on
hist(IRF_ID(2,:),50)
title('Inflation IR at Horizon 0')
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')

subplot(2,2,3)
hold on
hist(IRF_ID(3,:),50)
title('Interest Rate IR at Horizon 0')
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')
plot([CSet(1), CSet(1)],[0 20000],'b')
plot([CSet(2), CSet(2)],[0 20000],'b')


h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')

subplot(2,2,4)
hold on
hist(IRF_ID(4,:),50)
title('Money Balances IR at Horizon 0')
hline = refline(0,0);
set(hline,'Color','k','LineStyle','--')

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')

print -depsc2 VAR_IRF_Q2_7.eps


