%% DISCRETE VERSUS CONTINUOUS MODELS OF, EG, SAI, RECRUITMENT CURVES

% 1. abstract model
% x=input intensity from 0 to 1
% y=output intensity from 0 to 1
% discrete model = {0 from 0 to 0.5; 1 from 0.5 to 1}
% continuous model, x=y
% vary: true correlation between x and y
% measure: t, r, r-squared, p, datapoints required to reach 80% power and p<.05

% constants		
n=20;					% datapoints per sample
xlim=[0,1];				% limits of the input (intensity, time)
ylim=[0,1];				% limits of the output (amplitude)
is=10000;				% N iterations to simulate

summary=nan(is,6,2);			% iterations, parameters (c, m, r, r^2, t, p), models (continuous, discrete)

% set up simulation
X(:,1)=linspace(xlim(1),xlim(2),n)';	% x = equally-spaced points on the x-axis
X(:,2)=round(X)-(round(X(:,1))-0.5)./2;	% 0.25 and 0.75

% start simulation
for i=1:is
    for m=1:2
        Y=X(:,m)+rand(numel(X(:,m)),1)-0.5;			% add random variability from -0.5 to +0.5 to X
        Y=Y-(min(Y));						% subtract minimum value
        Y=Y./(max(Y)-min(Y));					% re-scale data to 0:1
        summary(i,[1,2],m)=[ones(length(X(:,m)),1) X(:,m)]\(Y);% simple linear regression (MUCH quicker than fitlm)
        r=corrcoef(X(:,m),Y);					% correlation matrix
        summary(i,3,m)=r(1,2);					% r-value
    end
end

% summary stats
summary(:,4,:)=summary(:,3,:).^2;				% compute r^2
summary(:,5,:)=summary(:,3,:)./sqrt((1-summary(:,4,:))./(n-2));	% convert r to t
summary(:,6,:)=2.*tcdf(summary(:,5,:),n-2,'upper');		% calculate p

% plot last example
figure(1);
for m=1:2
    subplot(2,3,1+((m-1).*3));		% raw data - continuous
    hold on;
    plot([min(X(:,m)),max(X(:,m))],[min(X(:,m)),max(X(:,m))],'k-');
    plot(X(:,m),Y,'ko');
    plot([min(X(:,m)),max(X(:,m))],[min(X(:,m)),max(X(:,m))].*summary(i,2,m)+summary(i,1,m),'r-');% plot best-fit line
    axis([xlim,ylim]);
    ylabel('Output (eg amplitude)');

    subplot(2,3,2+((m-1).*3));		% histogram of r^2 - continuous
    hold on;
    histogram(summary(:,4,m),100);
    plot([nanmean(summary(:,4,m)),nanmean(summary(:,4,m))],[0,500],'r-');
    axis([0,1,0,500]);

    subplot(2,3,3+((m-1).*3));		% histogram of p - continuous
    hold on;
    histogram(log(summary(:,6,m)),100);
    plot([nanmean(log(summary(:,6,m))),nanmean(log(summary(:,6,m)))],[0,500],'r-');
    axis([-20,0,0,500]);
    xticks(-20:5:0);
end

% final formatting
subplot(2,3,1);
title('example data & fit');
subplot(2,3,2);
title('distribution of r^2');
subplot(2,3,3);
title('distribution of log(p)');
subplot(2,3,4);
xlabel('Input (eg intensity/time)');
ylabel('Output (eg amplitude)');
subplot(2,3,5);
xlabel('Correlation, r^2');
subplot(2,3,6);
xlabel('Correlation, log(p)');
print('discrete_continuous.png','-dpng');
close(1);