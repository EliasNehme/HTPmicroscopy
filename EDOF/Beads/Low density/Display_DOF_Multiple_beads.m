% This code centers the beads according to their axial position, by finding
% the minimum of their 1D gaussian fit for sigma(z), and averages the
% sigma(z) values to estimate the PSF widths.

dz = 1; % axial step of the zstacks, in um (keep as 1).

% colors for plots
colorEDOF = [0.9686    0.2000    0.2843];
colorST = [0.4235    0.6280    0.9912];

% Load 2d gaussian fit results (fit per frame around a frame having the maximal pixel value)
load('fittedST');
% fittedST = fitted;
FS = F_ST;
zST = (-FS:FS)*dz;
load('fittedEDOF');
% fittedEDOF = fitted;
FE = F_EDOF;
zEDOF = (-FE:FE)*dz;

umPerPix = 5.86/5; % um per pixel

% reshape data - each bead will correspond a row
HS = reshape(fittedST(:,3),[(FS+FS+1),size(fittedST,1)/(FS+FS+1)])';
HE = reshape(fittedEDOF(:,3),[(FE+FE+1),size(fittedEDOF,1)/(FE+FE+1)])';

 % Remove beads having bad localizations
HS(max(HS,[],2)>13,:) = [];
HE(max(HE,[],2)>13,:) = [];

% fit 1D gaussian to each row
for i = 1:size(HS,1)
[xData, yData] = prepareCurveData( zST, HS(i,:) );

% Set up fittype and options.
ft = fittype( 'b-a*exp(-(z-z0)^2/s^2)', 'independent', 'z', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% Startpoint - [a, b, s, z0]
opts.StartPoint = [1 10 5 0];

% Fit model to data.
fitresultS = fit( xData, yData, ft, opts );
centerST(i) = fitresultS.z0;
end

% fit 1D gaussian to each row
for i = 1:size(HE,1)
[xData, yData] = prepareCurveData( zEDOF, HE(i,:) );

% Set up fittype and options.
ft = fittype( 'b-a*exp(-(z-z0)^2/s^2)', 'independent', 'z', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% Startpoint - [a, b, s, z0]
opts.StartPoint = [1 10 5 0];

% Fit model to data.
fitresultE = fit( xData, yData, ft, opts );
centerEDOF(i) = fitresultE.z0;
% figure(8); plot(zEDOF,HE(i,:),zEDOF,fitresultE(zEDOF));
% pause(.1)
end
%
HE1 = HE;
for i = 1:length(centerEDOF)
    shiftAmount = round(centerEDOF(i));
    HE1(i,:) = circshift(HE(i,:),shiftAmount);
    if shiftAmount<0
        HE1(i,1+end+shiftAmount:end) = nan;
    elseif shiftAmount>0
        HE1(i,1:shiftAmount) = nan;
    end
end

%Plot

figure(6); clf
s1 = mean(HS,'omitnan')*umPerPix; % s1 - ST PSF widths in microns
s2 = mean(HE1,'omitnan')*umPerPix; % s2 - EDOF PSF widths in microns
plot(zST,s1,'linewidth',4,'color',colorST); hold on
plot(zEDOF,circshift(s2,3),'linewidth',4,'color',colorEDOF);
legend('Standard PSF','EDOF')
xlim([-40 40]);
ylim([0 9]);
xlabel('z (\mum)');
xlabel('\sigma (\mum)');
set(gca,'fontsize',14)

save('multi_bead_fits_results','s1','s2','zST','zEDOF','HS','HE');

%% plot DOF extension per PSF width
load('multi_bead_fits_results.mat')

focalDefinition = min(s1+0.1):0.1:4;

ratio = zeros(1,length(focalDefinition));
for i = 1:length(focalDefinition)
    DOF1(i) = sum(s1<focalDefinition(i));
    DOF2(i) = sum(s2<focalDefinition(i));
end
ratio = DOF2./DOF1;
xlabel('z (\mum)')
ylabel('\sigma (\mum)')
set(gca,'fontsize',18);
figure(7); clf
plot(focalDefinition, ratio,'r.','MarkerSize',20); ylim([.5 2.5]);
xlabel('\sigma_{focus} (\mum)');
ylabel('DOF ratio (EDOF/no mask)');
set(gca,'fontsize',14)
grid on
%%
focalDefinition = min(s1+0.1):0.1:4;

ratio = zeros(1,length(focalDefinition));
for i = 1:length(focalDefinition)
    DOF1(i) = sum(s1<focalDefinition(i));
    DOF2(i) = sum(s2<focalDefinition(i));
end
ratio = DOF2./DOF1;
xlabel('z (\mum)')
ylabel('\sigma (\mum)')
set(gca,'fontsize',18);
figure(7); clf
plot(focalDefinition*(5.86/5), ratio,'r*'); ylim([.5 2.5]);
xlabel('\sigma_{focus} (\mum)');
ylabel('DOF ratio (EDOF/no mask)');
set(gca,'fontsize',14)
