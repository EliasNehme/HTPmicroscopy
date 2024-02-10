% This script finds isolated beads from each dataset (Standard PSF, EDOF),
% crops and localized them using a 2D gaussian.
%
% Per-bead crops are saved as matlab variables crops_XX, and results of 2D
% gaussian fits are saved as matlab variables fittedXX.
%
% This code takes a long time to run due to many thousands of gaussian
% fits. The results are provided as the .mat files.
%
% "ST" stands for Standard PSF data.

% Get data
SDir = dir('ST\*.tif');
EDir = dir('EDOF\*.tif');

% Define sizes for analysis
W = 10; % size of crop around the bead - size will be 2*W+1
F_ST = 40; % Analysis is performed over 'F_ST' frames in each direction, ST 
F_EDOF = 60; %  Analysis is performed over 'F_EDOF' frames in each direction

doSaveResults = 1;

if doSaveResults
    disp('Notice - cropped data and localizations are saved, watch for overrides!');
else
    disp('Notice - cropped data and localizations are not saved');
end

% Read all frames as stacks
if ~exist('EStack0','var')
    for i = 1:400
        SStack0(:,:,i) = imread([SDir(i).folder '\' SDir(i).name]);
        EStack0(:,:,i) = imread([EDir(i).folder '\' EDir(i).name]);
    end
end

% Read tables 
if ~exist('tS','var')
    % where there are beads (Thunderstorm over max projections)
tS = readtable('ThunStormMAXST.csv');
tE = readtable('ThunStormMAXEDOF.csv');
    % tables of which beads to exclude (after applying minimal density filter via
    % thunderstorm)
excludeS = readtable('ThunStormMAXST_close_locs.csv');
excludeE = readtable('ThunStormMAXEDOF_close_locs.csv');
end

nm2pix = 5/5.86*1e-3; % nanometer to pixel conversion factor

% Get pixel positions of the thunderstorm localizations (s - standard, e - EDOF)
xs = round(tS.x_nm_*nm2pix);
ys = round(tS.y_nm_*nm2pix);
xe = round(tE.x_nm_*nm2pix);
ye = round(tE.y_nm_*nm2pix);

% remove ROIs that are too close together (given by Thunderstorm density
% filter)
xs(ismember(tS.id,excludeS.id)) = [];
ys(ismember(tS.id,excludeS.id)) = [];
xe(ismember(tE.id,excludeE.id)) = [];
ye(ismember(tE.id,excludeE.id)) = [];

% remove ROIs near the boundary - standard PSF
removes = double(xs<10*W);
removes = removes + double(xs>1400-10*W);
removes = removes + double(ys<10*W);
removes = removes + double(ys>1040-10*W);
xs = xs(removes==0);
ys = ys(removes==0);

% remove ROIs near the boundary - EDOF
removes = double(xe<10*W);
removes = removes + double(xe>1400-10*W);
removes = removes + double(ye<10*W);
removes = removes + double(ye>1040-10*W);
xe = xe(removes==0);
ye = ye(removes==0);

% Crop Standard-PSF beads
totalST = zeros(0);
tic
for i = 1:length(xs)
    croppedS = SStack0(ys(i)-W:ys(i)+W,xs(i)-W:xs(i)+W,:);
    
    maxesS = squeeze(max(max(croppedS)));
    
    [~,frameS(i)] = max(maxesS);

% Add
if abs(frameS(i)-200)<(199-F_ST)
    croppedS = croppedS(:,:,frameS(i)-F_ST:frameS(i)+F_ST);
    totalST = cat(3,totalST,croppedS);
end

end
disp(['time for ST crop: ' num2str(toc) 's'])

% Crop EDOF beads
totalEDOF = zeros(0);
tic
for i = 1:length(xe)
    croppedE = EStack0(ye(i)-W:ye(i)+W,xe(i)-W:xe(i)+W,:);
    
    maxesE = squeeze(max(max(croppedE)));
    
    [~,frameE(i)] = max(maxesE);

% Add
if abs(frameE(i)-200)<(199-F_EDOF)
    croppedE = croppedE(:,:,frameE(i)-F_EDOF:frameE(i)+F_EDOF);
    totalEDOF = cat(3,totalEDOF,croppedE);
end

end
disp(['time for EDOF crop: ' num2str(toc) 's'])

% fit 2D gaussian per frame - ST
[X,Y] = meshgrid(-W:W);
figure(5); clf
for i = 1:size(totalST,3)
     ROI = double(totalST(:,:,i));
        
        % Initial guess for bg is the mean of corner photons
        getBGImage = fftshift(ROI);
        bg0 = mean(getBGImage(1:4,1:4),'all');
        
        % Initial guess for Nph is the sum of all photons, minus initial bg guess
        Nph0 = max(300,sum(ROI(:))-bg0*size(ROI,1)*size(ROI,2));
        
        centerRow = 0; centerCol = 0;
        
        % Define functions
        % PSF generation
        MakeModelIm = @(x,y,s,N,bg) N/2/pi/s^2*exp(-(((x-X).^2+(y-Y).^2)/2/s^2))+bg;
        % cost function for fmincon
        LSCost = @(gParams) sum((MakeModelIm(gParams(1),gParams(2),gParams(3),gParams(4),gParams(5))-ROI).^2,'all');
        
        % Define initial guess and limits for fmincon
        theta0 = [centerRow centerCol 8 Nph0 bg0];
        
        % Search bounds - defined by image size
        l1 = size(ROI,1)/10;
        l2 = size(ROI,2)/10;
        
        lb = [centerCol-round(l1) centerRow-round(l2) 0.5 Nph0/500 0];
        ub = [centerCol+round(l1) centerRow+round(l2) 15 Nph0*100 sum(ROI(:))];
        
        options = optimoptions(@fmincon,'Display','iter','algorithm','interior-point','tolX',1e-13,'TolFun',1e-13,'display','off');
        
        % fit gaussian
        [gFit] = fmincon(LSCost,theta0,[],[],[],[],lb,ub,[],options);
        
        fittedST(i,:) = gFit;
        if ~mod(i,100)
             plot(fittedST(:,3));
             title('Fitted \sigma per frame, Standard')
             drawnow
        end
end
% 
% 
% % fit 2D gaussian per frame - ST
% [X,Y] = meshgrid(-W:W);
% figure(5); clf
% for i = 1:size(totalST,3)
%      ROI = double(totalST(:,:,i));
%         
%         % Initial guess for bg is the mean of corner photons
%         getBGImage = fftshift(ROI);
%         bg0 = mean(getBGImage(1:4,1:4),'all');
%         
%         % Initial guess for Nph is the sum of all photons, minus initial bg guess
%         Nph0 = max(300,sum(ROI(:))-bg0*size(ROI,1)*size(ROI,2));
%         
%         centerRow = 0; centerCol = 0;
%         
%         % Define functions
%         % PSF generation
%         MakeModelIm = @(x,y,s,N,bg) N/2/pi/s^2*exp(-(((x-X).^2+(y-Y).^2)/2/s^2))+bg;
%         % cost function for fmincon
%         LSCost = @(gParams) sum((MakeModelIm(gParams(1),gParams(2),gParams(3),gParams(4),gParams(5))-ROI).^2,'all');
%         
%         % Define initial guess and limits for fmincon
%         theta0 = [centerRow centerCol 8 Nph0 bg0];
%         
%         % Search bounds - defined by image size
%         l1 = size(ROI,1)/10;
%         l2 = size(ROI,2)/10;
%         
%         lb = [centerCol-round(l1) centerRow-round(l2) 0.5 Nph0/500 0];
%         ub = [centerCol+round(l1) centerRow+round(l2) 15 Nph0*100 sum(ROI(:))];
%         
%         options = optimoptions(@fmincon,'Display','iter','algorithm','interior-point','tolX',1e-13,'TolFun',1e-13,'display','off');
%         
%         % fit gaussian
%         [gFit] = fmincon(LSCost,theta0,[],[],[],[],lb,ub,[],options);
%         
%         fittedST(i,:) = gFit;
%         if ~mod(i,200)
%              plot(fittedST(:,3));
%              drawnow
%         end
% end

% fit 2D gaussian per frame - EDOF
figure(5); clf
for i = 1:size(totalEDOF,3)
     ROI = double(totalEDOF(:,:,i));
        
        % Initial guess for bg is the mean of corner photons
        getBGImage = fftshift(ROI);
        bg0 = mean(getBGImage(1:4,1:4),'all');
        
        % Initial guess for Nph is the sum of all photons, minus initial bg guess
        Nph0 = max(300,sum(ROI(:))-bg0*size(ROI,1)*size(ROI,2));
        
        centerRow = 0; centerCol = 0;
        
        % Define functions
        % PSF generation
        MakeModelIm = @(x,y,s,N,bg) N/2/pi/s^2*exp(-(((x-X).^2+(y-Y).^2)/2/s^2))+bg;
        % cost function for fmincon
        LSCost = @(gParams) sum((MakeModelIm(gParams(1),gParams(2),gParams(3),gParams(4),gParams(5))-ROI).^2,'all');
        
        % Define initial guess and limits for fmincon
        theta0 = [centerRow centerCol 8 Nph0 bg0];
        
        % Search bounds - defined by image size
        l1 = size(ROI,1)/10;
        l2 = size(ROI,2)/10;
        
        lb = [centerCol-round(l1) centerRow-round(l2) 0.5 Nph0/500 0];
        ub = [centerCol+round(l1) centerRow+round(l2) 15 Nph0*100 sum(ROI(:))];
        
        options = optimoptions(@fmincon,'Display','iter','algorithm','interior-point','tolX',1e-13,'TolFun',1e-13,'display','off');
        
        % fit gaussian
        [gFit] = fmincon(LSCost,theta0,[],[],[],[],lb,ub,[],options);
        
        fittedEDOF(i,:) = gFit;
        if ~mod(i,200)
             plot(fittedEDOF(:,3));
             title('Fitted \sigma per frame, EDOF')
             drawnow
        end
end


if doSaveResults
    save('crops_ST','totalST','F_ST','F_EDOF','W');
    save('crops_EDOF','totalEDOF','F_ST','F_EDOF','W');
    % save('fittedST','fittedST','F_ST','F_EDOF','W','lb','ub');
    save('fittedST','fittedST','F_ST','W','lb','ub');
    save('fittedEDOF','fittedEDOF','F_ST','F_EDOF','W','lb','ub');
    Display_DOF_Multiple_beads
end
        
        
