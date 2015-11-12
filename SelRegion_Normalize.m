%% This program selects an elliple, defines its boundary, measures
%  the offset as the mean value at the boundary, subtracts the
%  offset and normalizes  w.r.t to the integral.
%  Also takes projections of normalized spectra onto emission and
%  absorption axes and plots them. 

clear all; clc;

Index = 6;
RootPath = strcat('.\3D',num2str(Index),'\2dfftmatrixS1\');
FName = 'MComplexT'; 
EmiFreqName = 'gEmiFreq.dat';
AbsFreqName = 'gAbsFreq.dat';

% Create output folder
if ~isdir(strcat('.\3D',num2str(Index),'\AnalysisResults\NormalizedData'))
    mkdir(strcat('.\3D',num2str(Index),'\AnalysisResults\'),'NormalizedData');
end
OutPath = strcat('.\3D',num2str(Index),'\AnalysisResults\NormalizedData\');

NContourLevels = 20;
FsHeNeFrg = 473.61338;          % HeNe frequency THz
TUdrSmplRatio = 32;             % # of fringes each step T moves
DelayTStep = TUdrSmplRatio / (FsHeNeFrg*2);

% Read the abs and emi frequency grid
gAbsFreq = dlmread(strcat(RootPath, AbsFreqName), '\t');
gEmiFreq = dlmread(strcat(RootPath, EmiFreqName), '\t');

% Define the diag line
DiagLine(1, :) = linspace(gEmiFreq(1, 1), gEmiFreq(1, end), 20);
DiagLine(2, :) = linspace(gAbsFreq(end, 1), gAbsFreq(1, 1), 20);

% The number of rows and cols of the matrix
NAbsDim = size(gAbsFreq, 1);
NEmiDim = size(gAbsFreq, 2);

% Define control parameters
NFiles = 25;            % Max file # in figure labels
T = zeros(NFiles+1,1);	% Define T values
EmiProj = zeros(NFiles+1,NEmiDim);    % Normalized emission projection
AbsProj = zeros(NFiles+1,NAbsDim);    % Normalized emission projection
flag = 1;               % Used to control execution

% Define the region for the ellipse
CEmi = 1544.35;    % Emission center
CAbs = -1544.5;   % Absorption center
DiagL = 3.2;    % Length along diagonal
XDiagL = 2.5;   % Length along cross-diagonal
PkReg = zeros(NAbsDim,NEmiDim);     % Peak region
BkReg = zeros(NAbsDim,NEmiDim);     % Peak boundary
% NPtsPk = 0;       % # of points in peak
% NPtsBk = 0;       % # of points in boundary

Thresh = ((gEmiFreq-CEmi)-(gAbsFreq-CAbs)).^2./(2*(DiagL/2)^2)...
    + ((gEmiFreq-CEmi)+(gAbsFreq-CAbs)).^2./(2*(XDiagL/2)^2) - 1;

% Define region within peak
for j = 1 : NAbsDim
    for k = 1 : NEmiDim
        if Thresh(j,k) <= 0
            % Add location to peak
            PkReg(j,k) = 1;
%             NPtsPk = NPtsPk + 1;
        end
    end
end
NPtsPk = sum(sum(PkReg));       % # of points in peak

% Define peak boundary
for j = 1 : NAbsDim
    for k = 1 : NEmiDim
        if PkReg(j,k) == 1
            if PkReg(j,k+1)==0 || PkReg(j,k-1)==0 || PkReg(j-1,k)==0 || PkReg(j+1,k)==0
                % Add location to boundary
                BkReg(j,k) = 1;
%                 NPtsBk = NPtsBk + 1;
            end
        end
    end
end
NPtsBk = sum(sum(BkReg));       % # of points in boundary

%% Check if the selected region is OK

ChkInd = 0;            % Index for spectra to check ellipse position
ChkName = strcat(RootPath,FName,num2str(ChkInd),'.dat');
ChkComplex = dlmread(ChkName, '\t');
ChkAmp = abs(ChkComplex);
ChkMax = max(max(ChkAmp));
ChkNorm = ChkAmp./ChkMax;

% Find the max position frequency
[r,c] = find(ChkNorm == 1);
PkEmiF = gEmiFreq(r,c);
PkAbsF = gAbsFreq(r,c);

MaxVal = max(max(ChkNorm));

figure(1);
set(gcf, 'Units', 'inch');
set(gcf, 'position', [0.5 1 5 5]);
contourf(gEmiFreq, gAbsFreq, ChkNorm, ...
    linspace(0, MaxVal, NContourLevels), 'LineStyle', 'none');
line(DiagLine(1,:), DiagLine(2,:), 'LineStyle', ':',...
    'Color', [0 0 0]);
colormap('Jet');
hold on;
contour(gEmiFreq, gAbsFreq, BkReg, ...
    linspace(0, MaxVal, NContourLevels), 'LineWidth', 1.5);
hold off;

% Check if the slice position is correct or not
HowProcessBtn = questdlg('Is the slice position correct?',...
    'Yes','No');
if strcmp(HowProcessBtn, 'No')|isempty(HowProcessBtn)
    flag =0;
end;

%% Perform analysis if the correct region is selected

if flag
    for j = 0 : NFiles
        T(j+1) = 0.2 + 60*j*DelayTStep;
        M2DComplex = dlmread(strcat(RootPath,FName,num2str(j),'.dat'),'\t');
        M2DAmpl = abs(M2DComplex);
        % Evaluate background
        M2DAmplBk = M2DAmpl.*BkReg;
        M2DBk = sum(sum(M2DAmplBk))./NPtsBk;
        % Subtract background
        M2DAmpl = M2DAmpl - M2DBk;
        % Select the peak
        M2DAmplCut = M2DAmpl.*PkReg;
        % Divide by the integral
        M2DAmplAvg = sum(sum(M2DAmplCut))./NPtsPk;
        M2DAmplNorm = M2DAmplCut./M2DAmplAvg;
%         % Normalize so that elements are in range [0,1]
%         M2DAmplNorm = M2DAmplNorm./max(max(M2DAmplNorm));
        % Save background corrected normalized spectrum
        dlmwrite(strcat(OutPath,'MNormAmplT',num2str(j),'.dat'),...
            M2DAmplNorm,'\t');
%         % Projection on emission axis
%         EmiProj(j+1,:) = sum(M2DAmplNorm,1);
%         % Projection on absorption axis
%         AbsProj(j+1,:) = sum(M2DAmplNorm,2);
    end
    
    
%     % Plot emission projection
%     figure(2);
%     VMax = max(max(EmiProj));
%     VMin = min(min(EmiProj));
%     contourf(gEmiFreq(1,:),T,EmiProj,...
%         linspace(VMin,VMax,NContourLevels),'LineStyle','none');
%     set(gca,'FontSize',16);
%     set(gca,'YTick',[1 10 20 30 40 50],...
%         'YTickLabel',[' 1'; '10'; '20'; '30'; '40'; '50']);
%     xlabel('Emission energy (meV)','FontSize',20);
%     ylabel('Population time (ps)','FontSize',20);
%     colorbar('Location','EastOutside','FontSize',14);
%     saveas(gcf, strcat(rootPath,'Emi_Projection'), 'emf');
%     
%     % Plot absorption projection
%     figure(3);
%     VMax = max(max(AbsProj));
%     VMin = min(min(AbsProj));
%     contourf(gAbsFreq(:,1),T,AbsProj,...
%         linspace(0,VMax,NContourLevels),'LineStyle','none');
%     set(gca,'FontSize',16);
%     set(gca,'YTick',[1 10 20 30 40 50],...
%         'YTickLabel',[' 1'; '10'; '20'; '30'; '40'; '50']);
%     xlabel('Absorption energy (meV)','FontSize',20);
%     ylabel('Population time (ps)','FontSize',20);
%     colorbar('Location','EastOutside','FontSize',14);
%     saveas(gcf, strcat(rootPath,'Abs_Projection'), 'emf');

end
