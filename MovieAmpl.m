%% Plot data and make video files
clear all; clc;

% The directory and file names of amplitude, real and imaginary data
ScanIndex = '6';
RootPath = strcat('.\3D',ScanIndex,'\');
S1Path = strcat(RootPath, '2dfftmatrixS1\');
% S2Path = strcat(rootPath, '2dfftmatrixS2\');
% Data file
FName = 'MComplexT';
% Frequency axes files
EmiFreqName = 'gEmiFreq.dat';
AbsFreqName = 'gAbsFreq.dat';

% Create output folder
if ~isdir(strcat(RootPath,'AnalysisResults\'))
    mkdir(RootPath,'AnalysisResults');
end
OutPath = strcat(RootPath,'AnalysisResults\');


% Plot variables
NContourLevels = 20;            % Contour levels for plots
SliceOffset = 0.3;              % in meV
FsHeNeFrg = 473.61338;          % HeNe frequency THz
TUdrSmplRatio = 32;             % # of fringes each step T moves

% Read the abs and emi frequency grid
gAbsFreq = dlmread(strcat(S1Path, AbsFreqName), '\t');
gEmiFreq = dlmread(strcat(S1Path, EmiFreqName), '\t');

% Define the diag line
DiagLine(1, :) = linspace(gEmiFreq(1, 1), gEmiFreq(1, end), 20);
DiagLine(2, :) = linspace(gAbsFreq(end, 1), gAbsFreq(1, 1), 20);

% The number of rows and cols of the matrix
NAbsDim = size(gAbsFreq, 1);
NEmiDim = size(gAbsFreq, 2);

fig1 = figure(1);
numframes = 25; 
A = moviein(numframes,fig1); % create the movie matrix 
set(fig1,'NextPlot','replacechildren') 

for m = 0 : numframes
% Read the 2d matrix
% M2DComplexS1 = dlmread(strcat(rootPath, FName, num2str(m), '.dat'), '\t');
% M2DComplexS2 = dlmread(strcat(rootPath, FName, num2str(m), '.dat'), '\t');

M2DComplex = dlmread(strcat(S1Path, FName, num2str(m), '.dat'), '\t');
M2DAbs = abs(M2DComplex);

% M2DComplexS1 = flipud(M2DComplexS1); % flip S1 upside down

% M2DCorr = M2DComplexS1 + M2DComplexS2; % correlation spectrum (S1+S2)

if m==0
    M2DMax = max(max(M2DAbs));     % Define normalization factor
    [RMax,CMax] = find(M2DAbs == M2DMax);
    MaxAbsFreq = gAbsFreq(RMax,1);
end
M2DAbs = M2DAbs/M2DMax;
% VMax = 1;
VMax = max(max(M2DAbs));

% Define time steps
DelayTStep = TUdrSmplRatio / (FsHeNeFrg*2);
T = 0.2 + 60*m*DelayTStep;
% if m < 11
%     T = 12*m*DelayTStep;
% else
%     T = T + 30*DelayTStep;
% end

% Plot 2D spectrum
% fig = figure(31);
set(gcf, 'Units', 'inch');
set(gcf, 'position', [5 1 6 6]);

% contour(gEmiFreq, gAbsFreq, M2DAbs, ...
%     linspace(0, VMax, 1*NContourLevels), 'LineWidth', 1.5);
contourf(gEmiFreq, gAbsFreq, M2DAbs, ...
    linspace(0, VMax, 1*NContourLevels), 'LineStyle', 'none');
line(DiagLine(1,:), DiagLine(2,:), 'LineStyle', ':', 'Color', [0 0 0]);
line([gEmiFreq(1) gEmiFreq(end)],[MaxAbsFreq MaxAbsFreq],'LineStyle','-','Color',[1 1 1]);
line([gEmiFreq(1) gEmiFreq(end)],[MaxAbsFreq-SliceOffset MaxAbsFreq-SliceOffset],...
    'LineStyle',':','Color',[1 1 1]);
line([gEmiFreq(1) gEmiFreq(end)],[MaxAbsFreq+SliceOffset MaxAbsFreq+SliceOffset],...
    'LineStyle','-.','Color',[1 1 1]);
% line(SliceLine(1,:), SliceLine(2,:), 'LineWidth', 2, 'LineStyle', ':', 'Color', 'b');

set(gca,'FontSize',16);
% set(gca,'YTick',[-1546 -1545 -1544 -1543],...
%     'YTickLabel',['-1546'; '-1545'; '-1544'; '-1543']);
xlabel('Emission energy (meV)','FontSize',20);
ylabel('Absorption energy (meV)','FontSize',20);
title(['T = ', num2str(T, 3), ' ps']); 
h = get(gca, 'title');
set(h, 'FontSize', 20, 'Color', [1 1 1], 'FontWeight', 'Bold');
pos = get(h,'Position');
pos(2) = pos(2)-1;
pos(1) = pos(1)+2.5;
colorbar('Location','NorthOutside','FontSize',14);
set(h, 'Position', pos);
colormap('Jet');

A(:,m+1) = getframe(fig1);

end

movie(fig1,A,1,15) % Play the MATLAB movie 
movie2avi(A,strcat(OutPath,'GaAs_4QW_LLLL.avi'),'fps',2.5,'quality',100);
% save movie.mat A % save the MATLAB movie to a file 
% mpgwrite(A,jet,'movie.mpg'); % Convert the movie to MPEG format 
% % Notice the MPEG file is about a quarter of the size of the MATLAB movie file 
% unix('mpeg_play movie.mpg') % Play the MPEG movie
