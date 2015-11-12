function [PkVal EmiFreq AbsFreq] = PeakValsFunc(RangeValues,Scan_Index,Xd)
%UNTITLED2 Summary of this function goes here
%   This is a modification of the PeakVals script.  Basically it finds
% the peaks at a given range of values N1, the scan index is indicated by
% Scan_index (so that you can look at lots of data quickly).

%  Currently the function is set to close all graphs (and some are
%  commented out to make the function run faster.   If you want to see them
%  uncomment them out and on the last line comment out "close all"

% % Define path for 3D scans
% Folder = '.\3D9\2dfftmatrixS1\';
% FName = strcat('MComplexT',num2str(Scan_Index),'.dat');
% % Define path for 2D scans
Folder = strcat('.\2D',num2str(Scan_Index),'\Output\');
FName = strcat('MAmpl1.dat');

gAbsFName = 'gAbsFreq.dat';
gEmiFName = 'gEmiFreq.dat';

gEmiFreq = dlmread(strcat(Folder,gEmiFName),'\t');
gAbsFreq = dlmread(strcat(Folder,gAbsFName),'\t');
MSize = size(gEmiFreq);
MDim = [gEmiFreq(1) gEmiFreq(end) -gEmiFreq(1) -gEmiFreq(end)];
% MDim give the boundaries of emission and absorption energy (meV) axes of the 2D matrix
% in form of [ Emission_low Emission_High Absorption_low Absorption_High ]
% (replace absorption w/ 2 quantum energy if S3)

if MDim(1)*MDim(3) < 0
    isRephasing = 1;   %Rephasing data
    %     MLowLmt = max([ MDim(1) abs(MDim(4)) ]);
    %     MHiLmt = min([ MDim(2) abs(MDim(3)) ]);
    AxisDim = [MDim(1) MDim(2) MDim(4) MDim(3)];
else
    isRephasing = 0;   %Non-Rephasing data
    %     MLowLmt = max([ MDim(1) abs(MDim(3)) ]);
    %     MHiLmt = min([ MDim(2) abs(MDim(4)) ]);
    AxisDim = [MDim(1) MDim(2) MDim(1) MDim(2)];
end

MData = dlmread(strcat(Folder, FName), '\t');
% MData = MData(:,1:MSize(2));
MData = (abs(MData));
%     MData2 = MData.^2;
VMax = max(max(MData));
%     VMax2 = max(max(MData2));

Xpinit1 = RangeValues(1);
Xpinit2 = RangeValues(2);
Ypinit1 = RangeValues(3);
Ypinit2 = RangeValues(4);

% number of contour lines to be used
NContourLevels = 20;
% thickness of contour lines
ContourlineWidth = 1.0;
% % the num of points representing diag or cross-diag line
% LineLength = 200;
% MLowLmt=1595;
% MHiLmt=1620;


% Plot Amp, real and imaginary parts of 2D data
NRow = size(MData, 1);  NCol = size(MData, 2);

%% Plot 2D Amplitude
opengl neverselect;
fig1 = figure(11);
set(fig1, 'Position', [ 50 200 700 660 ]);
%      MData(NRow, NCol-1)= -1; MData(NRow, NCol)= 1;
contour(gEmiFreq, gAbsFreq, MData, linspace(0, VMax, 1*NContourLevels), 'LineWidth', ContourlineWidth);
axis(AxisDim);
axis square;  colormap jet;

DiagX = linspace( MLowLmt, MHiLmt, LineLength );
DiagY = -DiagX;
hold on

%% Determine Integration Region
figure(fig1);
datacursormode on;
dcm_obj = datacursormode(fig1);
set(dcm_obj,'DisplayStyle','window','SnapToDataVertex','off');
% CurInfo = getCursorInfo(dcm_obj);
% [Xp Yp] = CurInfo.Position;
% mouse cursor point at any random position
if ~Xd
    [Xp, Yp, MouseButton]= ginput(4);
%     if MouseButton ~=1
%         break;
%     end
else
    Xp(1) = Xpinit1; Xp(4) = Xpinit1; Xp(3) = Xpinit2; Xp(2) = Xpinit2;
    Yp(1) = Ypinit1; Yp(2) = Ypinit1; Yp(3) = Ypinit2; Yp(4) = Ypinit2;
end

%% Calculate and plot lines with values above as corners
htoplength = linspace(Xp(1),Xp(2),LineLength);
hbottomlength = linspace(Xp(4),Xp(3),LineLength);
vleftlength = linspace(Yp(4),Yp(1),LineLength);
vrightlength = linspace(Yp(3),Yp(2),LineLength);

topslope = (Yp(2)-Yp(1))/(Xp(2)-Xp(1));
bottomslope = (Yp(3)-Yp(4))/(Xp(3)-Xp(4));
leftslope = (Yp(1)-Yp(4))/(Xp(1)-Xp(4));
rightslope = (Yp(3)-Yp(2))/(Xp(3)-Xp(2));

vtoplength = Yp(1)+topslope.*(htoplength-htoplength(1));
vbottomlength = Yp(4)+bottomslope.*(hbottomlength-hbottomlength(1));
hleftlength = Xp(4)+(vleftlength-vleftlength(1))./leftslope;
hrightlength = Xp(3)+(vrightlength-vrightlength(1))./rightslope;

line(htoplength, vtoplength, 'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);
line(hbottomlength, vbottomlength, 'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);
line(hleftlength, vleftlength, 'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);
line(hrightlength, vrightlength, 'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);
hold off

%% Determine indices in frequency and amplitude data for region within BOX
i=1;
while gEmiFreq(1,i) < Xp(1)
    i = i+1; % lower index for X integration
end

j=1;
while gEmiFreq(1,j) < Xp(2)
    j = j+1; % upper index for X integration
end

xindex = j-i;

k=1;
while gAbsFreq(k,1) < Yp(4)
    k = k+1; % lower index for Y integration
end

l=1;
while gAbsFreq(l,1) < Yp(1)
    l = l+1; % upper index for Y integration
end

yindex = l-k;

%% Meshgrid for region within box
[gEmiFreqInt, gAbsFreqInt] = meshgrid( single(linspace(Xp(1), Xp(2), xindex+1)), ...
    single(linspace(Yp(4), Yp(1), yindex+1)));

% Isolate region of amplitude data within box chosen by indices at
% beginning and plot
MDataInt = MData(k:l,i:j);
figure(14)
% MData(NRow, NCol-1)= -1; MData(NRow, NCol)= 1;
% VMax = 500;
contour(gEmiFreqInt, gAbsFreqInt, MDataInt, linspace(0, VMax, NContourLevels), 'LineWidth', ContourlineWidth);
axis([Xp(1) Xp(2) Yp(4) Yp(1)]);
axis square;  colormap jet;

%% Find position of max

i=1;
j=1;

while MDataInt(i,j) < max(max(MDataInt))
    for j=1:(size(MDataInt,2))
        if MDataInt(i,j) == max(max(MDataInt))
            break
        end
    end
    
    if MDataInt(i,j) == max(max(MDataInt))
        break
    end
    
    i=i+1;
    if i == (size(MDataInt,1))
        break
    end
end

format long
xpeakvalue = gEmiFreqInt(1,j);
ypeakvalue = gAbsFreqInt(i,1);

% Find position of exciton peak
i=1;
j=1;
while gEmiFreq(1,i)<xpeakvalue
    i=i+1;
end
while gAbsFreq(j,1)<ypeakvalue
    j=j+1;
end
PkX = i;
PkY = j;
PkVal = MData(PkY,PkX);
EmiFreq = gEmiFreq(PkY,PkX);
AbsFreq = gAbsFreq(PkY,PkX);