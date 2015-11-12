%% The program evaluates the change in signal at certain absorption
% band. This is done to map out redistribution of exciton population
% with T. We select a few horizontal slices around the frequency where
% the "new" peak appears. We find the change in the slice as a
% function of T.

clear all; clc;
Index = 6;
RootPath = strcat('.\3D',num2str(Index),'\AnalysisResults\');
NormDataPath = strcat(RootPath,'NormalizedData\');
DiffDataPath = strcat(RootPath,'AmpDiffData\');
% S1Path = strcat(rootPath, '2dfftmatrixS1\');
% S2Path = strcat(rootPath, '2dfftmatrixS2\');
FName = 'MNormAmplT'; 
EmiFreqName = 'gEmiFreq.dat';
AbsFreqName = 'gAbsFreq.dat';

NContourLevels = 20;
FsHeNeFrg = 473.61338;          % HeNe frequency THz
TUdrSmplRatio = 32;             % # of fringes each step T moves
DelayTStep = TUdrSmplRatio / (FsHeNeFrg*2);

% Read the abs and emi frequency grid
gAbsFreq = dlmread(strcat(NormDataPath, AbsFreqName), '\t');
gEmiFreq = dlmread(strcat(NormDataPath, EmiFreqName), '\t');

% Define the diag line
DiagLine(1, :) = linspace(gEmiFreq(1, 1), gEmiFreq(1, end), 20);
DiagLine(2, :) = linspace(gAbsFreq(end, 1), gAbsFreq(1, 1), 20);

% The number of rows and cols of the matrix
NAbsDim = size(gAbsFreq, 1);
NEmiDim = size(gAbsFreq, 2);

NFiles = 25;                         % Max file # in figure labels
SliceDiff1 = zeros(NFiles,NEmiDim);  % Variable to store differential slices
SliceDiff2 = zeros(NFiles,NEmiDim);  % Variable to store differential slices
T = zeros(NFiles,1);                % Define T values

% Define the slices of interest
SliceFreq1 = -1544.9;      % High absorption energy slice
SliceFreq2 = -1544.3;      % Low absorption energy slice
NSlice = 5;     % Average over multiple slices
Check = 1;      % Set to 1 to check slice position
flag = 1;

% Find the positions of the reference slices
AbsAxis = gAbsFreq(:,1);
% First slice
AbsDiff1 = abs(AbsAxis - SliceFreq1);   
Pos1 = find(AbsDiff1 == min(AbsDiff1));     % Matrixindex
% Second slice
AbsDiff2 = abs(AbsAxis - SliceFreq2);
Pos2 = find(AbsDiff2 == min(AbsDiff2));     % Matrixindex

%% Check if the reference line is correct or not
if Check
    % Check slice position
    ChkInd = 13;
    ChkName = strcat(NormDataPath,FName,num2str(ChkInd),'.dat');
    ChkComplex = dlmread(ChkName, '\t');
    ChkAmp = abs(ChkComplex);
    ChkMax = max(max(ChkAmp));
%%    
    figure(1);
    set(gcf, 'Units', 'inch');
    set(gcf, 'position', [0.5 1 5 5]);
%     contour(gEmiFreq, gAbsFreq, ChkAmp, ...
%         linspace(0, ChkMax, 1*NContourLevels), 'LineWidth', 1.5);
    contourf(gEmiFreq, gAbsFreq, ChkAmp, ...
        linspace(0, ChkMax, NContourLevels), 'LineStyle', 'none');
    line(DiagLine(1,:), DiagLine(2,:), 'LineStyle', ':',...
        'Color', [0 0 0]);
    xlabel('Emission energy (meV)','FontSize',16);
    ylabel('Absorbtion energy (meV)','FontSize',16);
    line([gEmiFreq(1) gEmiFreq(end)],[AbsAxis(Pos1) AbsAxis(Pos1)],...
        'LineStyle','-','Color',[0 0 0]);
    line([gEmiFreq(1) gEmiFreq(end)],[AbsAxis(Pos2) AbsAxis(Pos2)],...
        'LineStyle','-','Color',[0 0 0]);
%     set(gca,'XTick',1543:1:1546,'XTickLabel',...
%         {'1543','1544','1545','1546'},'FontSize',12);
%     set(gca,'YTick',-1546:1:-1543,'XTickLabel',...
%         {'-1546','-1545','-1544','-1543'},'FontSize',12);
    
    % Check if the slice position is correct or not
    HowProcessBtn = questdlg('Are the slice positions correct?',...
        'Yes','No');
    if strcmp(HowProcessBtn, 'No')|isempty(HowProcessBtn)
        flag =0;
    end;
end

if flag
    %% Get the slice average value for T = 0
    M2DComplex1 = dlmread(strcat(NormDataPath,FName,'0.dat'),'\t');
    M2DAmpl1 = abs(M2DComplex1);
    
    SliceCut1 = M2DAmpl1(Pos1-(NSlice-1)/2:Pos1+(NSlice-1)/2,:);    % Select region to average
    BkgdSliceAvg1 = sum(SliceCut1,1)./NSlice;       % Average slice

    SliceCut2 = M2DAmpl1(Pos2-(NSlice-1)/2:Pos2+(NSlice-1)/2,:);    % Select region to average
    BkgdSliceAvg2 = sum(SliceCut2,1)./NSlice;   % Average slice
    
    % figure(2);
    % plot(gEmiFreq(1,:),BkgdSliceAvg);
    
    %% Substract slice value for T=0 from the slices for different T's
    for j = 1:NFiles
        T(j) = 0.2 + 60*j*DelayTStep;
     
        M2DComplex = dlmread(strcat(NormDataPath,FName,num2str(j),'.dat'),'\t');
        M2DAmpl = abs(M2DComplex);
        SliceCut1 = M2DAmpl(Pos1-(NSlice-1)/2:Pos1+(NSlice-1)/2,:);
        SliceAvg1 = sum(SliceCut1,1)./NSlice;
        SliceDiff1(j,:) = SliceAvg1 - BkgdSliceAvg1;

        SliceCut2 = M2DAmpl(Pos2-(NSlice-1)/2:Pos2+(NSlice-1)/2,:);
        SliceAvg2 = sum(SliceCut2,1)./NSlice;
        SliceDiff2(j,:) = SliceAvg2 - BkgdSliceAvg2;
    end

    %%
    % Max and min values in "Difference" spectra
    VMax1 = max(max(SliceDiff1));
    VMin1 = min(min(SliceDiff1));
    VMax2 = max(max(SliceDiff2));
    VMin2 = min(min(SliceDiff2));
    
    if VMax1 > VMax2
        VMax = VMax1;
    else
        VMax = VMax2;
    end
    
    if VMin1 < VMin2
        VMin = VMin1;
    else
        VMin = VMin2;
    end
    
    SliceDiff1(1) = VMax;    SliceDiff1(end) = VMin;
    SliceDiff2(1) = VMax;    SliceDiff2(end) = VMin;
    
    figure(3);
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [2 2 17 9]);
    subplot (121);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [2 2 5 5]);
    contourf(gEmiFreq(1,:),T,SliceDiff1,...
        linspace(VMin,VMax,NContourLevels),'LineStyle','none');
    line([-SliceFreq1 -SliceFreq1], [T(1) T(end)],'Color',[0 0 0],'LineWidth',1);
    colormap('Jet');
    set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(gca,'YTick',[3 10 20 30 40 50],...
        'YTickLabel',['03'; '10'; '20'; '30'; '40'; '50']);
    xlabel('Emission energy (meV)');
    ylabel('Population time (ps)');
    title('High Absorption Energy Slice');
    c = colorbar('North');
%     set(c,'XTick',[0 0.03 0.06 0.09 1.20]);
    x=get(c,'Position');
    x(4)=0.5*x(4);
    x(2)=x(2)+0.25;
    set(c,'Position',x);
    hold on;
    
%     figure(3);
    subplot (122);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [10 2 5 5]);
    contourf(gEmiFreq(1,:),T,SliceDiff2,...
        linspace(VMin,VMax,NContourLevels),'LineStyle','none');
    line([-SliceFreq2 -SliceFreq2], [T(1) T(end)],'Color',[0 0 0],'LineWidth',1);
    set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(gca,'YTick',[3 10 20 30 40 50],...
        'YTickLabel',['03'; '10'; '20'; '30'; '40'; '50']);
    xlabel('Emission energy (meV)');
    ylabel('Population time (ps)');
    title('Low Absorption Energy Slice');
    c = colorbar('North');
%     set(c,'XTick',[0 0.03 0.06 0.09 1.20]);
    x=get(c,'Position');
    x(4)=0.5*x(4);
    x(2)=x(2)+0.25;
    set(c,'Position',x);
    hold off;
%     saveas(gcf, strcat(OutPath,'SliceDiff'), 'emf');
    
    %% Quantify asymmetry
    
    %  Offset frequency and flipping axis
    CenterF = -(SliceFreq1 + SliceFreq2)/2;
    Offset = gEmiFreq(1,:) - CenterF;
    OffsetFlip = fliplr(Offset);
    for j = 1 : NEmiDim-1
        if Offset(j)*Offset(j+1) <= 0
            Ref1 = j;
        end
        if OffsetFlip(j)*OffsetFlip(j+1) <= 0
            Ref2 = j;
        end
    end
    Shift = Ref1 - Ref2;
    ShiftMat = zeros(NFiles, abs(Shift));
    
    FlipMat = fliplr(SliceDiff2);
    
    % Plot Difference spectra as a fn of offset enrgy
    figure(5);
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [5 5 17 9]);
    subplot (121);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [2 2 5 5]);
    contourf(Offset,T,SliceDiff1,...
        linspace(VMin,VMax,NContourLevels),'LineStyle','none');
    colormap('Jet');
    set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(gca,'YTick',[3 10 20 30 40 50],...
        'YTickLabel',['03'; '10'; '20'; '30'; '40'; '50']);
    xlabel('Offset energy (meV)');
    ylabel('Population time (ps)');
    title('High Absorption Energy Slice');
    c = colorbar('North');
%     set(c,'XTick',[0 0.03 0.06 0.09 1.20]);
    x=get(c,'Position');
    x(4)=0.5*x(4);
    x(2)=x(2)+0.25;
    set(c,'Position',x);
    subplot (122);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [10 2 5 5]);
    contourf(OffsetFlip,T,FlipMat,...
        linspace(VMin,VMax,NContourLevels),'LineStyle','none');
%     line([-SliceFreq2 -SliceFreq2], [T(1) T(end)],'Color',[0 0 0],'LineWidth',1);
    set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(gca,'YTick',[3 10 20 30 40 50],...
        'YTickLabel',['03'; '10'; '20'; '30'; '40'; '50']);
    xlabel('Emission energy (meV)');
    ylabel('Population time (ps)');
    title('Low Absorption Energy Slice');
    c = colorbar('North');
%     set(c,'XTick',[0 0.03 0.06 0.09 1.20]);
    x=get(c,'Position');
    x(4)=0.5*x(4);
    x(2)=x(2)+0.25;
    set(c,'Position',x);
    set(gca,'XDir','Reverse');
%     hold off;
    
    if Shift > 0
        FinFlipMat = cat(2,ShiftMat,FlipMat(:,1:NEmiDim - Shift));
    else
        FinFlipMat = cat(2,FlipMat(:,abs(Shift)+1:end),ShiftMat);
    end
    SDAsymm = SliceDiff1 - FinFlipMat;
    TotAsymm = sum(sum(abs(SDAsymm)));
%     dlmwrite(strcat(OutPath,'Asymmetry.dat'),SDAsymm,'\t');
    VMax = max(max(SDAsymm));
    VMin = min(min(SDAsymm));
    %%
    figure(4);
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [10 10 8 9]);
    set(gca, 'Units', 'centimeters');
    set(gca, 'Position', [2 2 5 5]);
    contourf(Offset,T,SDAsymm,...
        linspace(VMin,VMax,NContourLevels),'LineStyle','none');
        set(gca,'FontSize',10);
    set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
    set(gca,'YTick',[3 10 20 30 40 50],...
        'YTickLabel',['03'; '10'; '20'; '30'; '40'; '50']);
    xlabel('Energy offset (meV)');
    ylabel('Population time (ps)');
    title('Asymmetry between slices');
    c = colorbar('North');
    colormap('Jet');
%     set(c,'XTick',[0 0.03 0.06 0.09 1.20]);
    x=get(c,'Position');
    x(4)=0.5*x(4);
    x(2)=x(2)+0.25;
    set(c,'Position',x);
%     saveas(gcf, strcat(OutPath,'Asymmetry'), 'emf');
    
    %% Plot for poster
%     figure(5);
%     set(gcf, 'Units', 'centimeters');
%     set(gcf, 'Position', [2 2 29 16]);
%     subplot (121);
%     set(gca, 'Units', 'centimeters');
%     set(gca, 'Position', [2 2 15 13]);
%     contourf(gEmiFreq(1,:),T,SliceDiff1,...
%         linspace(VMin1,VMax1,NContourLevels),'LineStyle','none');
%     set(gca,'FontSize',20);
%     set(gca,'YTick',[3 10 20 30 40 50],...
%         'YTickLabel',[3 10 20 30 40 50]);
%     set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
%     xlabel('Emission energy (meV)','FontSize',24);
%     ylabel('Population time (ps)','FontSize',24);
% %     title('High Absorption Energy Slice', 'FontSize', 28);
%     colorbar('Location','North','FontSize',20,'TickLength',[0.015 0.015],...
%         'XTick',[-3 -2 -1 0 1]);
%     hold on;
%     
%     figure(5);
%     subplot (122);
%     set(gca, 'Units', 'centimeters');
%     set(gca, 'Position', [15.4 2 15 13]);
%     contourf(gEmiFreq(1,:),T,SliceDiff2,...
%         linspace(VMin2,VMax2,NContourLevels),'LineStyle','none');
%     set(gca,'FontSize',20);
%     set(gca,'YTick',[3 10 20 30 40 50],...
%         'YTickLabel',[ ]);
%     set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
% %     title('Low Absorption Energy Slice');
%     colorbar('Location','North','FontSize',20,'TickLength',[0.015 0.015],...
%         'XTick',[-2 -1 0 1]);
%     hold off;
% %%    
%     figure(6);
%     set(gcf, 'Units', 'centimeters');
%     set(gcf, 'position', [2 2 17 16]);
%     set(gca,'Units','centimeters');
%     set(gca,'Position', [3.2 2 15 13]);
%     contourf(Offset,T,SDAsymm,...
%         linspace(VMin,VMax,NContourLevels),'LineStyle','none');
%         set(gca,'FontSize',12);
%     set(gca,'YTick',[3 10 20 30 40 50],...
%         'YTickLabel',[3 10 20 30 40 50], 'FontSize',20);
%     set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
%     xlabel('Energy offset (meV)','FontSize',24);
%     ylabel('Population time (ps)','FontSize',24);
%     colorbar('Location','North','FontSize',20,'TickLength',[0.015 0.015]);

%     %% Figure for publications
%     fig10 = figure(10);
%     set(gcf, 'Units', 'centimeters');
%     set(gcf, 'Position', [5 5 6 5.2]);
%     % VMax = 1;
%     set(gca,'Units','centimeters');
%     set(gca,'Position', [1.6 1.5 5.1 3.4]);
%     hFig = contourf(gEmiFreq(1,:),T,SliceDiff,...
%         linspace(VMin,VMax, NContourLevels),...
%         'LineStyle', 'none');
% %     axis( AxisRange );
%     xlabel('   Emission\newlineenergy (meV)','FontSize',10);
%     ylabel('T (ps)','FontSize',10);
%     colormap jet;
% %     set(gca,'YTick',[-1546 -1544 -1542],'YTickLabel',[-1546 -1544 -1542],...
% %         'XTick',[1542 1544 1546],'XTickLabel',[1542 1544 1546]);
%     set(gca,'YTick',[10 30 50],'YTickLabel',[10 30 50]);
%     set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
%     colorbar('Location','North','TickLength',[0.025 0.025],...
%         'XTick',[-1.5 -0.5 0.5]);
end