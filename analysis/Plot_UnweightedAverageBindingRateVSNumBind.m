%% Analysis_TransitionMatrix

%% lclemens@uci.edu
clear all;
close all;

%% Set parameters for model

% model name
modelName = 'LocalStructuring';

% find files
folder    = '../data';
filesubfolder = '';
filetitle     = 'TCRZeta.Stiffenrange';

% find file, define states for phosphostates plot
phosphostate_filename = 'OccupiediSitesZeta.txt';
phosphostates = 1:1:64;
phosphosites = 1:1:6;

% parameter for TCR zeta
locationTotal = 6;
NFilSweep = 1;
iSiteTotal(1:NFilSweep) = [6];

% local stiffening parameter
sweep = 5;
sweepParameter = 'StiffenRange';

% figure parameters
legendlabels = {'No Stiffening', 'StiffenRange = 0','StiffenRange = 1','StiffenRange = 2',...
                'StiffenRange = 3','StiffenRange = 4','StiffenRange = 5'};
colorIndices = sweep+2;
colors = flipud(cool(max(sweep)+2));
ms = 10;
lw = 2;
modificationLabel = '(Phosphorylated)';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Binding rates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables

sumRates = zeros(length(sweep),locationTotal,2);
avgRates = zeros(length(sweep),locationTotal,2);

%% Create transition matrices, calculate average binding rates
for nfSweep = 1:length(NFilSweep)

    % set NFil, ITAM locations for this iteration
    NFil = NFilSweep(nfSweep);

    % start parameter sweep
    for s = 1:length(sweep)

        % choose file
        filename = strcat(filetitle,'.',num2str(sweep(s)),'.cat');
        disp(filename);

        % initialize
        OccupiedLocations = zeros(2^locationTotal,1);
        OccupiedLocationsMatrix = zeros(2^locationTotal,locationTotal);
        OccupiedLocationsDecimal = zeros(2^locationTotal,1);
        POcc = zeros(2^locationTotal,locationTotal);
        PBind = zeros(2^locationTotal,locationTotal);
        PBindKinase = zeros(2^locationTotal,locationTotal);

        %% Read from File
        M = dlmread(fullfile(folder,filesubfolder, filename));
        ntMetropolis = M(:,1);

        OccupiedLocations = M(:,end);
        OccupiedLocationsMatrix(:,1:locationTotal) = M(:,(end-locationTotal):(end-1));

        % up to total number of iSites - 6 for TCR Zeta
        siteCounter = 1;

        % starting index - 8+2*(locationTotal+1) is output only once, 6+2 takes
        % us to the correct index in the filament output
        ind = 8+2*(locationTotal+1)+6+2;
        for nf = 1:NFil
            if(nf>1)
                ind = ind + (6 + 7*iSiteTotal(nf-1) + 2 + NFil + NFil);
            end

            for iy = 1:iSiteTotal(nf)
                POcc(:,siteCounter) = M(:,ind + 7*(iy-1));
                PBind(:,siteCounter) = M(:,ind + 7*(iy-1) + 1);
                siteCounter = siteCounter + 1;
            end
        end

        %% Convert binary to decimal
        for j=1:2^locationTotal
            binaryString = num2str(OccupiedLocations(j));
            OccupiedLocationsDecimal(j) = bin2dec(binaryString);
        end

        %% Create kinase transitionMatrix (forward binding) and Phosphatase transitionMatrix (reverse binding)
        PBindKinase = PBind.*(~OccupiedLocationsMatrix);
        PBindPhosphatase = PBind.*(OccupiedLocationsMatrix);
        for j=1:2^locationTotal
            PBindKCorrectOrder(OccupiedLocationsDecimal(j)+1,:) = PBindKinase(j,:);
            PBindPCorrectOrder(OccupiedLocationsDecimal(j)+1,:) = PBindPhosphatase(j,:);
        end
        KinaseTransition = fliplr(PBindKCorrectOrder);
        PhosphataseTransition = fliplr(PBindPCorrectOrder);

        %% Calculate total possible ways to transition from State i to State i+1
        for i=1:locationTotal
            totalRates(i) = nchoosek(locationTotal,i-1).*(locationTotal-(i-1));
        end

        %% Calculate average rates of transition from one state to next

        for j=1:2^locationTotal % for each state (i.e. 010010)
            totalOccupied(1) = size(find(KinaseTransition(j,:)==0),2); % find number of 0 entries (aka phosphorylated sites)
            totalOccupied(2) = size(find(PhosphataseTransition(j,:)==0),2); % find number of 0 entries (aka unphosphorylated sites)
            if(totalOccupied(1)<locationTotal) % if not completely occupied
                sumRates(s,totalOccupied(1)+1,1) = sumRates(s,totalOccupied(1)+1,1)+sum(KinaseTransition(j,:));
            end
            if(totalOccupied(2)<locationTotal)
                sumRates(s,totalOccupied(2)+1,2) = sumRates(s,totalOccupied(2)+1,2)+sum(PhosphataseTransition(j,:));
            end
        end

        % find average rates of transition from one state to another
        avgRates(s,:,1) = sumRates(s,:,1)./totalRates;
        avgRates(s,:,2) = sumRates(s,:,2)./totalRates;


    end
end

% state matrix - 0 or 1;
phosphoStates = dlmread(fullfile(folder,phosphostate_filename),'_');

% for plotting binding rates of last set
PBindPlot = (1-POcc).*(1-phosphoStates);
PBindPlot(PBindPlot==0) = NaN;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot phosphostate combinations

% Vertical - phosphostates on yaxis
figure('Name','Phosphostates'); clf;
hm = heatmap(phosphosites,phosphostates,phosphoStates);
hm.Colormap = flipud(gray);
hm.ColorbarVisible = 'off';
for l=1:64
    ydispLabels{l} = {''};
end
for l=1:6
    xdispLabels{l} = {''};
end
hm.CellLabelColor = 'None';
hm.YDisplayLabels = ydispLabels;
hm.XDisplayLabels = xdispLabels;
set(gcf,'units','inches','Position',[0.5 0.5 0.85 8.75]);


%% Plot binding rates of last sweep parameter
figure('Name','Binding rates'); clf;
hm = heatmap(phosphosites,phosphostates,PBindPlot);
hm.Colormap = pink;
hm.MissingDataLabel = 'Phosphorylated';

hm.ColorLimits = [0 ceil(max(PBindPlot(PBindPlot > 0))*100)/100];
for l=1:64
    ydispLabels{l} = {''};
end

hm.YDisplayLabels = ydispLabels;
hm.CellLabelFormat = '%0.4f';
hm.FontSize = 12;
hm.ColorbarVisible = 'on';

for l=1:6
    xdispLabels{l} = {''};
end
hm.XDisplayLabels = xdispLabels;
xlabel('Binding site');
ylabel('State');
set(gcf,'units','inches','Position',[1 1 6.5 11]);


%% Plot average binding rate unweighted - Phosphorylation
figure('Name','Average Binding Rate - Phosphorylation'); clf; hold on; box on;

% plot
for s = 1:length(sweep)
    plotColor = colorIndices(s);
    plotData = reshape(avgRates(s,:,1),[1 locationTotal]);
    plot(1:1:locationTotal,plotData,'-s','LineWidth',lw,'Color',colors(plotColor,:),'MarkerFaceColor',colors(plotColor,:),'MarkerSize',ms);
end

% set axis labels and scales
set(gca,'XTick',1:1:locationTotal);
set(gca,'XTickLabel',{'0 -> 1', '1 -> 2', '2 -> 3', '3 -> 4','4 -> 5', '5 -> 6', '6 -> 7', '7 -> 8', '8 -> 9', '9 -> 10'});
xlabel1 = {['Number of Modified Sites'],modificationLabel};
ylabel1 = {['Average Binding Rate'],['(per free space binding)']};
title1 = 'Average Binding Rate vs Total Modified Sites';

% set axes limits
ylim([0.01,0.035]); % TCRZeta Membrane On

% set second position and show labels
pos = get(gcf, 'position');
set(gcf,'units','centimeters','position',[[1,1],30,25]);
set(gca,'FontName','Arial','FontSize',18);
xlabel(xlabel1,'FontName','Arial','FontSize',18);
ylabel(ylabel1,'FontName','Arial','FontSize',18);
title(title1,'FontName','Arial','FontSize',18);

% set colorbar parameters based on model
set(gcf,'Colormap',cool)
colormap cool;
h = colorbar;
h = colorbar('Ticks',[0 1],'TickLabels',{'',''},'YDir','reverse');
set(h,'ylim',[0 1]);

%% Plot average binding rate unweighted - Dephosphorylation
figure('Name','Average Binding Rate - Dephosphorylation'); clf; hold on; box on;

% plot
for s = 1:length(sweep)
    plotColor = colorIndices(s);
    plotData = reshape(avgRates(s,:,2),[1 locationTotal]);
    plot(1:1:locationTotal,plotData,'-s','LineWidth',lw,'Color',colors(plotColor,:),'MarkerFaceColor',colors(plotColor,:));
end

% set axis labels and scale
set(gca,'XTick',1:1:locationTotal);
set(gca,'XTickLabel',{'6 -> 5','5 -> 4','4 -> 3','3 -> 2','2 -> 1','1 -> 0'});
title1 = 'Average Binding Rate vs Total Modified Sites';
xlabel1 = {['Number of Modified Sites'],modificationLabel};
ylabel1 = {['Average Binding Rate'],['(per free space binding)']};

% set second position and show labels
pos = get(gcf, 'position');
set(gcf,'units','centimeters','position',[[1,1],30,25]);
set(gca,'FontName','Arial','FontSize',18);
xlabel(xlabel1,'FontName','Arial','FontSize',18);
ylabel(ylabel1,'FontName','Arial','FontSize',18);
title(title1,'FontName','Arial','FontSize',18);

% set axes limits
ylim([0,0.6]); % TCRZeta Membrane On

% set colorbar parameters
set(gcf,'Colormap',cool)
colormap cool;
h = colorbar;
h = colorbar('Ticks',[0 1],'TickLabels',{'',''},'YDir','reverse');
set(h,'ylim',[0 1]);
