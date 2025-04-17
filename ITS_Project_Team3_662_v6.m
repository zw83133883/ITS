%% ITS Project
% this script first fits candidate distributions to the speed data from T_speed_data 
% and creates histogram, runs ks tests, and displays qq plots.
% then it pauses and prompts for simulation parameters before calling the simulation function,
% which uses the best-fit distributions to simulate baseline vs predictive routing.
%
% part 1: distribution selection
clear; clc; close all;

% load the speed data 
load('Supporting_Data_Team_02.mat');  %this creates two primary datasets, t_speed data x4 and traffic mode

% build a little helper cell so you can keep using dataCell{i}
dataCell = cell(1,4);
for i = 1:4
  dataCell{i} = T_speed_data{:,i};
end

% interactively inspect each series in dfittool
dataNames = {'S65','S50','S40','S15'};
for i = 1:4
    fprintf('Launching dfittool for %s – close window when done\n', dataNames{i});
    dfittool(T_speed_data{:,i});
end

% now continue with scripted candidate‐distribution fitting...
candTypes = {'Normal','Lognormal','Gamma','Weibull'};
bestPD    = struct('name',cell(1,4),'pd',cell(1,4));
for i = 1:4
    d    = T_speed_data{:,i};

end

% candidate distributions to test
candTypes = {'Normal','Lognormal','Gamma','Weibull'};%can we re-code this the long way using DFITTOOL homework methods
dataNames = {'S65','S50','S40','S15'};
bestPD = struct('name',cell(1,4),'pd',cell(1,4));  % to store best fit for each column

% candidate distributions to test
candTypes = {'Normal','Lognormal','Gamma','Weibull'};
dataNames = {'S65','S50','S40','S15'};
bestPD    = struct('name',cell(1,4),'pd',cell(1,4));

%all S## “Fits” in one 2×2 subplot
figure('Name','All Speed Fits','NumberTitle','off');
for i = 1:4
    d       = dataCell{i};                    %original d = T_speed_data{:,i}
    pVals   = zeros(1,numel(candTypes));     % original pVals
    pdFits  = cell(1,numel(candTypes));      % original pdFits
    x_rng   = linspace(min(d),max(d),200);
    cols    = lines(numel(candTypes));

    subplot(2,2,i);
    histogram(d,'Normalization','pdf'); hold on;
    for j = 1:numel(candTypes)
        pdFits{j} = fitdist(d,candTypes{j});           % original fit
        [~,pVals(j)] = kstest(d,'CDF',pdFits{j});      % original KS test
        plot(x_rng, pdf(pdFits{j},x_rng), ...
             'Color',cols(j,:),'LineWidth',1.5);
    end
    title(['Fits for ', dataNames{i}]);
    xlabel('speed (mph)'); ylabel('pdf');
    legend(candTypes,'Location','best','FontSize',6);
    hold off;

    % store best
    [~,bestIdx]      = max(pVals);
    bestPD(i).name  = candTypes{bestIdx};
    bestPD(i).pd    = pdFits{bestIdx};
end

%all QQplots in one 2×2 subplot
figure('Name','All QQ Plots','NumberTitle','off');
for i = 1:4
    d    = dataCell{i};               % same as before
    nPts = numel(d);

    subplot(2,2,i);
    qqplot(d, random(bestPD(i).pd,nPts,1)); 
    title(['qq plot for ', dataNames{i}, ' best fit: ', bestPD(i).name]);
end

    switch bestPD(i).name %is this part nesescary, can instead we just list all the results out in the command window and provide selection criteria
        case 'Normal'
            fprintf('  estimated mu = %.4f, sigma = %.4f\n', bestPD(i).pd.mu, bestPD(i).pd.sigma);
        case 'Lognormal'
            fprintf('  estimated mu = %.4f, sigma = %.4f\n', bestPD(i).pd.mu, bestPD(i).pd.sigma);
        case 'Gamma'
            fprintf('  estimated a (shape) = %.4f, b (scale) = %.4f\n', bestPD(i).pd.a, bestPD(i).pd.b);
        case 'Weibull'
            fprintf('  estimated A (scale) = %.4f, B (shape) = %.4f\n', bestPD(i).pd.A, bestPD(i).pd.B);

    end
    % % show qq plot for best candidate
    % figure;
    % qqplot(d, random(bestPD(i).pd, nPts, 1));
    % title(['qq plot for ', dataNames{i}, ' best fit: ', bestPD(i).name]);

%% part 2: run sim with user inputs
disp(' ');
disp('entering simulation phase baby!');
nT    = input('enter number of trips to simulate: ');
confLev = input('enter confidence level (suggested .90-.99999): ');
mu    = input('enter lognormal mu for trip distances (default to 1): ');%need code to display these distribution perameters in section 1
sigma = input('enter lognormal sigma for trip distances (default to 1.6): ');%need code to display these perameters in section 1
uniq  = input('enter oversample parameter for genlogntrips (suggested 1-20): ');% used to control how many candidate trips are generated before selecting a subset that matches the desired lognormal distribution %default to 20, but use less when testing the code

simulateITS(nT, confLev, mu, sigma, uniq, bestPD, G, T_roadcond_data);

%% simulate ITS
function simulateITS(nT, confLev, mu, sigma, uniq, bestPD, G, T_roadcond_data)
    % this function runs the routing simulation using the best fit speed distributions
    % bestPD is a structure array with fields:
    %   bestPD(1): best fit for S65 (nominal interstate speeds)
    %   bestPD(2): best fit for S50 (nominal highway speeds)
    %   bestPD(3): best fit for S40 (construction conditions)
    %   bestPD(4): best fit for S15 (accident conditions)
    
    mpH = 60;  % minutes per hour
    
    %% parallel pool
    delete(gcp('nocreate')) %clear out the pool in case theres soemthing idle in t here
    parpool('local',6);%change this to whatever your computer and license will allow
    
    %% load data
    load('EastCoast.mat');      % lvariable G
    load('Supporting_Data_Team_02.mat'); % for road condition data
    
    %% generate trips 
    trips = genlogntrips(G, nT, confLev, mu, sigma, uniq);
    
    %% nominal routing
    Gb = G;  % baseline graph
    numE = height(Gb.Edges);
    spdBase = zeros(numE,1);

    for e = 1:numE
        if Gb.Edges.Speed(e) == 65
            spdBase(e) = random(bestPD(1).pd);
        else
            spdBase(e) = random(bestPD(2).pd);
        end
    end
    Gb.Edges.Weight = Gb.Edges.Distance ./ spdBase * mpH;
    
    tBase = zeros(nT,1);
    tPred = zeros(nT,1);
    rerouted = false(nT,1);
    
    %% predictive routing 
    allCond = [];
    for col = 1:width(T_roadcond_data)
        allCond = [allCond; T_roadcond_data{:,col}];
    end
    condStr = string(allCond);
    pNorm = sum(condStr=="normal") / numel(condStr);
    pAcc  = sum(condStr=="accident") / numel(condStr);
    pCons = sum(condStr=="construction") / numel(condStr);
    if pNorm==0 && pAcc==0 && pCons==0
        pNorm = 0.90; pAcc = 0.05; pCons = 0.05;
    end
    cProbs = [pNorm, pAcc, pCons];

    % normal- if edge speed is 65 -> bestPD(1), if 50 -> bestPD(2)
    % accident- use bestPD(4) (S15)
    % construction- use bestPD(3) (S40)
    
    %% actual simulation loop
    for i = 1:nT
        sn = trips(i,1);  
        en = trips(i,2);  
        [pB, tB] = shortestpath(Gb, sn, en, 'Method','positive');
        if isempty(pB)
            tBase(i) = NaN; tPred(i) = NaN;
            continue;
        end
        tBase(i) = tB;
        
        Gp = G;
        ne = height(Gp.Edges);
        rv = rand(ne,1);
        cp = cumsum(cProbs);
        spdEff = zeros(ne,1);
        for e = 1:ne
            % find nominal distribution based on edge speed
            if Gp.Edges.Speed(e)==65
                nomPD = bestPD(1).pd;
            else
                nomPD = bestPD(2).pd;
            end
            if rv(e) <= cp(1)
                % normal- use nominal distribution
                spdEff(e) = random(nomPD);
            elseif rv(e) <= cp(2)
                % accident use S15 distribution
                spdEff(e) = random(bestPD(4).pd);
            else
                % construction- use S40 distribution
                spdEff(e) = random(bestPD(3).pd);
            end
        end
        Gp.Edges.Weight = Gp.Edges.Distance ./ spdEff * mpH;
        
        [pP, tP] = shortestpath(Gp, sn, en, 'Method','positive');%from word doc
        if isempty(pP)
            tPred(i) = NaN;
        else
            tPred(i) = tP;
        end
        
        if length(pP) ~= length(pB) || any(pP ~= pB)
            rerouted(i) = true;
        end
    end
    
    %% complie results
    valid = ~isnan(tBase) & ~isnan(tPred);
    tBase = tBase(valid);
    tPred = tPred(valid);
    rerouted = rerouted(valid);
    nValid = length(tBase);
    
    savings = (tBase - tPred) ./ tBase * 100;
    mSave = mean(savings);
    sdSave = std(savings);
    alpha = 1 - confLev;
    tval = tinv(1 - alpha/2, nValid-1);
    err = tval * sdSave / sqrt(nValid);
    ciLo = mSave - err;
    ciHi = mSave + err;
    
    %display results
    fprintf('\npredictive routing savings:\n');
    fprintf('mean: %.2f%%, 95%% ci: [%.2f%%, %.2f%%]\n', mSave, ciLo, ciHi);
    if ciLo > 5
        fprintf('=> >5%% saving (95%% conf)\n');
    else
        fprintf('=> not >5%% saving (95%% conf)\n');
    end
    
    fracRer = sum(rerouted)/nValid;
    fprintf('fraction rerouted: %.2f%%\n', fracRer*100);
    
    figure;
    histogram(savings);
    xlabel('time savings (%)'); ylabel('count');
    title('predictive routing time savings distribution');
end
