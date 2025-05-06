%% ITS Project
%team 2
%spring 2025
%SYST 662
% this script first fits candidate distributions to the speed data from T_speed_data 
% and creates histogram, runs ks tests, and displays qq plots.
% then it pauses and prompts for simulation parameters before calling the simulation function,
% which uses the best-fit distributions to simulate baseline vs predictive routing.
%
% part 1: distribution selection
clear; clc; close all;

% load the data sets
% this creates two primary datasets, t_speed data x4 and traffic mode
load('Supporting_Data_Team_02.mat');  
load('EastCoast.mat');

dataCell = cell(1,4);
for i = 1:4
    % format the data in a simpler way for later
  dataCell{i} = T_speed_data{:,i};
end

% interactively inspect each series in dfittool un-comment to activate the
% pop-up
dataNames = {'S65','S50','S40','S15'};
% % for i = 1:4
% %     fprintf('Launching dfittool for %s – close window when done\n', dataNames{i});
% %     dfittool(T_speed_data{:,i});
% % end

candTypes = {'Normal','Lognormal','Gamma','Weibull'};
bestPD = struct('name',cell(1,4),'pd',cell(1,4));

%all S## “Fits” in one 2×2 subplot
figure('Name','All Speed Fits','NumberTitle','off');
for i = 1:4
    d       = dataCell{i};                   % original d = T_speed_data{:,i}
    pVals   = zeros(1,numel(candTypes));     % original pVals
    pdFits  = cell(1,numel(candTypes));      % original pdFits
    x_rng   = linspace(min(d),max(d),200);
    cols    = lines(numel(candTypes));

    subplot(2,2,i);
    % increased bin size to 42 to make it snappier
    thebigoleoneoffourhistograms= histogram(d, 42, 'Normalization','pdf');
    % handle visibility should correct our legend issue
    set(thebigoleoneoffourhistograms,'HandleVisibility','off'); 
    hold on; 
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

    % Store best
    [~,bestIdx]      = max(pVals);
    bestPD(i).name  = candTypes{bestIdx};
    bestPD(i).pd    = pdFits{bestIdx};
    
    % Display all candidate functions parameters along with KS p‐value
fprintf('\nParameters for %s:\n', dataNames{i});
for j = 1:numel(candTypes)
    pdj = pdFits{j};
    switch candTypes{j}
      case 'Normal'
        fprintf('  %-8s  mu = %.3f, sigma = %.3f,  KS p = %.3f\n', ...
                candTypes{j}, pdj.mu, pdj.sigma, pVals(j));
      case 'Lognormal'
        fprintf('  %-8s  mu = %.3f, sigma = %.3f,  KS p = %.3f\n', ...
                candTypes{j}, pdj.mu, pdj.sigma, pVals(j));
      case 'Gamma'
        fprintf('  %-8s  a = %.3f, b = %.3f,       KS p = %.3f\n', ...
                candTypes{j}, pdj.a, pdj.b, pVals(j));
      case 'Weibull'
        fprintf('  %-8s  A = %.3f, B = %.3f,       KS p = %.3f\n', ...
                candTypes{j}, pdj.A, pdj.B, pVals(j));
    end
end
fprintf('  → selected: %s\n', bestPD(i).name);
end

%all QQplots in one 2×2 subplot
figure('Name','All QQ Plots','NumberTitle','off');
for i = 1:4
    d    = dataCell{i};               % same as before
    nPts = numel(d);

    subplot(2,2,i);
    qqplot(d, random(bestPD(i).pd,nPts,1)); 
    title(['qq plot for ', dataNames{i}, ' best fit: ', bestPD(i).name]);

    % is this part nesescary, can instead 
    % we just list all the results out in the command window 
    % and provide selection criteria
    switch bestPD(i).name 
        case 'Normal'
            fprintf('  estimated mu = %.4f, sigma = %.4f\n', ...
                bestPD(i).pd.mu, bestPD(i).pd.sigma);
        case 'Lognormal'
            fprintf('  estimated mu = %.4f, sigma = %.4f\n', ...
                bestPD(i).pd.mu, bestPD(i).pd.sigma);
        case 'Gamma'
            fprintf('  estimated a (shape) = %.4f, b (scale) = %.4f\n', ...
                bestPD(i).pd.a, bestPD(i).pd.b);
        case 'Weibull'
            fprintf('  estimated A (scale) = %.4f, B (shape) = %.4f\n', ...
                bestPD(i).pd.A, bestPD(i).pd.B);

    end
    % % show qq plot for best candidate
    % figure;
    % qqplot(d, random(bestPD(i).pd, nPts, 1));
    % title(['qq plot for ', dataNames{i}, ' best fit: ', bestPD(i).name]);
end
%% Part 2: Run simulation with user inputs
disp(' ');
disp('Entering simulation phase! Press Ctrl+C to abort at any time.');

% Prompt for number of trips
nT = input('Enter number of trips to simulate [1–9999, default = 1000]: ');
if isempty(nT), nT = 1000; end
if nT < 1 || nT > 9999
    error('Number of trips must be between 1 and 9999.');
end

% Prompt for confidence level
confLev = input('Enter confidence level (0.90–0.99999, default = 0.95): ');
if isempty(confLev), confLev = 0.95; end
if confLev < 0.90 || confLev > 0.99999
    error('Confidence level must be between 0.90 and 0.99999.');
end

% Prompt for lognormal mu
mu = input('Enter lognormal mu for trip distances (default = 1): ');
if isempty(mu), mu = 1; end

% Prompt for lognormal sigma
sigma = input('Enter lognormal sigma for trip distances (default = 1.6): ');
if isempty(sigma), sigma = 1.6; end

% Prompt for oversampling parameter
uniq = input('Enter oversample parameter for genlogntrips [1–20, default = 1]: ');
if isempty(uniq), uniq = 1; end
if uniq < 1 || uniq > 20
    error('Oversample parameter must be between 1 and 20.');
end

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
    
    
    %% generate trips 
    trips = genlogntrips(G, nT, confLev, mu, sigma, uniq);
    
    %Verify trip length is 10 mi on avg, with sigma 33miles (requirement 11)

        d = trips(:,3)/100;  % convert back from hundredths of a mile
    fprintf('Empirical trip distances: mean = %.2f miles, std = %.2f miles\n', ...
            mean(d), std(d));

    %to prove we meet the mu and sigma target
    figure;
    histogram(d, 20,'Normalization','pdf');
    hold on;
    x = linspace(min(d),max(d),200);
    theoretical = pdf( makedist('Lognormal','mu',1,'sigma',1.6), x );
    plot(x, theoretical, 'r-','LineWidth',1.5);
    hold off;
    xlabel('Distance (miles)');
    ylabel('PDF');
    title('Trip‐lengths: empirical vs. LogNormal(\mu=1,\sigma=1.6)');
    legend('Empirical','Theoretical');
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
    
    %Diagnostic output baseline speed draws
    figure('Name','Diagnostic: Speed Draws','NumberTitle','off');
    subplot(2,1,1);
    histogram(spdBase,'Normalization','pdf');
    xlabel('speed (mph)');
    ylabel('pdf');
    title('Baseline speed draws (spdBase)');

    % preallocate for collecting every predictive-speed draw
    % clear it out in case we keep re-running and muck up the memory
    spdEffAll = zeros(numE, nT);
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

        spdEffAll(:,i) = spdEff;%collect diagnostic data
        Gp.Edges.Weight = Gp.Edges.Distance./spdEff*mpH;
        
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
     %Diagnostic- show predictive-speed draws
    subplot(2,1,2);
    histogram(spdEffAll(:),'Normalization','pdf');
    xlabel('speed (mph)');
    ylabel('pdf');
    title('Predictive speed draws (all spdEff)');
    %% compile results
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
    histogram(savings, 10);%change bin size here as needed
    xlabel('time savings (%)'); ylabel('count');
    title('predictive routing time savings distribution');
end

%% Part 4: validation using Monte Carlo replications
nVal = input('Enter # validation replications (default 100): ');
if isempty(nVal), nVal = 100; end
validateITS(nT, confLev, mu, sigma, uniq, bestPD, G, T_roadcond_data, nVal)

function validateITS(nT, confLev, mu, sigma, uniq, bestPD, G, Tcond, nVal)
    
    % Monte Carlo validation: replications of mean %saving
    valSave = zeros(nVal,1);
    rng('default'); 
    for k=1:nVal
        rng(k);  %different seed each rep
        [tB,tP] = runOne(nT, confLev, mu, sigma, uniq, bestPD, G, Tcond);
        valSave(k) = mean((tB-tP)./tB*100);
    end
    valSave = valSave 
    mVS = mean(valSave); sVS = std(valSave);
    alpha = 1 - confLev;
    tcrit = tinv(1-alpha/2,nVal-1);
    err = tcrit * sVS / sqrt(nVal);
    ci = mVS + [-1,1]*err;

    fprintf(['\nValidation (%d reps): mean saving=%.2f%%, ' ...
        'CI=[%.2f%%,%.2f%%]\n'], nVal, mVS, ci);  

    figure; histogram(valSave, 10, 'Normalization','pdf');
    hold on
    pdN = fitdist(valSave,'Normal');
    xs = linspace(min(valSave),max(valSave),200);
    plot(xs,pdf(pdN,xs),'LineWidth',1.5);
    title('Validation savings dist vs Normal fit');
    legend('Empirical','Normal fit'); hold off

    figure; qqplot(valSave);
    title('QQ plot of validation %savings');

    %diag plot of cumulative mean to see stability of montecarlo section
    %code
    cumMean = cumsum(valSave) ./ (1:nVal)';
    figure;
    plot(1:nVal, cumMean, 'LineWidth',1.5);
    xlabel('Replication #');
    ylabel('Cumulative mean % saving');
    title('Convergence of Monte Carlo estimate');
    yline(cumMean(end),'r--','Final mean');

%plot CI for the montecarlo for diagnostic of the montecarlo code by
%sensitivity
    ciWidth = 2 * tinv(1 - (1-confLev)/2, (1:nVal)'-1) .* ...
          arrayfun(@(k) std(valSave(1:k)),1:nVal) ./ sqrt(1:nVal)';
    figure;
    plot(1:nVal, ciWidth);
    xlabel('Replications');
    ylabel('CI Width (% saving)');
    title('How CI narrows with more reps');

end

function [tBase, tPred] = runOne(nT, confLev, mu, sigma, uniq, bestPD, G, Tcond)
    mpH = 60;

    %generate the same trips as in simulateITS
    trips = genlogntrips(G, nT, confLev, mu, sigma, uniq);

    %probabilities
    allC = [];
    for col = 1:width(Tcond)
        allC = [allC; Tcond{:,col}];
    end
    pNorm = sum(allC=="normal")/ numel(allC);
    pAcc  = sum(allC=="accident")/ numel(allC);
    pCons = sum(allC=="construction")/ numel(allC);
    if pNorm + pAcc + pCons == 0
        pNorm = 0.90; pAcc = 0.05; pCons = 0.05;
    end
    cp = cumsum([pNorm, pAcc, pCons]);

    E = height(G.Edges);
    %draw one - mixture of speeds for every edge
    rv     = rand(E,1);
    spdEff = zeros(E,1);
    for e = 1:E
        if rv(e) <= cp(1)
            %nominal
            if G.Edges.Speed(e) == 65
                pd = bestPD(1).pd;
            else
                pd = bestPD(2).pd;
            end
        elseif rv(e) <= cp(2)
            %accident
            pd = bestPD(4).pd;
        else
            %construction
            pd = bestPD(3).pd;
        end
        spdEff(e) = random(pd);
    end

    %baseline route and time under spdEff
    %first, we pick the baseline path by nominal speeds
    Gb = G;
    spdNom = zeros(E,1);
    for e = 1:E
        if G.Edges.Speed(e) == 65
            nomPD = bestPD(1).pd;
        else
            nomPD = bestPD(2).pd;
        end
        spdNom(e) = mean(nomPD); %average of the nominal fit
    end
    Gb.Edges.Weight = Gb.Edges.Distance ./ spdNom * mpH;
    
    %compute tBase per trip
    tBase = nan(nT,1);
    for i = 1:nT
        sn = trips(i,1);
        en = trips(i,2);
        [pB, ~] = shortestpath(Gb, sn, en, 'Method','positive');
        if isempty(pB)
            tBase(i) = NaN;
        else
            tBase(i) = sum(G.Edges.Distance(pB)./spdEff(pB)) * mpH;
        end
    end
    
    %predictive route and time under the spdEff
    Gp = G;
    Gp.Edges.Weight = G.Edges.Distance ./ spdEff * mpH;
    
    % compute tPred per trip
    tPred = nan(nT,1);
    for i = 1:nT
        sn = trips(i,1);
        en = trips(i,2);
        [pP, ~] = shortestpath(Gp, sn, en, 'Method','positive');
        if isempty(pP)
            tPred(i) = NaN;
        else
            tPred(i) = sum(G.Edges.Distance(pP)./spdEff(pP)) * mpH;
        end
    end

    %drop NaNs
    ok = ~isnan(tBase) & ~isnan(tPred);
    tBase  = tBase(ok);
    tPred  = tPred(ok);
end
