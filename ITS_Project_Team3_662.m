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

dataNames = {'S65','S50','S40','S15'}; %titles
candTypes = {'Normal','Lognormal','Gamma','Weibull'}; %data under comparison
dataCell  = cell(1,4); %holds the speed data
for i=1:4
    dataCell{i} = T_speed_data{:,i};
end

%store the KS‐best and LL‐best fits
bestPD_KS = struct('name',cell(1,4),'pd',cell(1,4));
bestPD_LL = struct('name',cell(1,4),'pd',cell(1,4));

% plot all four histograms + fits in one figure
figure('Name','All Speed Fits','NumberTitle','off');
for i = 1:4
    d = dataCell{i};
    x_rng = linspace(min(d),max(d),200); %x values for plot
    cols = lines(numel(candTypes)); %change colors
    pVals = zeros(1,numel(candTypes)); %hold values
    logLikes = zeros(1,numel(candTypes));%hold values
    pdFits = cell(1,numel(candTypes));  %hold objects
    
    subplot(2,2,i);
    h = histogram(d,42,'Normalization','pdf'); hold on;
    set(h,'HandleVisibility','off'); %correct display issue
    for j=1:numel(candTypes)
        pdFits{j} = fitdist(d,candTypes{j}); %fit data
        [~,pVals(j)] = kstest(d,'CDF',pdFits{j}); %ks test
        logLikes(j) = sum(log(pdf(pdFits{j},d))); %lkog l
        plot(x_rng, pdf(pdFits{j},x_rng), 'Color',cols(j,:),'LineWidth',1.5); % overlay pdf fit

    end

    %Formatting
    hold off;
    title(['Fits for ',dataNames{i}]);
    xlabel('speed (mph)'); ylabel('pdf');
    legend(candTypes,'Location','best','FontSize',6);
    
    % pick the KS‐best and LL‐best
    [~,idxKS] = max(pVals); %index of max KS p values
    [~,idxLL] = max(logLikes); %index of max logL values
    bestPD_KS(i) = struct('name',candTypes{idxKS}, 'pd',pdFits{idxKS});
    bestPD_LL(i) = struct('name',candTypes{idxLL}, 'pd',pdFits{idxLL});
    
    % print out
    fprintf('\n%s: KS→%-8s (p=%.3f),  LL→%-8s (logL=%.1f)\n', ...
            dataNames{i}, ...
            candTypes{idxKS}, pVals(idxKS), ...
            candTypes{idxLL}, logLikes(idxLL));
    for j=1:numel(candTypes)
        pdj = pdFits{j};
        switch candTypes{j}
          case {'Normal','Lognormal'}
            fprintf('   %-8s mu=%.3f σ=%.3f  p=%.3f  logL=%.1f\n', ...
                    candTypes{j}, pdj.mu, pdj.sigma, pVals(j), logLikes(j));
          case 'Gamma'
            fprintf('   %-8s a=%.3f b=%.3f  p=%.3f  logL=%.1f\n', ...
                    candTypes{j}, pdj.a, pdj.b, pVals(j), logLikes(j));
          case 'Weibull'
            fprintf('   %-8s A=%.3f B=%.3f  p=%.3f  logL=%.1f\n', ...
                    candTypes{j}, pdj.A, pdj.B, pVals(j), logLikes(j));
        end
    end
end

%plot QQ‐plots X2
figure('Name','QQ Plots (KS winners)','NumberTitle','off');
for i=1:4
    subplot(2,2,i);
    qqplot(dataCell{i}, random(bestPD_KS(i).pd,numel(dataCell{i}),1));
    title([dataNames{i} ' KS winner: ' bestPD_KS(i).name]);
end

figure('Name','QQ Plots (LL winners)','NumberTitle','off');
for i=1:4
    subplot(2,2,i);
    qqplot(dataCell{i}, random(bestPD_LL(i).pd,numel(dataCell{i}),1));
    title([dataNames{i} ' LL winner: ' bestPD_LL(i).name]);
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

%Prompt for fit choice ks or LogL
fitChoice = input('Use KS fits or LL fits? Enter ''KS'' (default) or ''LL'' (Fred''s Favorite): ', 's' );
if isempty(fitChoice) || strcmpi(fitChoice,'KS')
    bestPD = bestPD_KS;
else strcmpi(fitChoice,'LL');
    bestPD = bestPD_LL;
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
uniq = input('Enter oversample parameter for genlogntrips [1–20, default = 20]: ');
if isempty(uniq), uniq = 20; end % default 20 ~87% unique
if uniq < 1 || uniq > 20
    error('Oversample parameter must be between 1 and 20.');
end

nVal = input('Enter # validation replications [1–10000, default = 100]: ');
if isempty(nVal), nVal = 100; end
if nVal < 1 || nVal > 10000
    error('# of validation replications must be between 1 and 10000.');
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
    
    %% parallel pool
    %%%delete(gcp('nocreate')) %clear out the pool in case theres soemthing idle in t here
    %%%parpool('local',6);%change this to whatever your computer and license will allow
    c = parcluster('local');
    c.NumWorkers = 6; %Increase its maximum workers to n     
    c.saveProfile; %lock it in    
    %Restart any existing pool for reliability 
    %start a new one with n workers
    delete(gcp('nocreate'));
    parpool(c, 6);
    
    %% generate trips 
    trips = genlogntrips(G, nT, confLev, mu, sigma, uniq);
    
    % Verify trip length is 10 mi on avg, with sigma 33miles (requirement 11)

    d = trips(:,3)/100;  % convert back from hundredths of a mile
    fprintf('Empirical trip distances: mean = %.2f miles, std = %.2f miles\n', ...
            mean(d), std(d));

    % to prove we meet the mu and sigma target
    figure;
    histogram(d, 45,'Normalization','pdf');
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
    mpH = 60;  % minutes per hour
    Gb = G;  % baseline graph
    numE = height(Gb.Edges);
    % Gb.Edges.Weight = Gb.Edges.Distance ./ Gb.Edges.Speed * mpH; % weight by time
    
    %diagnostic basline speed spread
    spdBaseAll = zeros(numE,nT);
    for run=1:nT
      for e=1:numE
        % pick the correct PD for edge e (65 or 50)
        if    G.Edges.Speed(e)==65, nomPD = bestPD(1).pd;
        else                      nomPD = bestPD(2).pd;
        end
        spdBaseAll(e,run) = random(nomPD);
      end
    end
    
    % preallocate for collecting every predictive-speed draw
    % clear it out in case we keep re-running and muck up the memory
    spdEffAll = zeros(numE, nT);
    tBase = zeros(nT,1);
    tPred = zeros(nT,1);
    rerouted = false(nT,1);
    
    %% predictive routing       

    % normal- if edge speed is 65 -> bestPD(1), if 50 -> bestPD(2)
    % accident- use bestPD(4) (S15)
    % construction- use bestPD(3) (S40)
    ne = height(G.Edges);
    % nScenarios = 10;
    % spdEffScenarios = zeros(ne, nScenarios);

    allCond = [];
    for col = 1:width(T_roadcond_data)
        allCond = [allCond; T_roadcond_data{:,col}];
    end
    
    %kill case sensitivity!!
    condStr = lower(string(allCond));
    
    %pull data for each type
    pNorm = sum(condStr=="normal")/ numel(condStr);
    pAcc  = sum(condStr=="accident")/ numel(condStr);
    pCons = sum(condStr=="construction")/ numel(condStr);
    if pNorm==0 && pAcc==0 && pCons==0
        %Diagnostic - tell user if the data isnt matching
        pNorm = 0.90; pAcc = 0.05; pCons = 0.05;
                fprintf('*** Diagnostic: no road-condition labels matched, using defaults: pNorm=%.2f, pAcc=%.2f, pCons=%.2f\n', ...
                 pNorm, pAcc, pCons);
    end
    cProbs = [pNorm, pAcc, pCons];
    mpH = 60;  % minutes per hour
    Gp = G;
    ne = height(Gp.Edges);
    rv = rand(ne,1);
    cp = cumsum(cProbs);
    spdEff = zeros(ne,1);
    %% actual simulation loop
    for i = 1:nT
        rv = rand(numE,1);
        spdEff = zeros(numE,1);
        for e = 1:numE
            if G.Edges.Speed(e) == 65
                nomPD = bestPD(1).pd;
            else
                nomPD = bestPD(2).pd;
            end

            if rv(e) <= cp(1)
                spdEff(e) = random(nomPD);
            elseif rv(e) <= cp(2)
                spdEff(e) = random(bestPD(4).pd);
            else
                spdEff(e) = random(bestPD(3).pd);
            end
        end

        spdEffAll(:, i) = spdEff; % Store per-trip draw
        Gp.Edges.Weight = G.Edges.Distance ./ spdEff * mpH;
        sn = trips(i,1);  
        en = trips(i,2);  
        [pB, tB] = shortestpath(Gb, sn, en, 'Method','positive');
        if isempty(pB)
            tBase(i) = NaN; tPred(i) = NaN;
            continue;
        end
        tBase(i) = tB / spdBaseAll(e) * 60; % assuming tB is really distance
        
        [pP, tP] = shortestpath(Gp, sn, en, 'Method','positive');%from word doc
        if isempty(pP)
            tPred(i) = NaN;
        else
            tPred(i) = tP / spdEff(e) * 60; % assuming tP is really distance
        end
        
        if length(pP) ~= length(pB) || any(pP ~= pB)
            rerouted(i) = true;
        end
    end
    % Diagnostic- show predictive-speed draws
        xlimrange = [0, 100];
    figure('Name','Speed Draws Comparison','NumberTitle','off');
    tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
    
    %raw nominal edge speeds
    ax1 = nexttile;
    histogram(Gb.Edges.Speed,'Normalization','pdf');
    title(ax1,'Raw Edge Speeds (Gb.Edges.Speed)');
    xlabel(ax1,'Speed (mph)');
    ylabel(ax1,'PDF');
    xlim(ax1,xlimrange);
    
    %baseline draws from fitted distributions
    ax2 = nexttile;
    histogram(spdBaseAll(:),'Normalization','pdf');
    title(ax2,'Baseline Speed Draws');
    xlabel(ax2,'Speed (mph)');
    ylabel(ax2,'PDF');
    xlim(ax2,xlimrange);
    
    %predictive draws (accidents/construction mixture)
    ax3 = nexttile;
    histogram(spdEffAll(:),'Normalization','pdf');
    title(ax3,'Predictive Speed Draws');
    xlabel(ax3,'Speed (mph)');
    ylabel(ax3,'PDF');
    xlim(ax3,xlimrange);
    
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

    %verify sim LogNormal mu and sigma
    pdTarget = makedist('Lognormal','mu',1,'sigma',1.6');

    % Kolmogorov–Smirnov test
    [h,p] = kstest(savings,'CDF',pdTarget);
    fprintf('Savings ~ LogN(1,1.6)? h=%d (0=pass), p=%.3f\n', h, p);

    %histogram + theoretical PDF overlay
    figure;
    histogram(savings, 45, 'Normalization','pdf');
    hold on;
      x = linspace(min(savings), max(savings), 200);
      plot(x, pdf(pdTarget,x), 'r-','LineWidth',1.5);
    hold off;
    title('Savings: empirical vs LogNormal(\mu=1,\sigma=1.6)');
    xlabel('Time savings (%)');
    ylabel('PDF');
    legend('Empirical','LogNormal fit');

    %QQ-plot against the target distribution
    figure;
    qqplot(savings, random(pdTarget,numel(savings),1));
    title('QQ plot: savings vs LogNormal(\mu=1,\sigma=1.6)');
    
    
    % display results
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
    histogram(savings, 10);% change bin size here as needed
    xlabel('time savings (%)'); ylabel('count');
    title('predictive routing time savings distribution');

    % plot the road condition data
    figure;
    road_data = categorical(T_roadcond_data{:,1});
    histogram(road_data)
    title('Road Condition Distribution')
    xlabel('Condition')
    ylabel('Frequency')
end

%% Part 4: validation using Monte Carlo replications
validateITS(nT, confLev, mu, sigma, uniq, bestPD, G, T_roadcond_data, nVal)

function validateITS(nT, confLev, mu, sigma, uniq, bestPD, G, Tcond, nVal)
    
    % Monte Carlo validation: replications of mean %saving
    valSave = zeros(nVal,1);
    rng('default'); 
    for k=1:nVal
        rng(k);  % different seed each rep
        [tB,tP] = runOne(nT, confLev, mu, sigma, uniq, bestPD, G, Tcond);
        valSave(k) = mean((tB-tP)./tB*100);
    end
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

    % diag plot of cumulative mean to see stability of montecarlo section
    % code
    cumMean = cumsum(valSave) ./ (1:nVal)';
    figure;
    plot(1:nVal, cumMean, 'LineWidth',1.5);
    xlabel('Replication #');
    ylabel('Cumulative mean % saving');
    title('Convergence of Monte Carlo estimate');
    yline(cumMean(end),'r--','Final mean');

    % plot CI for the montecarlo for diagnostic of the montecarlo code by
    % sensitivity
    ciWidth = 2 * tinv(1 - (1-confLev)/2, (1:nVal)'-1) .* ...
          arrayfun(@(k) std(valSave(1:k)),1:nVal) ./ sqrt(1:nVal)';
    figure;
    plot(1:nVal, ciWidth);
    xlabel('Replications');
    ylabel('CI Width (% saving)');
    title('How CI narrows with more reps');

end

function [tBase, tPred] = runOne(nT, confLev, mu, sigma, uniq, bestPD, G, Tcond)
    time_factor = 60; % hours to minutes

    % generate the same trips as in simulateITS
    trips = genlogntrips(G, nT, confLev, mu, sigma, uniq);

    % probabilities
    allC = [];
    for col = 1:width(Tcond)
        allC = [allC; Tcond{:,col}];
    end
    condStr = lower(string(allC));

    pNorm = sum(condStr=="normal")/ numel(allC);
    pAcc  = sum(condStr=="accident")/ numel(allC);
    pCons = sum(condStr=="construction")/ numel(allC);
    if pNorm + pAcc + pCons == 0
        pNorm = 0.90; pAcc = 0.05; pCons = 0.05;
        fprintf(['*** Diagnostic: no road-condition labels matched, ' ...
            'using defaults: pNorm=%.2f, pAcc=%.2f, pCons=%.2f\n'], ...
                 pNorm, pAcc, pCons);
    end
    cp = cumsum([pNorm, pAcc, pCons]);

    E = height(G.Edges);
    % draw one - mixture of speeds for every edge
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
            % accident
            pd = bestPD(4).pd;
        else
            % construction
            pd = bestPD(3).pd;
        end
        spdEff(e) = random(pd);
    end

    % baseline route and time under spdEff
    % first, we pick the baseline path by nominal speeds
    Gb = G;
    spdNom = zeros(E,1);
    for e = 1:E
        if G.Edges.Speed(e) == 65
            nomPD = bestPD(1).pd;
        else
            nomPD = bestPD(2).pd;
        end
        spdNom(e) = mean(nomPD); % average of the nominal fit
    end
    Gb.Edges.Weight = Gb.Edges.Distance ./ spdNom * time_factor;
    
    % compute tBase per trip
    tBase = nan(nT,1);
    for i = 1:nT
        sn = trips(i,1);
        en = trips(i,2);
        [pB, ~] = shortestpath(Gb, sn, en, 'Method','positive');
        if isempty(pB)
            tBase(i) = NaN;
        else
            tBase(i) = sum(Gb.Edges.Distance(pB)./Gb.Edges.Speed(pB)) * time_factor;
        end
    end
    
    % predictive route and time under the spdEff
    Gp = G;
    Gp.Edges.Weight = G.Edges.Distance ./ spdEff * time_factor;
    
    % compute tPred per trip
    tPred = nan(nT,1);
    for i = 1:nT
        sn = trips(i,1);
        en = trips(i,2);
        [pP, ~] = shortestpath(Gp, sn, en, 'Method','positive');
        if isempty(pP)
            tPred(i) = NaN;
        else
            tPred(i) = sum(Gp.Edges.Distance(pP)./spdEff(pP)) * time_factor;
        end
    end

    % drop NaNs
    ok = ~isnan(tBase) & ~isnan(tPred);
    tBase  = tBase(ok);
    tPred  = tPred(ok);
end
