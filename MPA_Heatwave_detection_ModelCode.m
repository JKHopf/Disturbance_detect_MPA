% MPA_Heatwave_detection_v6.m
%
% Author Jess K Hopf
% Feb 2023
%
% The purpose of this model is:
% 1) to assess how the 2014-16 marine heatwave in central California
%    affected our ability to detect MPA effects. 
%    - adults and juvs
% 2) to assess how reserves affect our ability to detect resilience.
%
% Model breakdown:
% - Intergral projection model IPM. 
%   (see Easterling et al. 2000, Ecology, White et al. 2016, Ecol Apps &
%   Nickols et al. 2019)
% - Single-species, multiple size 'classes'
%       - Case study spp is Blue rockfish (SMYS)
% - Single popualtion model, two protection scenarios:
%       1) Reserve popualtion (fishing stops)
%       2) Fished populaiton (fishing continues at same rate as before)
% - open population with stochastic recrtuitment
% - density independent model (as recruitment is decoupled from local 
%   population size)
% 
% See v0 for the White et al script where sections came from 
%
% ***This version is designed to build a table of simulation outputs which 
% is then read into R for plotting
%
% Builds on v7:
%   - updates code for parallel computing
%   - population options: open (D-indep, external recruitment) and   
%                         closed (DD, 10% reserve) 
%   - brings back process error into adults
%   - larval variability is the same for both pops (regardless of open or
%     closed population)

%% Model Run ----------------------------------------------------
tic
% Top-level stuff --------------------

    % clear all
    clear
    
    % add path
    addpath('Functions\')
    addpath('Outputs\')
    
    % Run Model Function:
    RMfunc = str2func('runmodel_DD_v7'); 
    
% ------------ Scenarios ---------------
% these need to be manually changed for the relevant scenario - sorry, Im
% not that good at coding!

% Basic model parameters:
    % number of runs/replicates                                      
    meta.RR = 500; %  2;% 

    % Proportion of area in reserves
    meta.A = 0.1; %

% Demographic:
    % is the population open or closed? 
    meta.popOC = 'open'; % 'closed'; %     

    % variability in recruitment
    % variability in recruitment to pop (YOY data) = 5.4*10^3
    % set to 0 for no variability
    meta.Rstdv = 5.4*10^3;  

    % process error in adults (variability in adult mortlaity year to year)
    % sigma (stddev) value = SE = 0.014 in Dick et al. (2017; Table 51 pg 121) (recent stock assessment)
    % set to 0 for no variability
    meta.PES = 0.014; 

    % Natural mortality (see fixparm below)
    M = 0.119; %  0.119*0.5; % 


% Disturbance:      
    % who is affected by the disturbance?
    meta.whoaff = 'adults'; %  'juvs'; % 

    % Length of disturbance (years)
    LHW_vec =  [1,3,5]; % 2; %
    
    % Time since MPA established (years; 'age of MPA+1')
    THW_vec = [1,5,9];%  2; %
    


% ---------------- Fixed model stuff ------------------

% Recruitment statistics/data ----------------------
% estimated from data in Nickols et al. 2019
% units in cm, and TL (total length)

    % mean of 'recruits' size (estimated empirically from the data)
    Recruits.meansize = 7.75;
    
    % sd of 'recruits' size
    Recruits.sdsize = 1.15; 
    
    % max size of recruits in data (YOY)
    meta.Rsize = 10; 
    
    % Blue rockfish (SMYS):
    %        [Linf      k     x0     M    Lfish   Lmat   Lvar   c   d];
    % values in cm
    % convert, x0 = Linf*(1-exp(K*t0))
    % Lvar = value of CV (for adults size) from stock assessment.
    % Lmat = length at 50% maturity
    % For notes on parameter sources, see
    % Life_History_Paramter_for_RF_models.xlsx 
    % (Created by K. Nickols, see also appendix S1, Table 1A, Nickoles 2019)
    % c & d = fecundity at length (mm) paras from Dick et 
    %         al. (2017, Fish Res), Table 6
    % M updated to match Dick et al. (2017; Table 51 pg 157) (recent stock assessment)
    fixparm = [38.150142 0.172 6.2533 M 21.0295 27.086 0.1 exp(-15.561) 4.816];

    % Fishing rate (Pt Lobos, Nickols 2019)
    F = 0.19; 
    
        
% IPM parameters ------------------------
% See White et al (2016) Ecol Apps, specifically appendix S3.

    % number of grids
    meta.meshno = 100;
    % min mesh size (min length for fish)
    meta.meshmin = 0;
    % max mesh size (larger than any fish is likely to grow)
    meta.meshmax = fixparm(1)*2; 
    
    % build mesh (size groupings = lengths)
    meta.x = linspace(meta.meshmin,meta.meshmax,meta.meshno);
    % change in x (mesh/grid size)
    meta.dx = meta.x(2)-meta.x(1);
       

% Recrtuiment ------------------------

    % Make pdf vector of recruits (rho)
    meta.Ro = normpdf(meta.x,Recruits.meansize,Recruits.sdsize)';  
        
% Funcudity at length
    % x10 to convert cm -> mm
    meta.Fun = Func_fecunds(meta.x.*10, fixparm(8), fixparm(9));
    % set age at maturiy
    meta.Fun(1:find(meta.x>=fixparm(6),1)) = 0; 

% Open:
    % mean number of recruits for each year (arbitary/scaling factor)
    % doesn't affect post-reserve build-up or post-disturbance recovery dynamics
    % dictates abundance of unfished/no variability population (~10xRmean @ 100% resserves)
    meta.Rmean = 10^4; % 0.65; % 


% Closed: Beverton-Holt Para values
    % calculuate life-time egg production 
        % get fecundities for lengths < Linf  <-------------------------------------
        F_real = meta.Fun(meta.x<=fixparm(1));
        % calc LEP
        LEP = sum(cumprod([1;repmat(exp(-fixparm(4)), length(meta.Fun)-1, 1)]).*meta.Fun);
    
    % a = slope at/near zero/origin
    % set so that if population drops below 25% of the unfished LEP then
    % population declines
    meta.BHa = 1/(0.25*LEP); % 1/(0.25*10); % 
    % b = max density of settlers
    meta.BHb = meta.Rmean; % 10; 


% ------------------------ Intial distribution --------------------  

% Make intial distribution base on steady state ---------------------
    % starting distribution
    meta.Nint1 = repmat(10,meta.meshno,1);
    meta.Nint2 = repmat(10,meta.meshno,1);
    
    % timescale to run model
    meta.T = 50;  
    
    % fishing pressures in both pops (res, fished)
    meta.F = [F,F]; % [0,0]; % 
    
    % vector of added mortality due to heatwave
    % dims = mesh number, run time at least, number of populations (2)
    % one = no effect
    % set no disturbance
    HeatWL = ones(meta.meshno,100,2);
    
    % run model   
    [NRint1, NFint1] = RMfunc(meta, fixparm, HeatWL);
    [NRint2, NFint2] = RMfunc(meta, fixparm, HeatWL);
    
    % get init vector
    % one vector per replicate (so different starting conditions)
    % (only fished pop, as pre-res)
    meta.Nint1 = squeeze(NRint1(:,end,:));
    meta.Nint2 = squeeze(NFint1(:,end,:));
    
    % check outputs
%     figure(32)
%     hold on
%     plot(squeeze(sum(NFint1,1)),'b')
%     plot(squeeze(sum(NRint1,1)),'g')
%     ylabel('total pop abundance')
%     xlabel('time')

 
% --------------------- Model Scenarios w/reserves -----------------------

% Fishing pressure post-reserves:
    % fishery squeeze with reserve est
    F = F/(1- meta.A);
    % fishing pressures in both pops (res, fished)
    meta.F = [0,F];

% re-set timescale to run model
% use 30, 60 is only if need to run longer 
meta.T = 30; %  60;% 

% intensity vector
% (0.1 = reduction to 10%, 1 = no reduction/heatwave)
% less intervals for SA
switch meta.whoaff 
    case 'juvs'
    MHW_vec = 0.0:0.05:1; % 0.2;% 
    case 'adults'
    MHW_vec = 0.2:0.1:1; % 0.2:0.05:1; % 0.2; % 
end

% pre-assign variables
% dim = mesh length, time, reps, length of MR age vector,
%       length of distlength vector, length of intensity vector
NR = nan(length(meta.x),meta.T,meta.RR,length(THW_vec),...
         length(LHW_vec),length(MHW_vec));
NF = NR;


% run model
for a = 1:length(THW_vec)  % for all MPA ages
    THW = THW_vec(a);
    
    for l = 1:length(LHW_vec)  % for all dist lengths
        LHW = LHW_vec(l); 
    
       parfor i = 1:length(MHW_vec)  % for all dist intensities 
            % vector of added mortality due to heatwave
            % dims = mesh number, run time at least, number of populations (2)
            % one = no effect
            % set up with no disturbance
            HeatWL = ones(meta.meshno,100,2);
            
            % add disturbance
            % affects both reserve and fished area equally
            HeatWL(:,THW:(THW+LHW-1),:) = MHW_vec(i);

            % population abundances
            [NR(:,:,:,a,l,i), NF(:,:,:,a,l,i)] = RMfunc(meta, fixparm, HeatWL);
            
        end
    end 
end
        
%           %  check outputs
%             figure(9)
%             hold on
%             plot(repmat((-20:(meta.T-1))',1,meta.RR),...
%                 [squeeze(sum(NRint1(:,(end-20):(end-1),:),1));squeeze(sum(NR,1))],...
%                 'Color',[0.72, 0.05, 0.28, 0.8], 'LineWidth', 1)
% %             plot(repmat((-20:(meta.T-1))',1,meta.RR),...
% %                 [squeeze(sum(NFint1(:,(end-20):(end-1),:),1));squeeze(sum(NF,1))],...
% %                 'Color',[0.40, 0.75, 0.76, 0.8], 'LineWidth', 1)
%             plot(repmat((-20:(meta.T-1))',1,meta.RR),...
%                 [squeeze(sum(NFint1(:,(end-20):(end-1),:).*(meta.A/(1-meta.A)),1));...
%                  squeeze(sum(NF.*(meta.A/(1-meta.A)),1))],...
%                 'b', 'LineWidth', 1)            
%             xline(0, '--k')
%             grid(gca,'minor')
%             grid on
            
            
%             figure(1234)
%             hold on
%             plot(repmat((1:meta.T)',1,meta.RR)-1,...
%                 squeeze(sum(NR,1))./7443.4,...%./squeeze(sum(NR(:,1),1)),...
%                 'Color','g', 'LineWidth', 0.1)
%             plot(repmat((1:meta.T)',1,meta.RR)-1,...
%                 squeeze(sum(NF,1)).*(meta.A/(1-meta.A))./7443.4,...%./squeeze(sum(NF(:,1),1)),...
%                 'Color','c', 'LineWidth', 0.1)
%             xline(0, '--k')
%             grid(gca,'minor')
%             grid on



% calculate Biomass
NR_bio = NR.*meta.x';
NF_bio = NF.*meta.x';

% Save variables (for testing/working AUCs)
%  save(datestr(now, 'yyyy-mm-dd') + "Variables_V7_results_" + "_" + meta.whoaff + "_" + meta.popOC + "OnePopTest" )
toc


% ------------------- Calculate AUCs & add to table -----------------------
% find proportion of false positives (FP) and true positives (TP)   
% calculate area under the curve (AUC)

tic

% include measurement error -----------
% affects all years post reserves 
% (we'll choose sampling years in the AUC step)
% see function script for more info
    % aggregation para (degree of clumping, small = more clumping)
k = 1; 
NRsamp = Func_MeasureError(NR, meta, k);
NFsamp = Func_MeasureError(NF, meta, k);


%     % check output
%     figure
%     hold on
%     plot(repmat((1:meta.T)',1,meta.RR),...
%         squeeze(sum(NR,1)),... %;./squeeze(sum(NR(:,1),1)),...
%         'Color',[0.72, 0.05, 0.28, 0.8], 'LineWidth', 1)
%     plot(repmat((1:meta.T)',1,meta.RR),...
%         squeeze(sum(NF.*(meta.A/(1-meta.A)),1)),... %;./squeeze(sum(NF(:,1),1)),...
%         'Color',[0.40, 0.75, 0.76, 0.8], 'LineWidth', 1)
%     plot(repmat((1:meta.T)',1,meta.RR),...
%         squeeze(sum(NRsamp,1)),... %;./squeeze(sum(NRsamp(:,1),1)),...
%         'Color',[0.72, 0.05, 0.28, 0.8], 'LineWidth', 1, 'LineStyle','--')
%     plot(repmat((1:meta.T)',1,meta.RR),...
%         squeeze(sum(NFsamp.*(meta.A/(1-meta.A)),1)),... %;./squeeze(sum(NFsamp(:,1),1)),...
%         'Color',[0.40, 0.75, 0.76, 0.8], 'LineWidth', 1, 'LineStyle','--')
%     xline(0, '--k')
%     ylabel('Normalised abundance')
%     xlabel('Years since MPA est')
%     grid(gca,'minor')
%     grid on
   

% Calculate AUC -----------

% set sampling years (post disturbance)
sampYrs = 0:(meta.T-4);

% abundance or biomass?
Measure = 'Biomass'; % 'Abund'; % 

switch Measure
    case 'Biomass'
    NR_AUC = NRsamp.*(1:(max(meta.x)+1))'; % NR_bio; %  
    NF_AUC = NFsamp.*(1:(max(meta.x)+1))'; % NF_bio; % 
    case 'Abund'
    NR_AUC = NR; % 
    NF_AUC = NF; % 
end

% normalise fished area to reseves 
NF_AUC = NF_AUC.*(meta.A/(1-meta.A));

% Create empty table
% NOTE: samp_time = time since disturbance (0 = yr of disturbance)
%       Yr_sampled = time since MPA established (0 = year established)
TabSize = [length(THW_vec)*length(LHW_vec)*length(sampYrs)*length(MHW_vec),15];
varNames = ["Open_Closed","WhoAff","Var_R","Var_A","M","K",...
            "Yr_Dist_Start","Dist_Length","Yr_Dist_End","Yr_Sampled",...
            "Samp_time","Dist_Impact","IO","BA", "BACI"];
varTypes = [repmat("string",1,2),repmat("double",1,13)];
Tab = table('Size', TabSize,'VariableTypes',varTypes,'VariableNames',varNames);
% Tab(row#,col#)

% count itnerations
n = 1;

for a = 1:length(THW_vec)  % for all MPA ages
    THW = THW_vec(a);
            
    for l = 1:length(LHW_vec)  % for all dist lengths
        LHW = LHW_vec(l); 

        for i = 1:length(MHW_vec) % for each intensity
            % inside/outside
                
            for s = 1:length(sampYrs) % for each sampling year
            % (years since end of heatwave; must be at least 1yr since)
            SampT = sampYrs(s); % THW+LHW-1+sampYrs(s);  
            
                NFs = squeeze(sum(NF_AUC(:,s,:,a,l,i),1));
                NF1 = squeeze(sum(NF_AUC(:,1,:,a,l,i),1));
                NRs = squeeze(sum(NR_AUC(:,s,:,a,l,i),1));
                NR1 = squeeze(sum(NR_AUC(:,1,:,a,l,i),1)); 
            
                % inside/outside
                [FP1,TP1] = ROC(NFs,NRs);
                [~,Tind1] = unique(FP1);
                AUC1 = -trapz(FP1,TP1);
                
                % before/after
                % NOTE: when initial ditribution is set, it does not have a 
                % variance (i.e., no distribution) making BA comparisons with 
                % AUC not possible.  
                [FP2,TP2] = ROC(NR1,NRs);
                [~,Tind2] = unique(FP2);
                AUC2 = -trapz(FP2,TP2);

                %BACI
                % (outside-after/outside-before),(inside-after/inside-before)    
                [FP3,TP3] = ROC(NR1-NF1,NRs-NFs); 
                [~,Tind3] = unique(FP3);
                AUC3 = -trapz(FP3,TP3);
                
                
%                 % Plot ROCs
%                 figure(325)
%                 hold on
%                 
%                 subplot(3,1,1)
%                 hold on
%                 plot(FP1,TP1)
%                 title('IO')
%                 plot([0,1],[0,1],":k",'linewidth',1)
%                 
%                 subplot(3,1,2)
%                 hold on
%                 plot(FP2,TP2)
%                 ylabel('True Positive')
%                 title('BA')
%                 plot([0,1],[0,1],":k",'linewidth',1)
%                 
%                 subplot(3,1,3)
%                 hold on
%                 plot(FP3,TP3)
%                 xlabel('False Positive')
%                 title('BACI')
%                 plot([0,1],[0,1],":k",'linewidth',1)
                
                
                % assign to table
                % "Open_Closed","WhoAff",
                % "Var_R","Var_A", "M","K",
                % "Yr_Dist_Start","Dist_Length","Yr_Dist_End","Yr_Sampled",
                % "Samp_time","Dist_Impact","IO","BA", "BACI"
                Tab(n,:) = {meta.popOC, meta.whoaff, ...
                            meta.Rstdv, meta.PES, M, k,...
                            THW-1, LHW, THW+LHW-1, sampYrs(s),...
                            sampYrs(s)-(THW+LHW-1), MHW_vec(i), ...
                            AUC1, AUC2, AUC3};
                
                % update interation count
                n = n+1;
                                                
            end
        end 
    end
end

toc 

% Save csv
% writetable(Tab, datestr(now, 'yyyy-mm-dd') + "_V7results_" + meta.popOC + "_" + meta.whoaff + "_BASELINE.csv")



        

