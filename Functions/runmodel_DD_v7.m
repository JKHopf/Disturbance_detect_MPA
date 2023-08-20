function [N_R, N_F] = runmodel_DD_v7(meta, fixparm, HeatWL)
% builds on v6a: vectorises replicates

% unpack meta structure
T = meta.T;
RR = meta.RR;
F = meta.F;
A = meta.A;
Rsize = meta.Rsize;
Rmean = meta.Rmean;
x = meta.x;
dx = meta.dx;
Ro = meta.Ro;
% HeatWL = meta.HeatWL;
Nint1 = meta.Nint1;
Nint2 = meta.Nint2;
whoaff = meta.whoaff;
popOC = meta.popOC;
Fun = meta.Fun;
BHa = meta.BHa;
BHb = meta.BHb;
PES = meta.PES;
Rstdv = meta.Rstdv;


% Model run set-up --------------------
    % repilcated rho vector
    Ro = repmat(Ro,1,1,RR);

    %set up matrix to hold model runs 
    % N_R and N_F is size y(size classes),t(time),RR(reps)
    % reserve scenario
    N_R = nan(length(x),T,RR); 
    % fished scenario 
    N_F = N_R;  
    

    % set intial conditions (size distributions)
    if size(Nint1,2) == 1
        N_R(:,1,:) = repmat(Nint1,1,RR);
        N_F(:,1,:) = repmat(Nint2,1,RR);
    else
        N_R(:,1,:) = Nint1;
        N_F(:,1,:) = Nint2;
    end

% Model run ---------------
    for t = 2:T %loop over years        
     
        % open or closed population
        switch popOC

            % juv 
            
            case'closed'
            % closed population with DD, recruits are then a function of adult
            % fecundity and other recruits (BH DD = recruits compete) 
            
            % incoming larvae (to the metapopulation)
            L = sum(Fun.*(N_R(:,t-1,:)+N_F(:,t-1,:)));  
            
            % --- recrtuiment variability before DD (same for both pops) ---
            % NOT USING
            % (comment out for testing with no variability, or if moving var to after DD)
%             L = max([zeros(1,1,3),normrnd(L,Rstdv)]); 
           
            % number that survive DD
            % DD survival - Beverton-Holt function (recruits affect each other)
            % This is the same for both pops under a well mixed assumption
            R = (BHa./(1+(BHa.*L)./(BHb))).*L;
            
            % --- recrtuiment variability after DD (same for both pops) ---
            % (comment out for testing with no variability, or if moving var to before DD)
            R = max([zeros(1,1,RR),normrnd(R,Rstdv)]); 

            % multiple by prob density function
            Radd = R.*Ro;
            
            case 'open'                   
            % Add variability to recruitment (same for both pops)
              Radd = Ro.*max([zeros(1,1,RR),normrnd(Rmean,Rstdv,1,1,RR)]);
            % testing - no recruitment
%               Radd = zeros(1,1,RR);

            
        end
              
        % create the kernel that caculates the probability of 
        % growing and moving from size x to y
        % (this is where fishing is included)
        % size = mesh size x mesh size x reps
        % the stochasticity (process error) for adults is embedded 
        % in this function
        kmatR = mkkern_vectorized(x, F(1), fixparm, RR, PES);
        kmatF = mkkern_vectorized(x, F(2), fixparm, RR, PES);
               
        % weight by midpoint rule (even weighting)
        kmatR = dx.*kmatR;
        kmatF = dx.*kmatF;
        
    % Heatwave disturbance
        % temporary abundances for the year 
        Ntemp_R = N_R(:,t-1,:);
        Ntemp_F = N_F(:,t-1,:); 

        switch whoaff
            case 'juvs'
            % Recruitment affected (post DD) - multiply number recruiting 
            % by the proportion that survive (HeatWL)
            % ---Juvs
            Radd_R = Radd.*A.*repmat(HeatWL(:,t-1,1),1,1,RR);
            Radd_F = Radd.*(1-A).*repmat(HeatWL(:,t-1,2),1,1,RR);            
    
            case 'adults'
            % adults affected - multiply sizes above min adult size by 
            % the proportion that survive (HeatWL)
            % ---Juvs
            Radd_R = Radd.*A;
            Radd_F = Radd.*(1-A);
            % ---Adults
            % find min adult size, based on max R size
            min_adt = find(x>Rsize,1);
            Ntemp_R(min_adt:end,:,:) = N_R(min_adt:end,t-1,:).*repmat(HeatWL(min_adt:end,t-1,1),1,1,RR);
            Ntemp_F(min_adt:end,:,:) = N_F(min_adt:end,t-1,:).*repmat(HeatWL(min_adt:end,t-1,2),1,1,RR);       
        end
            

        % Run for each protection scenario
        % keep >= 0 (issue from + var recruits)
        
        % RESERVE
        N_R(:,t,:) = pagemtimes(kmatR,Ntemp_R) + Radd_R;
        N_R(:,t,:) = max(0,N_R(:,t,:)); 

        % FISHING
        N_F(:,t,:) = pagemtimes(kmatF,Ntemp_F) + Radd_F;
        N_F(:,t,:) = max(0,N_F(:,t,:)); 

    end