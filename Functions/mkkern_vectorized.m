function kxy = mkkern_vectorized(x,F,fixparm,RR, PES)
% This function creates the IPM kernel.

% Adapted from Easterling et al. (2001) & Nickols et al. (2019)
% updated from mkkern to allow for vectorisation 
% (i.e., mutiple simultaneous replicates)
% + add process error to adult (natural mortality varies)

% Set up the integration mesh kernel ----------------------
    % For MLPA monitoring model
    y = x;
    %this creates a vector (y) that is equal to x

    [x,y] = meshgrid(x,y); % Matlab built in function
    % x is an array (original x by original x) with each
    % row equal to originalvector x
    % y is an array (original y by original y) with each 
    % column equal to original vector y
    % X corresponds to size at time t
    % Y corresponds to size at time t+1
    % RR = replicates

% Populate the kernel ---------------------
    % Define which sizes can be fished:
    % is 1-prob of being an adult
    isjuv = 1 - normcdf(x,fixparm(5),fixparm(7)); 

    % SURVIVAL PART OF KERNEL

    % natural mortality rate
    % process error
        % process error
        M = max([zeros(1,1,RR),normrnd(fixparm(4),PES,1,1,RR)]);

    % Mortality matrix (natural + fishing if fished age), size x by x by reps
    m = ones([size(x),RR]).*M + repmat((1-isjuv).*F,1,1,RR); 
    
    % convert mortality rate to survivorship
    % p1 size = length x by length y by reps
    p1 = exp(-m); 

    % GROWTH PART OF KERNEL
    Linf = fixparm(1);
    k = fixparm(2);
    pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)
    % add variability around von Bertalanffy growth
    psig1 = pmean1*fixparm(7); 

    % evaluate growth part of kernel
    % p2 size = length x by length y
    p2 = normpdf(y, pmean1, psig1);

    % COMBINE GROWTH AND SURVIVAL
    % make sure no negatives
    p1 = max(0,p1);    
    p2 = max(0,p2);    
     
    % final output matrix
    kxy = p1.*repmat(p2,1,1,RR);



