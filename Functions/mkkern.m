function kxy = mkkern(x,F,fixparm,T)
% This function creates the IPM kernel.

% Adapted from Easterling et al. (2001) & Nickols et al. (2019)

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

% Populate the kernel ---------------------
    % Define which sizes can be fished:
    % is 1-prob of being an adult
    isjuv = 1 - normcdf(x,fixparm(5),fixparm(7)); 

    % SURVIVAL PART OF KERNEL
    % natural mortality rate
    M = fixparm(4); 
    % Mortality matrix, size x
    m = ones(size(x)).*M + (1-isjuv).*F; 
    
    % convert mortality rate to survivorship, iterate over time steps
    p1 = exp(-m*T); 

    % GROWTH PART OF KERNEL
    Linf = fixparm(1);
    k = fixparm(2);
    pmean1=Linf - (Linf - x).*exp(-k); % (do not add in x0 for the one-step growth)
    % add variability around von Bertalanffy growth
    psig1 = pmean1*fixparm(7); 

    % evaluate growth part of kernel
    p2 = normpdf(y, pmean1, psig1);

    % COMBINE GROWTH AND SURVIVAL
    % make sure no negatives
    p1 = max(0,p1);    
    p2 = max(0,p2);    
     
    % final output matrix
    kxy = p1.*p2;



