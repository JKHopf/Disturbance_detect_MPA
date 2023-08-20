function [N_init] = load_initN(meta, site_name)
   
% define edges of mesh
edges = meta.x - meta.dx/2;

% Get site specific info
% MPAnew specifies if the location (within site) was assigned as an MPA in
% 2007
switch site_name
    case 'Big_Creek'
        load('Init_Dist_Data\SMYS_Big_Creek_pre2007_data_sorted.mat')
        MPAnew = logical([0 1 0 0]); 
        
    case 'White_Rock'
        load('Init_Dist_Data\SMYS_White_Rock_pre2007_data_sorted.mat')
        MPAnew = logical([1 1 1 1 0 0]); 

    case 'Pt_Lobos'
        load('Init_Dist_Data\SMYS_Pt_Lobos_pre2007_data_sorted.mat')
        MPAnew = logical([0 0 1 0 0 0 0]);
end

% Get data into a histogram 
% dims (lengths, sites, years)
N = IPM_histo(D_str.('SMYS'),Years,Site_Names,edges);


%divide total fish numbers by number of transects
for i = 1:length(Years)
    for j = 1:length(Site_Names)
        NT = D_str.('SMYS').(Site_Names{j})(i).data.numtrans;
        N(:,j,i) = N(:,j,i)./NT;
    end
end

%start with starting distribution from 2007 for site that is a new MPA
N_init1 = N(:,MPAnew,find(Years == 2007)); 

% Determine size classes to be included (account for size selectivity of sampling)
% see White et al. 2016 Ecological Applications
    % normal cdf parameters specifying
    % ogive giving probability of actually being observed in the kelp forest. 
    % For blue rockfish this is from onto_ogive.m base on Rick Starr's data.  
    Ogive_b = [28 7];                             
                                
    if isnan(Ogive_b(1)) 
        % if there is no ogive (for spp other than blue rockfish)
        Ogive = ones(size(x));
    else
        % probability of observation in the kelp forest
        Ogive = 1-normcdf(meta.x,Ogive_b(1),Ogive_b(2)); 
    end
    
    % probability of being observed
    OKlen = Ogive; 
    
    % update initial distribution
    N_init = N_init1./OKlen';


% plot distribution
% figure
% hold on
% plot(N_init, 'b', 'LineWidth', 1.5)
% ylim([0,4])
% xlim([0,50])



