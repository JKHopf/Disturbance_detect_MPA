function [SMURFdata] = Func_ReadSmurfData(Data, Spp)
% read in SMURF data, cleans, and selects for given species
% Spp must be a char string (e.g. 'PCLA' not "PCLA")
%
% Jess Hopf
% Jan 2023
% Based on same file/function name, but in Hopf et al. 2022 
% larval periodicity paper
%
% Reads from Data\ALL_SMURF_DATA_THROUGH_2016_Raw.csv, provided by White
% July 2020 AND Reads from PISCO_UCSB_subtidal_recruitment_fish_data.1.2.csv
% from https://search.dataone.org/view/doi%3A10.6085%2FAA%2FPISCO_UCSB_Fish_Recruitment.1.3#doi%3A10.6085%2FAA%2FPISCO_UCSB_subtidal_recruitment_fish_data.1.2.csv

% % for testing:
% Spp = 'SMYS';
% Data = "PISCO_UCSB_subtidal_recruitment_fish_data.1.2.csv";

%older data
if Data == "ALL_SMURF_DATA_THROUGH_2016_Raw.csv"

    SMURFdf = readtable(Data,'TextType','string');

    % replace missing data with zeros & store in new variable
    % (only applies to classcode 'NF', = no fish??)
    SMURFdf.NumFish = fillmissing(SMURFdf.Number_fish, 'constant', 0);

    % Sum the number of fish, for the species of interest, grouped by:
    % island x site x month x day x year x NumDays x N_smurfs x lat x long
    SMURFdf_sp = SMURFdf(SMURFdf.Species == Spp,{'Island','SITE','Year',...
                  'Month','Day','NumDays','N_Smurfs', 'Long','Lat','NumFish'});

    SMURFsp = grpstats(SMURFdf_sp, {'Island','SITE','Year','Month','Day',...
                           'NumDays','N_Smurfs','Long','Lat'},{'sum'}); 
    SMURFsp = renamevars(SMURFsp, 'sum_NumFish', 'NumFish');                   

    % standardise by number of days and smurfs
    SMURFsp.NumFishSt = SMURFsp.NumFish./SMURFsp.NumDays./SMURFsp.N_Smurfs;

    % Remove NaN (from /0) as these were the days where no smurfs were deployed
    % (Note: this wasn't done in the original file)
    SMURFsp(isnan(SMURFsp.NumFishSt),:) = [];

    % Sum up over the year, grouped by Island x site x year x lat x long)
    SMURFdata = grpstats(SMURFsp(:,{'Island','SITE','Year','Long','Lat','NumFish','NumFishSt'}),...
                        {'Island','SITE','Year','Long','Lat'},'sum');


    % add species names to column names for stats
    SMURFdata = renamevars(SMURFdata, {'sum_NumFish','sum_NumFishSt'},...
                           {['Num_' Spp], ['NumSt_' Spp]});


    % rename sites that were renamed and slightly relocated
    Old_names = {'HAZ-WEST','PEL-WEST','WIL','MORSE'};
    New_names = {'HAZ','PEL','VALLEY','GULL'};
    for i = 1:length(Old_names)
        SMURFdata.SITE(SMURFdata.SITE == Old_names(i),:) = New_names(i);
    end
    
end

% newer data
% site: https://search.dataone.org/view/doi%3A10.6085%2FAA%2FPISCO_UCSB_Fish_Recruitment.1.3#doi%3A10.6085%2FAA%2FPISCO_UCSB_subtidal_recruitment_fish_data.1.2.csv
if Data == "PISCO_UCSB_subtidal_recruitment_fish_data.1.2.csv"
    
    SMURFdf = readtable(Data,'TextType','string');
    
    % Rename some variables
    SMURFdf = renamevars(SMURFdf,...
                {'site_code','num_smurfs','soak_days','classcode','year','total_fish_collected'},...
                {'SITE','N_Smurfs','NumDays','Species','Year','NumFish'});
    
    % There are no missing data
%     sum(ismissing(SMURFdf))
    
    % Add island names
        % read in old data
        SMURFdfOld = readtable("ALL_SMURF_DATA_THROUGH_2016_Raw.csv",...
            'TextType','string');
    
        % get the island x site names
        IsSite = unique(SMURFdfOld(:,{'SITE','Island'}));
        
        % "ANA-LAND","PUR", "SMI-BAY", "COJO", "JAL",
        % "ELLWOOD", & "NAPLES" are missing islands
        % ANA-SOUTH is Anacapa
        % SMI-BAY is San Miguel
        % ELLWOOD & NAPLES are in naples MR which is just off the mainland
        % COJO is mainland in Cojo Bay
        % add these to the island-site table
        % PUR & JAL is mainland, but north of Pt Conception, leave these off
        IsSite = [IsSite;
                  table(["ANA-LAND";"COJO";"SMI-BAY"],...
                        ["Anacapa";"Mainland";"San Miguel"],...
                        'VariableNames',["SITE","Island"])];
                    
        % Add site to the df
        SMURFdf = innerjoin(SMURFdf,IsSite);

    % Original way
    % Sum the number of fish, for the species of interest, grouped by:
    % island x site x year x NumDays x N_smurfs 
    SMURFdf_sp = SMURFdf(SMURFdf.Species == Spp,{'Island','SITE','Year',...
                            'Species','NumDays','N_Smurfs','NumFish'});

    SMURFsp = grpstats(SMURFdf_sp, {'Island','SITE','Year',...
                           'Species','NumDays','N_Smurfs'},{'sum'}); 
    SMURFsp = renamevars(SMURFsp, 'sum_NumFish', 'NumFish');                   

    % standardise by number of days and smurfs
    SMURFsp.NumFishSt = SMURFsp.NumFish./SMURFsp.NumDays./SMURFsp.N_Smurfs;
  
    % Remove NaN (from /0) as these were the days where no smurfs were deployed
    % (Note: this wasn't done in the original file)
    SMURFsp(isnan(SMURFsp.NumFishSt),:) = [];

    % Sum up over the year, grouped by Island x site x year x lat x long)
    SMURFdata = grpstats(SMURFsp(:,{'Island','SITE','Year','Species','NumFish','NumFishSt'}),...
                        {'Island','SITE','Year','Species'},'mean');
    
    % add species names to column names for stats
    SMURFdata = renamevars(SMURFdata, {'mean_NumFish','mean_NumFishSt'},...
                            {'NumFish','NumFishSt'});
%                            {['Num_' Spp], ['NumSt_' Spp]});

    
    
    
%     SMURFdf_sp = SMURFdf(SMURFdf.Species == Spp,{'Island','SITE','Year',...
%                             'Species','NumDays','N_Smurfs','NumFish'});
%     
%     % convert to total number of sampling units
%     SMURFdf_sp.SampUnits = SMURFdf_sp.NumDays.*SMURFdf_sp.N_Smurfs;
% 
%     SMURFsp = grpstats(SMURFdf_sp, {'Island','SITE','Year',...
%                            'Species'},{'sum'}); 
% %     SMURFsp = renamevars(SMURFsp, 'sum_NumFish', 'NumFish'); 
% 
%     SMURFsp.NumFishSt = SMURFsp.sum_NumFish./SMURFsp.sum_SampUnits;
%     SMURFdata = SMURFsp(:,{'Island','SITE','Year','Species','NumFishSt'});
%     
    



    % rename sites that were renamed and slightly relocated
    Old_names = {'HAZ-WEST','PEL-WEST','WIL','MORSE'};
    New_names = {'HAZ','PEL','VALLEY','GULL'};
    for i = 1:length(Old_names)
        SMURFdata.SITE(SMURFdata.SITE == Old_names(i),:) = New_names(i);
    end
    
end

    end





