% function that adds measurement error to a population
function Nsampled = Func_MeasureError(N, meta, k)

    % bin sizes into more discrete size classes
    % the values of the elements in each bin are always less than the bin value
    % i.e. 3.5 will go into 4
    edges = 0:(max(meta.x)+1);
    disc_x = discretize(meta.x, edges, edges(2:end))';
    
    % summarise densities to discrete counts
    % resize to 2d
%     Nre = reshape(N,meta.meshno,meta.T*meta.RR);
%     Nsum = groupsummary(Nre,repmat(disc_x,1,meta.T*meta.RR),'sum');
%     Nsum = reshape(Nsum, size(Nsum,1), meta.T, meta.RR);

    
    Nre = reshape(N,size(N,1), prod(size(N,(2:6))));
    Nsum = groupsummary(Nre,repmat(disc_x,1, prod(size(N,(2:6)))),'sum');
    Nsum = reshape(Nsum, size(Nsum,1), size(N,2), size(N,3),...
                         size(N,4), size(N,5), size(N,6));
    
    % sample with poisson distribution
    % (sampling at each binned size with mean for that bin size)
%     Nsampled = poissrnd(Nsum);

    % sample with neg bionomial
    % k = meaure of dispersion, small = more clumping, lrg = more dispersed
    Nsampled = nbinrnd(k, k/(Nsum + k));

end