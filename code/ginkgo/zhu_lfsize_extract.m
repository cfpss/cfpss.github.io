% The zhu version is based on Chen Ming's version. 
% added some tuning parameters. 

% To sample leaf radii from truncated normal
% Created by Chen Ming, Date:  11/01/2016
function leaf_dia = zhu_lfsize_extract(leaf_param, n_leaf, sd_ratio,left_trun, right_trun)
% Input: 
% sd_ratio: chen ming's default 0.1.
% left_trun, right_trun: chen ming's defaul (0,1)

leaf_size_dist = leaf_param.leafsize_dist;
mleaf_size = leaf_param.leafsize_supp;
leaf_size_sup = leaf_param.leafsize_supp;

switch leaf_size_dist
    case 'uni'
        
        % Generate the leaf sizes from a uniform distribution
        
        % parameters for uniform distribution
        
        if (length(leaf_size_sup) == 2)
            a_leaf = leaf_size_sup(1);
            b_leaf = leaf_size_sup(2);
        else
            error('Wrong input');
        end
        
        leaf_dia = random('Uniform',a_leaf,b_leaf,n_leaf,1);

    case 'gauss'
        
        
        %%% Generate the leaf sizes from a gaussian distribution
        
        % parameters for normal distribution
        
            leaf_rds_mean = mleaf_size;
            leaf_rds_sigma = sd_ratio*mleaf_size;

        
         pd = makedist('Normal','mu',leaf_rds_mean,'sigma',leaf_rds_sigma);
          tr = truncate(pd,left_trun,right_trun);
          leaf_dia = random(tr, n_leaf,1);
end
end