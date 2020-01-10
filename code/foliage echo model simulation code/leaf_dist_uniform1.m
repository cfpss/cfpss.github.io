function [Sphcoord_sonar,Sphcoord_leaf,Sphcoord_leafnorm,leaf_dia] = leaf_dist_uniform1(pos_param,leaf_param,sonar_param, alphax)
%LEAF_DIST_UNIFORM: Draws the leaf locations from a uniform distribution &
%returns the sonar,leaf centers,leaf normal spherical coordinates, leaf
%diameters, distance/range between sonar and leaf.
% Inputs:
%        pos_param, the struct contains the info of leaf locations.
%        pos_param.bound: x,y,z-axes support (bounding box dimensions, meters)
%
%        leaf_param, the struct contains the info about the distribution of
%        leaves
%        leaf_param.
%            leafloc_dist, the distribution of leaf locations
%            leafsize_dist, the distribution of leaf radii
%            leafsize_supp, the leaf radius support for uniform distribution
%
%        sonar_param, the struct contains the info about sonar, details of
%        each field please see comments in Sim_leaf_scattering.m
%        alphax, the mean leaf orientation angle, standard deviation: 5
%        degrees
% Outputs:
%         TH_sonar (Azimuth),PHI_sonar (Elevation),R_sonar (Range): Spherical coordinates of sonar in radians, Sphcoord_sonar.
%         TH_leaf,PHI_leaf,R_leaf: Spherical coordinates of leaf centers in
%             radians, Sphcoord_leaf.
%         TH_leaf_norm,PHI_leaf_norm: Spherical coordinates of leaf normals
%             in radians, Sphcoord_leafnorm.
%         leaf_dia: Vector containing leaf diameters in meters.
%Created by:Anupam Kumar Gupta,Date: 1/12/2015, Modified by Chen Ming, Date
%10/01/2015


%%% Bounding Box params
ax = pos_param.bound(1,1);
bx = pos_param.bound(1,2);
ay = pos_param.bound(2,1);
by = pos_param.bound(2,2);
az = pos_param.bound(3,1);
bz = pos_param.bound(3,2);
%%% Leaf params
leaf_density = leaf_param.density;
leaf_size_dist = leaf_param.leafsize_dist;
mleaf_size = leaf_param.leafsize_supp;

%%% Sonar params
sonar_loc = sonar_param.loc;

%%% Find number of leaf  from the envelope volume and leaf density
Env_vol = (bx-ax)*(by-ay)*(bz-az);
n_leaf = round(Env_vol * leaf_density);


%%% Leaf locations:Uniform or Gaussian
        
        xleafloc = random('Uniform',ax,bx,n_leaf,1);
        yleafloc = random('Uniform',ay,by,n_leaf,1);
        zleafloc = random('Uniform',az,bz,n_leaf,1);
z_for_gazebo = zleafloc + 3;
%%% Generate the leaf sizes (radius)

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
            leaf_rds_sigma = 0.1*mleaf_size;

          pd = makedist('Normal','mu',leaf_rds_mean,'sigma',leaf_rds_sigma);
          tr = truncate(pd,0,1);
          leaf_dia = random(tr, n_leaf,1);
end

%%% Sonar location
if (length(sonar_loc) == 3)
    
    x_sonar = sonar_loc(1);
    y_sonar = sonar_loc(2);
    z_sonar = sonar_loc(3);
else
    error('Wrong input');
end

%%% Translate the coordinate system origin to sonar center
x_sonartrans = x_sonar - x_sonar;
y_sonartrans = y_sonar - y_sonar;
z_sonartrans = z_sonar - z_sonar;

LeafcenM = [xleafloc(:) yleafloc(:) zleafloc(:)];
Transmatrix = repmat([-x_sonar -y_sonar -z_sonar],n_leaf,1);
Leafcentrans = LeafcenM+Transmatrix; % Transmatrix is 0 matrix here.


% Convert the leaf, leaf surface normals & sonar coordinates from cartesian to spherical
[TH_leaf_sonar,PHI_leaf_sonar,R_leaf_sonar] = cart2sph([Leafcentrans(:,1);x_sonartrans],[Leafcentrans(:,2);y_sonartrans],[Leafcentrans(:,3);z_sonartrans]);


%%% Leaf
TH_leaf = TH_leaf_sonar(1:end-1);
PHI_leaf = PHI_leaf_sonar(1:end-1);
R_leaf = R_leaf_sonar(1:end-1);

%%% Leaf normal vectors
mu = alphax;
sigma = 5; % the std of the angle between pulse direction and leaf normal vector, in degrees
pd_alpha = makedist('Normal','mu',mu,'sigma',sigma);
tr_alpha = truncate(pd_alpha,0,90);
alpha_deg = random(tr_alpha,n_leaf,1);
alpha_rad = deg2rad(alpha_deg);

theta = 2*pi*rand(n_leaf,1);
X = ones(n_leaf,1);
Y = tan(alpha_rad).*sin(theta);
Z = tan(alpha_rad).*cos(theta);

[TH_leaf_norm, PHI_leaf_norm, ~] = cart2sph(X,Y,Z);

R_leaf_norm = ones(n_leaf,1);
%%% Sonar
TH_sonar = TH_leaf_sonar(end);
PHI_sonar = PHI_leaf_sonar(end);
R_sonar = R_leaf_sonar(end);

%%% Outputs
Sphcoord_leaf= [TH_leaf(:) PHI_leaf(:) R_leaf(:)];
Sphcoord_leafnorm.TH_leaf_norm = TH_leaf_norm;
Sphcoord_leafnorm.PHI_leaf_norm = PHI_leaf_norm;
Sphcoord_leafnorm.R_leaf_norm = R_leaf_norm;
Sphcoord_sonar = [TH_sonar PHI_sonar R_sonar];

end

