function [Sphcoord_leaf_filt,Sphcoord_leafnorm_filt,R_sonar_leaf_filt,leaf_dia_filt,Sonargain_filt,Incident_angles] = filter_leaves_beamG(echo_param,Sphcoord_sonar,Sphcoord_leaf,Sphcoord_leafnorm,leaf_dia,sonar_param)
% Filter_Leaves_BeamG: Takes in the uniformly distributed leaf locations in
% a box and filters them based on a incoming (sonar) amplitude threshold
% Inputs: (Before gain based filtering)
%        TH_sonar (Azimuth),PHI_sonar (Elevation),R_sonar (Range): Spherical coordinates of sonar in radians, Sphcoord_sonar.
%        TH_leaf,PHI_leaf,R_leaf: Spherical coordinates of leaf centers in
%            radians, Sphcoord_leaf.
%        TH_leaf_norm,PHI_leaf_norm: Spherical coordinates of leaf normals
%            in radians, Sphcoord_leafnorm.
%        leaf_dia: Vector containing leaf diameters in meters.
%        echo_param: General information about echo, details please see
%        Sim_leaf_scattering.m
%        
%        freqorbmwdth: 1) 'incident_freq' for entering frequency of incidence
%        b/w 60-80 kHz in kHz to estimate spread(sigma) of gauss beampattern 
%        based on horseshoe data, 2) 'beamwidth' for directly entering
%        beamwidth in degrees and estimate spread(sigma) of gauss
%        beampattern for the same
%        freqorbmwdthVal: frequency (kHz) b/w 60-80 kHz for 'incident_freq' in freqorbmwdth
%        and beamwidth (degrees) for 'beamwidth' in freqorbmwdth. 
%        
% Outputs: (After gain based filtering)
%        TH_sonar (Azimuth),PHI_sonar (Elevation),R_sonar (Range): Spherical coordinates of sonar in radians.
%        TH_leaf_filt,PHI_leaf_filt,R_leaf_filt: Spherical coordinates of leaf centers in
%            radians.
%        TH_leaf_norm,PHI_leaf_norm: Spherical coordinates of leaf normals
%            in radians.
%        leaf_dia_filt: Vector containing leaf diameters in meters.
%        Sonargain_filt: Vector containing sonar amplitude/gain at each leaf
% Created by: Anupam Kumar Gupta, Date:08/27/2015, Modified by Chen Ming,
% Date 11/01/2015

%%% Sonar params
TH_sonar = Sphcoord_sonar(1);
PHI_sonar = Sphcoord_sonar(2);
R_sonar = Sphcoord_sonar(3);
Sonaramp = sonar_param.peakAmp;
freqorbmwdth = sonar_param.freqorbmwdth;
freqorbmwdthVal = sonar_param.freqorbmwdthVal;
beamcenter = sonar_param.beamcenter;

%%% Leaf params
TH_leaf = Sphcoord_leaf(:,1);
PHI_leaf = Sphcoord_leaf(:,2);
R_leaf = Sphcoord_leaf(:,3);

%%% Leaf normal params
TH_leaf_norm = Sphcoord_leafnorm.TH_leaf_norm;
PHI_leaf_norm = Sphcoord_leafnorm.PHI_leaf_norm;
R_leaf_norm = Sphcoord_leafnorm.R_leaf_norm;

%%%% Find relative location of sonar w.r.t leaves and vice versa

%%% Relative distance between leaf center and sonar (range)

R_sonar_leaf = dist_points_sph_coord(R_leaf,TH_leaf,PHI_leaf,R_sonar*ones(size(R_leaf)),TH_sonar*ones(size(TH_leaf)),PHI_sonar*ones(size(PHI_leaf)));

%%% Angle between leaf center and sonar center 
Az_sonar_leaf = TH_sonar - TH_leaf;
El_sonar_leaf = PHI_sonar - PHI_leaf;

%%% Filter leaves based on beam gain
Zgain_sonar_emission = gauss_beam_sonar(Sonaramp,freqorbmwdth,freqorbmwdthVal,beamcenter,sonar_param.noise); % To understand

Zgain_sonar = Zgain_sonar_emission.^2; % Here, counting for the reception beampattern.
el = linspace(-pi/2,pi/2,181); % elevation
az = linspace(-pi,pi,361); % azimuth

%%% Sonar beampattern gain
[~,elloc_sonar] = min(abs(bsxfun(@minus,El_sonar_leaf,el)),[],2);
[~,azloc_sonar] = min(abs(bsxfun(@minus,Az_sonar_leaf,az)),[],2);
idx = sub2ind(size(Zgain_sonar),elloc_sonar,azloc_sonar);
sonarbp_gain = Zgain_sonar(idx); % sonar gain


p_ref = 2*1E-5; % reference sound pressure in air, in pascal
p0 = 2; % in pascal
r0 = 0.1; % in m, reference distance.

dis_b2 = r0^2*p0*sonarbp_gain/10^(echo_param.Gthresh/20)/p_ref; % Hassan related to equation 5 in summary
dis_b = sqrt(dis_b2);
idxfilter = find(R_sonar_leaf <= dis_b); %find the index where gain is greater than threshold

%%% Extract leaves & parameters that pass gain threshld
TH_leaf_filt = TH_leaf(idxfilter);
PHI_leaf_filt = PHI_leaf(idxfilter);
R_leaf_filt = R_leaf(idxfilter);

 
TH_leaf_norm_filt = TH_leaf_norm(idxfilter,:);
PHI_leaf_norm_filt = PHI_leaf_norm(idxfilter,:);
R_leaf_norm_filt = R_leaf_norm(idxfilter,:);

leaf_dia_filt = leaf_dia(idxfilter);
R_sonar_leaf_filt = R_sonar_leaf(idxfilter);

Sonargain_filt = sonarbp_gain(idxfilter);    

%%% Angle between leaf axis and leaf normal 
Vec1 = [R_sonar*ones(size(R_leaf_filt)) TH_sonar*ones(size(TH_leaf_filt)) PHI_sonar*ones(size(PHI_leaf_filt)) R_leaf_filt TH_leaf_filt PHI_leaf_filt];
Vec2 = [ R_leaf_norm_filt TH_leaf_norm_filt PHI_leaf_norm_filt];

Incident_angles = angle_vec_sph_coord(Vec1,Vec2);
Sphcoord_leaf_filt= [TH_leaf_filt(:) PHI_leaf_filt(:) R_leaf_filt(:)];
Sphcoord_leafnorm_filt.TH_leaf_norm_filt = TH_leaf_norm_filt;
Sphcoord_leafnorm_filt.PHI_leaf_norm_filt = PHI_leaf_norm_filt;
Sphcoord_leafnorm_filt.R_leaf_norm_filt = R_leaf_norm_filt;

end










