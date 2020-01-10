%%% Sim scattering: Run the leaf scattering simulation and record echoes
%Created by: Anupam Kumar Gupta, Date:09/01/2015, Modified by Chen Ming,
%Date: 10/01/2015

%%% Support for volume within which leaves are distributed uniformally
%pos_param.bound = [1 10;-2 2;-2 2];% x,y,z-axes support (bounding box dimensions, meters)
pos_param.bound = [1 3;-2 -1;-0.5 0.5];
%%% Leaf distribution parameters

%leafdensity = 50; % number of leaves per cubic meter, N/m^3
leafdensity = 25;
leafsize = 2E-2; % mean leaf radius, m
alphax = 85; % mean leaf orientation, in degrees
leaf_param.leafsize_dist = 'gauss'; % can choose from 'uni' (uniform) and 'gauss' (Gaussian), but leafsize_supp has to be [upper, lower] limit
leaf_param.leafsize_supp = []; % leaf radius range in m, single value for 'gauss', [a,b] form for 'uni'
% the angle between x axis and leaf normal vector in degrees

%%% Sonar parameters
sonar_param.loc = [0 0 0];
sonar_param.peakAmp = 1;% Peak amplitude of sonar beampattern
sonar_param.freqorbmwdth = 'beamwidth';
sonar_param.freqorbmwdthVal = 30;% beamwidth in degrees (elevation and azimuth same)
sonar_param.beamcenter = [0 0];% beam center coordinates
sonar_param.noise = [];

%%% Echo parameters
n_echoes = 1;  % no of echoes to generate per paramater combination
echo_param.n_echoes = n_echoes;
%echo_param.Gthresh = 20; % gain value for thresholding (dB)
echo_param.Gthresh = 5;   % should be less than 20
%echo_param.Fs = 400; % sampling frequency (kHz)
echo_param.Fs = 100; % sampling frequency (kHz)
%echo_param.Nbins = 24E3; % no of bins b/w 60-80 kHz
echo_param.Nbins = 250; % no of bins b/w 60-80 kHz
s = zeros(n_echoes, echo_param.Nbins); % the impulse responses
for i = 1:n_echoes
                leaf_param.density = leafdensity;
                leaf_param.leafsize_supp = leafsize;
                [Sphcoord_sonar,Sphcoord_leaf,Sphcoord_leafnorm,leaf_dia] = leaf_dist_uniform1(pos_param,leaf_param,sonar_param,alphax);
                [~,~,R_sonar_leaf_filt,leaf_dia_filt,Sonargain_filt,Incident_angles] = filter_leaves_beamG(echo_param,Sphcoord_sonar,Sphcoord_leaf,Sphcoord_leafnorm,leaf_dia,sonar_param);
                s(i,:) = time_domain_impulse1(Sonargain_filt,R_sonar_leaf_filt,leaf_dia_filt,Incident_angles,echo_param);
end
