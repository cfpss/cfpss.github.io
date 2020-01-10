%%% Record echoes from L-system ginkgo trees
% Created by Chen Ming, Date: 01/01/2017
clc;clear all;
close all;
addpath('ginkgoLsys');
addpath('simuroutines');

leafsize = 2.5E-2; % in meter


t = (0:24E3-1)/400E3*1E3;

%%% Parameters you don't change
leaf_param.leafloc_dist = 'uni';
leaf_param.leafsize_dist = 'gauss';
leaf_param.leafsize_supp = []; % leaf radius range in m
% the angle between x axis and leaf normal vector in degrees

%%% Sonar parameters

peakAmp = 1;% Peak amplitude of sonar beampattern
freqorbmwdth = 'beamwidth';

%sonar_param.beamcenter = [0 0];% beam center coordinates
beam_rot = 0;% beam rotation/orientation
%sonar_param.noise = [];

%%% Echo parameters
n_echoes = 100;
echo_param.n_echoes = n_echoes;
echo_param.Gthresh = 20; % gain value for thresholding (dB)
echo_param.Fs = 400; % sampling frequency (kHz)
echo_param.Nbins = 24E3; % no of bins b/w 60-80 kHz


leaf_param.leafsize_supp = leafsize;

%%% Parameters used in building L-system tree
nReps = 4; % growth level
bRatio = 0.62; % contraction ratio of branches other than center branch
angle = 50; % branching angle
Rb = 2; % the length of initial branch
lfnode = 4; % the number of leaf nodes per branch
Raa = 0.5; % lift the copy half length of the initial branch
Rab = 0.826; % Contraction ratio of center branch

%%% Echo parameters
bwidth = 10; % sonar beamwidth changing from 10 degrees to 50 degrees
y = zeros(length(bwidth),24000); % echoes


ys = 5;
zs = 6;
lfnum = zeros(length(bwidth),1);

for j = 1:length(bwidth)
    figure
    j
    freqorbmwdthVal = bwidth(j);
    sonar_loc = [0  ys zs];
    boxcenter = [0, 0, zs];

    [BR, lf, ntip, Itip] = GMT3(nReps, bRatio, angle, Rb, lfnode);

    %%% The second set of branches (tree)
    SecondBR = zeros(size(BR,1), size(BR,2));
    
    %%% Plot branches of the "two" trees
    RotM = rotz(90);
    N = size(BR,1);
    
        for n = 1:N
            hold on
            endrotated = RotM*BR(n, 2:4)';
            endrotated(3) = endrotated(3) + Raa*Rab*Rb; 
            tiprotated = RotM*BR(n, 5:7)';
            tiprotated(3) = tiprotated(3) + Raa*Rab*Rb;
      %      cylinder2P(BR(n, 8), 10, BR(n, 2:4), BR(n, 5:7));
       %     cylinder2P(BR(n, 8), 10, endrotated', tiprotated');   
            SecondBR(n,:) = [BR(n,1) endrotated' tiprotated' BR(n,8)];
            
        end
    newBR = [BR; SecondBR];  % new branches coordinates from two trees
    newItip = [Itip; Itip+length(Itip)/4];

    %%% Acquire echoes from the leaves of the "two" trees
    secondntip = RotM*transpose(ntip);
    secondntip = secondntip';
    secondntip(:,3) = secondntip(:,3) + Raa*Rab*Rb;
    newntip = [ntip; secondntip];
    treeLeaf2_loc = RotM*lf(:,1:3)';
    treeLeaf2_norm = RotM*lf(:,4:6)';
    treeLeaf2_loc_t = treeLeaf2_loc';
    treeLeaf2_loc_t(:,3) = treeLeaf2_loc_t(:,3) + Raa*Rab*Rb;
    newLeaf = [lf; [treeLeaf2_loc_t treeLeaf2_norm']];

    xleafloc = newLeaf(4:6,1);
    yleafloc = newLeaf(4:6,2);
    zleafloc = 5 + newLeaf(4:6,3);
    Xlfnorm = newLeaf(4:6,4);
    Ylfnorm = newLeaf(4:6,5);
    Zlfnorm = newLeaf(4:6,6);
    n_leaf = size(newLeaf(4:6,:),1);
    leaf_dia = lfsize_extract(leaf_param, n_leaf);

    R_leaf_norm = ones(n_leaf,1);
    [TH_leaf_norm, PHI_leaf_norm, ~] = cart2sph(Xlfnorm, Ylfnorm, Zlfnorm);


    hold on, waterdrop_sonarbeam(sonar_loc, boxcenter, freqorbmwdthVal);

    [Sphcoord_sonar, Sphcoord_leaf, beamcenter] = transcor(sonar_loc, boxcenter, xleafloc, yleafloc, zleafloc);
    [idxfilter, ~,~,R_sonar_leaf_filt,leaf_dia_filt,~,Sonargain_filt,Incident_angles] = filter_leaves_beamG(echo_param,Sphcoord_sonar,Sphcoord_leaf, TH_leaf_norm, PHI_leaf_norm, R_leaf_norm,leaf_dia,peakAmp, freqorbmwdth,beam_rot,freqorbmwdthVal , beamcenter);
    lfnum(j) = length(idxfilter);
    imp = time_domain_impulse1(Sonargain_filt,R_sonar_leaf_filt,leaf_dia_filt,Incident_angles,echo_param);
 
    %%% Plot leaves that lied in the beam

    % Prepare the coordinates of a sphere to represent a cluster (bundle)
    % of ginkgo leaves
    [xx,yy,zz] = sphere(10);
    xx = 5E-2*xx;
    yy = 5E-2*yy;
    zz = 5E-2*zz;
    
    in_node = newItip(idxfilter);
    in_node = unique(in_node);



    %%% Plot leaves that were not illuminated by the beam
    newItip(idxfilter,:) = [];
    out_node = unique(newItip);
    color = [83/255, 132/255, 10/255];
     %for greencolor   
    xtmp = xleafloc;
    ytmp = yleafloc;
    ztmp = zleafloc;
   
           for t = 1:length(idxfilter)
        gg = idxfilter(t);
        hold on
        surf(xx+ xtmp(gg), yy + ytmp(gg), zz + ztmp(gg), 'EdgeColor', 'r', 'FaceColor', 'r');
    end
    
    xtmp(idxfilter) = [];
    ytmp(idxfilter) = [];
    ztmp(idxfilter) = [];
    
 
    for g = 1:length(xtmp)
        %tt = out_node(g);
        surf(xx+ xtmp(g), yy + ytmp(g), zz + ztmp(g), 'EdgeColor', color, 'FaceColor', color);
        
    end
    
    
%     for g = 1:length(out_node)
%         tt = out_node(g);
%         surf(xx+ newntip(tt,1), yy + newntip(tt,2), zz + newntip(tt,3), 'EdgeColor', color, 'FaceColor', color);
%     end

    view(-140,20);axis equal;
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    xlim([-4,4]);
    ylim([-5.0,5]);
    zlim([0,8]);

    % Plot impulse response
    axes('Position',[.75, .75, .2, .2]);
    figure()
    plot((1:24000)./400E3*1E3,imp/max(imp),'-k');
    xlim([10,50])
    ylim([-1,1])
    xlabel('time (ms)','fontsize',15);
    ylabel('Amplitude','fontsize',15);
    set(gca, 'YTick',-1:2:1, 'YTickLabel',{'-1';' 1'});
    set(gcf, 'Renderer', 'zbuffer');
    % M(j) = getframe(gcf);
    y(j,:) = imp;
    % clf;
end

% save('ginkgo_7_twotrees_beamexpand.mat', 'y','lfnum', 'r', 'boxcenter', 'sonar_loc', 'bwidth', 'Rb', 'lfnode','nReps');
% movie2avi(M,'ginkgo7twotrees_beamexpand.avi','fps',1,'Compression','None');
                