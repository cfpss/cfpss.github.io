%%% Record echoes from L-system ginkgo trees
% Created by Chen Ming, Date: 01/01/2017
clc;
clear all;
close all;
addpath('ginkgoLsys');
addpath('simuroutines');

leafsize = 2.5E-2; % in meter


%t = (0:24E3-1)/400E3*1E3;
%t = (0:12324-1)/400E3*1E3;

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
%echo_param.Nbins = 24E3; % no of bins b/w 60-80 kHz



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
bwidth = 10;%10:20:50; % sonar beamwidth changing from 10 degrees to 50 degrees
%y = zeros(length(bwidth),12324); % echoes


%ys = -10;
%zs = 32;
%ys = 126;
%zs = 168;
lfnum = zeros(length(bwidth),1);

for j = 1:length(bwidth)
    figure
    j
    freqorbmwdthVal = bwidth(j);
%     sonar_loc = [10  ys zs];
%     boxcenter = [0, 0, 4.75];

    [BR, lf, ntip, Itip] = GMT3(nReps, bRatio, angle, Rb, lfnode);

    %%% The second set of branches (tree)
    SecondBR = zeros(size(BR,1), size(BR,2));
    
    %%% Plot branches of the "two" trees
    RotM = rotz(90);
    N = size(BR,1);
    
   % for brach I did line 75 to 85
    
        for n = 1:N
            hold on
            endrotated = RotM*BR(n, 2:4)';
            endrotated(3) = endrotated(3) + Raa*Rab*Rb; 
            tiprotated = RotM*BR(n, 5:7)';
            tiprotated(3) = tiprotated(3) + Raa*Rab*Rb;
            cylinder2P(BR(n, 8), 10, BR(n, 2:4), BR(n, 5:7));
            cylinder2P(BR(n, 8), 10, endrotated', tiprotated');   
            SecondBR(n,:) = [BR(n,1) endrotated' tiprotated' BR(n,8)];
            
        end
    newBR = [BR; SecondBR];  % new branches coordinates from two trees
    % newItip = [Itip; Itip+length(Itip)/4];

    %%% Acquire echoes from the leaves of the "two" trees
    figure
    secondntip = RotM*transpose(ntip);
    secondntip = secondntip';
    secondntip(:,3) = secondntip(:,3) + Raa*Rab*Rb;
    newntip = [ntip; secondntip];
    treeLeaf2_loc = RotM*lf(:,1:3)';
    treeLeaf2_norm = RotM*lf(:,4:6)';
    treeLeaf2_loc_t = treeLeaf2_loc';
    treeLeaf2_loc_t(:,3) = treeLeaf2_loc_t(:,3) + Raa*Rab*Rb;
    %newLeaf = [lf; [treeLeaf2_loc_t treeLeaf2_norm']];

%newLeaf = modified(); 
%newLeaf = modified_1(); 
%newLeaf = modified_2(newBR); 
%newLeaf = modified_3(newBR); 
newLeaf = Paper_modified_3_A(newBR); 
%newLeaf = newLeaf(1:100,:);

ntip = newLeaf(:,1:3);
secondntip = RotM*transpose(ntip);
    secondntip = secondntip';
    secondntip(:,3) = secondntip(:,3) + Raa*Rab*Rb;
    newntip = [ntip; secondntip];
    
    
newItip = [];
counter = 1;
for i = 1:4:length(newLeaf)
    for j = 0:3
        newItip(i+j) = counter;
    end
counter = counter + 1;    
end
newItip = newItip'

%xleafloc = newLeaf(:,1)*1E-2;
    xleafloc = newLeaf(:,1);
    xmin = min(xleafloc);
    xmax = max(xleafloc);
    
%     for i=1:length(xleafloc)
%     xleafloc(i) = xleafloc(i)*rand(1);
%       end
    
    zleafloc = newLeaf(:,3);
    zmin = min(zleafloc);
    zmax = max(zleafloc);
    
%       for i=1:length(zleafloc)
%     zleafloc(i) = zleafloc(i)*rand(1);
%       end
%       
    
    yleafloc = newLeaf(:,2);
    ymin = min(yleafloc);
    ymax = max(yleafloc);
    
%         for i=1:length(yleafloc)
%     yleafloc(i) = yleafloc(i)*rand(1);
%       end

%for sonar = 0:2:10
    xs = 2;
    ys = 6; 
    zs = 15; %31
    
    sonar_loc = [xs  ys zs];
    
    boxcenter = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];  
    
    T1 = floor(2*sqrt((xmax-xs)^2+(ymax-ys)^2)/340*1000);
    
    n_sample = (400E3*T1)/1000;

     t = (0:n_sample-1)/400E3*1E3;
    echo_param.Nbins = n_sample; % no of bins b/w 60-80 kHz
    
    Xlfnorm = newLeaf(:,4);
    Zlfnorm = newLeaf(:,6);
    Ylfnorm = newLeaf(:,5);
    n_leaf = size(newLeaf,1);
    leaf_dia = lfsize_extract(leaf_param, n_leaf);

    R_leaf_norm = ones(n_leaf,1);
    [TH_leaf_norm, PHI_leaf_norm, ~] = cart2sph(Xlfnorm, Ylfnorm, Zlfnorm);

    hold on
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

       for t = 1:length(in_node)
        gg = in_node(t);
        hold on
        surf(xx+ newntip(gg,1), yy + newntip(gg,2), zz + newntip(gg,3), 'EdgeColor', 'r', 'FaceColor', 'r');
    end
    
 %red ones in the tree plot   
    for t = 1:50:length(idxfilter)
        gg = idxfilter(t);
        hold on
        surf(xx+ xleafloc(gg), yy + yleafloc(gg), zz + zleafloc(gg), 'EdgeColor', 'r', 'FaceColor', 'r');
    end

    %%% Plot leaves that were not illuminated by the beam
    newItip(idxfilter,:) = [];
    out_node = unique(newItip);
    color = [83/255, 132/255, 10/255];

    %for greencolor   
    xtmp = xleafloc;
    ytmp = yleafloc;
    ztmp = zleafloc;
    xtmp(idxfilter) = [];
    ytmp(idxfilter) = [];
    ztmp(idxfilter) = [];
    
    for g = 1:50:length(xtmp)
        %tt = out_node(g);
        surf(xx+ xtmp(g), yy + ytmp(g), zz + ztmp(g), 'EdgeColor', color, 'FaceColor', color);
        
    end
    
%     for g = 1:length(out_node)
%         tt = out_node(g);
%         surf(xx+ newntip(tt,1), yy + newntip(tt,2), zz + newntip(tt,3), 'EdgeColor', color, 'FaceColor', color);
%     end

    view(-140,20);
    axis equal;
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
%     xlim([-2,2]);
%     zlim([0.02,9.5]);
%     ylim([-2.2,2.2]);

  xlim([xmin,xmax]);
    zlim([zmin,zmax]);
    ylim([ymin,ymax]);
    axis equal;
    %set(gca,'FontSize',14);
drawnow

%end
%     hold on;
% [V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/47-mapletree/MapleTreeStem.obj');
% 
% trisurf(F,V(:,1)*0.25,V(:,3)*0.25,V(:,2)*0.25,'FaceColor',[0.26,1,0.33 ]);
% 
% light('Position',[-1.0,-1.0,100.0],'Style','infinite');
% lighting phong;
% hold off;

    % Plot impulse response
    axes('Position',[.75, .75, .2, .2]);
    
    figure()

    %plot((1:24000)./400E3*1E3,imp/max(imp),'-k');
    plot((1:n_sample)./400E3*1E3,imp/max(imp),'-k');
    xlim([0,35])
    ylim([-1,1])
    xlabel('time (ms)','fontsize',15);
    ylabel('Amplitude','fontsize',15);
    set(gca, 'YTick',-1:2:1, 'YTickLabel',{'-1';' 1'});
    set(gcf, 'Renderer', 'zbuffer');
    % M(j) = getframe(gcf);
    y(j,:) = imp;
    drawnow
    % clf;
end


% hold on
% [X,Y,Z] = cylinder(1)
% surf(X,Y,Z)

% save('ginkgo_7_twotrees_beamexpand.mat', 'y','lfnum', 'r', 'boxcenter', 'sonar_loc', 'bwidth', 'Rb', 'lfnode','nReps');
% movie2avi(M,'ginkgo7twotrees_beamexpand.avi','fps',1,'Compression','None');
                
