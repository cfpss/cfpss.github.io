% This script is modified based on Chen's code, just to plot the tree
% branching and tree with leaves. 
clc;clear;
addpath('ginkgoLsys');
addpath('simuroutines');

leafsize = 2.5E-2; % in meter
%%% Parameters you don't change
leaf_param.leafloc_dist = 'uni';
leaf_param.leafsize_dist = 'gauss';
leaf_param.leafsize_supp = leafsize;

%%% Parameters used in building L-system tree
nReps = 4; % growth level
bRatio = 0.62; % contraction ratio of branches other than center branch
angle = 50; % branching angle
Rb = 2; % the length of initial branch
lfnode = 4; % the number of leaf nodes per branch
Raa = 0.5; % lift the copy half length of the initial branch
Rab = 0.826; % Contraction ratio of center branch

zs = 6;
boxcenter = [0, 0, zs];
% Generate branching structure
[BR, lf, ntip, Itip] = GMT3(nReps, bRatio, angle, Rb, lfnode);

% Plot the branching structure
figure()
SecondBR=zhu_branch_rotate_plot(BR, Raa, Rab, Rb);

newBR = [BR; SecondBR];  % new branches coordinates from two trees
% Get the new leaf parameters.
RotM = rotz(90);
newItip = [Itip; Itip+length(Itip)/4];
secondntip = RotM*transpose(ntip);
secondntip = secondntip';
secondntip(:,3) = secondntip(:,3) + Raa*Rab*Rb;
newntip = [ntip; secondntip];
treeLeaf2_loc = RotM*lf(:,1:3)';
treeLeaf2_norm = RotM*lf(:,4:6)';
treeLeaf2_loc_t = treeLeaf2_loc';
treeLeaf2_loc_t(:,3) = treeLeaf2_loc_t(:,3) + Raa*Rab*Rb;
newLeaf = [lf; [treeLeaf2_loc_t treeLeaf2_norm']];

xleafloc = newLeaf(:,1);
yleafloc = newLeaf(:,2);
zleafloc = newLeaf(:,3);
Xlfnorm = newLeaf(:,4);
Ylfnorm = newLeaf(:,5);
Zlfnorm = newLeaf(:,6);
n_leaf = size(newLeaf,1);
%leaf_dia = lfsize_extract(leaf_param, n_leaf);
leaf_dia = zhu_lfsize_extract(leaf_param, n_leaf, 0.5, 0.01, 0.2);
figure()
hist(leaf_dia)
% Prepare the coordinates of a sphere to represent a cluster (bundle)
% of ginkgo leaves

% Note: In this code, Chen Ming used the identifical balls as trees. No 
% scale according to the leaf_dia meter and no effort of making it like a
% disk. 
[xx,yy,zz] = sphere(10);
xx = 5E-2*xx;
yy = 5E-2*yy;
zz = 5E-2*zz;
figure()
surf(xx, yy, zz)

all_node = unique(newItip);
color = [83/255, 132/255, 10/255];
for g = 10:length(all_node)
    tt = all_node(g);
    surf(xx + newntip(tt,1), yy + newntip(tt,2), zz + newntip(tt,3), 'EdgeColor', color, 'FaceColor', color);
end

shading interp; colormap jet(256); camlight;

%     view(-140,20);axis equal;
%     xlabel('x (m)');
%     ylabel('y (m)');
%     zlabel('z (m)');
%     xlim([-4,4]);
%     ylim([-5,6]);
%     zlim([0,10]);

    % Plot impulse response
    %axes('Position',[.75, .75, .2, .2]);
    %plot((1:24000)./400E3*1E3,imp/max(imp));
    %xlim([0,60])
    %ylim([-1,1])
    %xlabel('time (ms)');
    %ylabel('amplitude (norm.)');
    %set(gca, 'YTick',-1:2:1, 'YTickLabel',{'-1';' 1'});
    %set(gcf, 'Renderer', 'zbuffer');
    % M(j) = getframe(gcf);
    %y(j,:) = imp;
    % clf;
%end

% save('ginkgo_7_twotrees_beamexpand.mat', 'y','lfnum', 'r', 'boxcenter', 'sonar_loc', 'bwidth', 'Rb', 'lfnode','nReps');
% movie2avi(M,'ginkgo7twotrees_beamexpand.avi','fps',1,'Compression','None');
                