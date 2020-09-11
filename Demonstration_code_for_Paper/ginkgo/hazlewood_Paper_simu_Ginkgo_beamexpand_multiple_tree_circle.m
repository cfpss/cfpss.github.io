%%% Record echoes from L-system ginkgo trees
% Created by Chen Ming, Date: 01/01/2017
clc;
clear all;
close all;
addpath('ginkgoLsys');
addpath('simuroutines');

leafsize = 2.5E-2; % in meter



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

IPP = IPP_points();
% IPP(2,1) = IPP(2,1) -2;
%    IPP(1,1) = IPP(1,1) + 2;

figure

newntip = [];
secondntip = [];
newBR = [];
lfnum = zeros(length(bwidth),1);
etas = [];
for j = 1:length(bwidth)
    
    for unige = 1:length(IPP(:,1))
        
        %figure
        j
        freqorbmwdthVal = bwidth(j);
        %     sonar_loc = [10  ys zs];
        %     boxcenter = [0, 0, 4.75];
        
        [BR, lf, ntip, Itip] = Paper_multiple_GMT3(nReps, bRatio, angle, Rb, lfnode);
        etas = [etas;size(BR,1)];
        %%% The second set of branches (tree)
        SecondBR = zeros(size(BR,1), size(BR,2));
        
        %%% Plot branches of the "two" trees
        RotM = rotz(90);
        N = size(BR,1);
        
        % for branch I did line 75 to 85
        
        for n = 1:N
            hold on
            endrotated = RotM*BR(n, 2:4)';
            endrotated(3) = endrotated(3) + Raa*Rab*Rb;
            tiprotated = RotM*BR(n, 5:7)';
            tiprotated(3) = tiprotated(3) + Raa*Rab*Rb;
            %cylinder2P(BR(n, 8), 10, BR(n, 2:4), BR(n, 5:7));
            %cylinder2P(BR(n, 8), 10, endrotated', tiprotated');
            SecondBR(n,:) = [BR(n,1) endrotated' tiprotated' BR(n,8)];
            
        end
        newbr = [BR; SecondBR];  % new branches coordinates from two trees
        % newItip = [Itip; Itip+length(Itip)/4];
        
        %%% Acquire echoes from the leaves of the "two" trees
        % figure
        secondNtip = RotM*transpose(ntip);
        secondNtip = secondNtip';
        secondNtip(:,3) = secondNtip(:,3) + Raa*Rab*Rb;
        newNtip = [ntip; secondNtip];
        treeLeaf2_loc = RotM*lf(:,1:3)';
        treeLeaf2_norm = RotM*lf(:,4:6)';
        treeLeaf2_loc_t = treeLeaf2_loc';
        treeLeaf2_loc_t(:,3) = treeLeaf2_loc_t(:,3) + Raa*Rab*Rb;

        
    
        newntip = [newntip;newNtip];
        secondntip = [secondntip;secondNtip];
        newBR = [newBR;newbr];
        hold on;
        
    end
    newLeaf = hazlewood_Paper_modified_2(newBR,IPP,etas);
    
    
    xn = 0;
    yn = 0;
    
    for i = 1:length(IPP(:,1))
        if i == 4
            continue
        end
        xn = xn + IPP(i,1);
        yn = yn + IPP(i,2);
        
    end
    
    xn = xn / (length(IPP(:,1)) -1);
    yn = yn / (length(IPP(:,2)) -1);
    
    d = 0;
    for i = 1:length(IPP(:,1))
        if i == 4
            continue
        end
        d1 = sqrt((xn-IPP(i,1))^2 +(yn-IPP(i,2))^2);
        
        if d1 > d
            x1 = IPP(i,1);
            y1 = IPP(i,2);
            d = d1;
        end
    end
    radius = 2 + sqrt((x1 - xn)^2 + (y1 - yn)^2);
    tic
    
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
    

    
    zleafloc = newLeaf(:,3);
    zmin = min(zleafloc);
    zmax = max(zleafloc);
    
 
    
    yleafloc = newLeaf(:,2);
    ymin = min(yleafloc);
    ymax = max(yleafloc);
    
 
    xs = 0;
    ys = 0;
    zs = 5;
    div = 6;
    
    idxf = {};
    impf = {};
    ctr = 1;
    
    Rad = sqrt(xs^2 + ys^2);
    st = 1;
    tx = 0;
    tx1 = 0;
    
    loop_counter = 1;
    
    for theta = 0:2*pi/div:2*pi-0.2    % for circle
        for ii = 1:2
            %bwidth1 = [5,7,10,13,20];
            bwidth1 = [50,20,25];
            freqorbmwdthVal = bwidth1(1+mod(st,length(bwidth1)));
            st = st+1;
           
            
            if ii == 2
                zs = 5;
                xs = IPP(4,1) + 6 * cos(theta);
                ys = IPP(4,2) + 6 * sin(theta);
                %if xs-16 < 1 && ys - 13 < 1
                %  continue
                % end
                
            else
                zs = 5;%25;
                xs = xn + radius * cos(theta);
                ys = yn + radius * sin(theta);
            end
            %ys = theta;
            YS(ctr) = ys;
            XS(ctr) = xs;
            
            sonar_loc = [xs  ys zs];
            
            boxcenter = [(xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2];
            
            T1 = floor(2*sqrt((xmax-xs)^2+(ymax-ys)^2)/340*1000);
            
            n_sample(ctr) = 24000;%(400E3*T1)/1000;
            
            t = (0:n_sample(ctr)-1)/400E3*1E3;
            echo_param.Nbins = n_sample(ctr); % no of bins b/w 60-80 kHz
            
            Xlfnorm = newLeaf(:,4);
            Zlfnorm = newLeaf(:,6);
            Ylfnorm = newLeaf(:,5);
            n_leaf = size(newLeaf,1);
            leaf_dia = lfsize_extract(leaf_param, n_leaf);
            
            R_leaf_norm = ones(n_leaf,1);
            [TH_leaf_norm, PHI_leaf_norm, ~] = cart2sph(Xlfnorm, Ylfnorm, Zlfnorm);
            
             hold on, waterdrop_sonarbeam(sonar_loc, boxcenter, freqorbmwdthVal);
            
            [Sphcoord_sonar, Sphcoord_leaf, beamcenter] = transcor(sonar_loc, boxcenter, xleafloc, yleafloc, zleafloc);
            [idxfilter, ~,~,R_sonar_leaf_filt,leaf_dia_filt,~,Sonargain_filt,Incident_angles] = filter_leaves_beamG(echo_param,Sphcoord_sonar,Sphcoord_leaf, TH_leaf_norm, PHI_leaf_norm, R_leaf_norm,leaf_dia,peakAmp, freqorbmwdth,beam_rot,freqorbmwdthVal , beamcenter);
            lfnum(j) = length(idxfilter);
            imp = time_domain_impulse1(Sonargain_filt,R_sonar_leaf_filt,leaf_dia_filt,Incident_angles,echo_param);
            
            idxf{ctr} = idxfilter;
            %    impf{ctr} = imp;
            ctr = ctr + 1;
            
            %end
            
            toc
            
            %%% Plot leaves that lied in the beam
            
            % Prepare the coordinates of a sphere to represent a cluster (bundle)
            % of ginkgo leaves
            [xx,yy,zz] = sphere(10);
            xx = 6E-2*xx;
            yy = 6E-2*yy;
            zz = 6E-2*zz;
            
            %for vt = 1:length(idxf)
            % idxfilter = idxf{vt};
            %  sonar_loc = [XS(vt) YS(vt) zs];
            % in_node = newItip(idxfilter);
            % in_node = unique(in_node);
            %    hold on, waterdrop_sonarbeam(sonar_loc, boxcenter, freqorbmwdthVal);
            
            
            %        for t = 1:length(in_node)
            %         gg = in_node(t);
            %         hold on
            %         %surf(xx+ newntip(gg,1), yy + newntip(gg,2), zz + newntip(gg,3), 'EdgeColor', 'r', 'FaceColor', 'r');
            %     end
            
            %red ones in the tree plot
            for t = 1:3:length(idxfilter)
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
            
            if loop_counter == 1
                for g = 1:11:length(xtmp)
                    % %   tt = out_node(g);
                    % surf(xx+ xtmp(g), yy + ytmp(g), zz + ztmp(g), 'EdgeColor', color, 'FaceColor', color);
                    
                end
                loop_counter = 0;
            end
            
            %      for g = 1:5:length(out_node)
            %          tt = out_node(g);
            %          surf(xx+ newntip(tt,1), yy + newntip(tt,2), zz + newntip(tt,3), 'EdgeColor', color, 'FaceColor', color);
            %             end
            
            %  view(2);
            %view(-140,20);
            %view(-164,30);
            %axis equal;
            xlabel('x (m)');
            ylabel('y (m)');
            zlabel('z (m)');
            %     xlim([-2,2]);
            %     zlim([0.02,9.5]);
            %     ylim([-2.2,2.2]);
            
            xlim([-6,30]);
            zlim([0,zmax]);
            ylim([-10,30]);
            
            %   xlim([xmin,xmax]);
            %     zlim([0,zmax]);
            %     ylim([ymin,ymax]);
                      
            
%             
%             figure()
%             
%             plot((1:24000)./400E3*1E3,imp/max(imp),'-k');
%             xlim([0,35])
%             ylim([-1,1])
%             xlabel('time (ms)','fontsize',15);
%             ylabel('Amplitude','fontsize',15);
%             set(gca, 'YTick',-1:2:1, 'YTickLabel',{'-1';' 1'});
%             set(gcf, 'Renderer', 'zbuffer');
%             y(j,:) = imp;
%             drawnow
        end
    end
end


