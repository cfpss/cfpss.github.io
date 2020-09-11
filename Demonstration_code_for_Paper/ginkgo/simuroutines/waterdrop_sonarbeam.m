% plot beampattern with L-system trees
% Created by Chen Ming, Date: 11/01/2016
function waterdrop_sonarbeam(sonar_loc, boxcenter, freqorbmwdthVal, n_sample, split_1)

xs = sonar_loc(1);
ys = sonar_loc(2);
zs = sonar_loc(3);

Vecbeam = boxcenter - sonar_loc;
[TH_beam, PHI_beam, ~] = cart2sph(Vecbeam(1), Vecbeam(2), Vecbeam(3));

N = 2000;
freqorbmwdth = 'beamwidth';
%freqorbmwdthVal = 10;% beamwidth in degrees (elevation and azimuth same)
sigma_x = predict_gauss_sigma(freqorbmwdth,freqorbmwdthVal); % get sigma for beamwidth corresponding to frequency of incidence
sigma_y = predict_gauss_sigma(freqorbmwdth,freqorbmwdthVal); % get sigma for beamwidth corresponding to frequency of incidence
x0 = TH_beam;
y0 = PHI_beam;
[az, el] = meshgrid(linspace(x0-pi/2,x0+pi/2,N), linspace(y0-pi/2,y0+pi/2,N));

theta = 0;
A = 1;

hold on
dB = 20;
r0 = 0.1;

a = ((cos(theta))^2)/(2*(sigma_x^2)) + ((sin(theta))^2)/(2*(sigma_y)^2);
b = (-sin(2*theta))/(4*(sigma_x)^2) + (sin(2*theta))/(4*(sigma_y)^2);
c = ((sin(theta))^2)/(2*(sigma_x)^2) + ((cos(theta))^2)/(2*(sigma_y)^2);
    


   %Z = A*exp( -((a*(X-x0).^2) + (2*b*(X-x0).*(Y-y0)) + (c*(Y-y0).^2)));
Z = A*exp(-((a*(az-x0).^2) + (2*b*(az-x0).*(el-y0)) + (c*(el-y0).^2)));

%sonarbp = exp(-( az - az_mean ).^2./2./gauss_sigma^2);
RX_square = 10^(-dB/20)*r0^2*1E5*(Z).^2;%
RX = sqrt(RX_square);

%idxfilter = find(R_sonar_leaf <= dis_b); %find the index where gain is greater than threshold

%%% Extract leaves & parameters that pass gain threshld
% TH_leaf_filt = TH_leaf(idxfilter);
% PHI_leaf_filt = PHI_leaf(idxfilter);
% R_leaf_filt = R_leaf(idxfilter);





x = xs + RX.*cos(el).*cos(az);
y = ys + RX.*cos(el).*sin(az);
z = zs + RX.*sin(el);
hSurface = surf(x,y,z,'EdgeColor','none','FaceLighting','phong');
set(hSurface,'FaceColor','r','FaceAlpha',0.2);
%camlight; lighting gouraud;
hold on, scatter3(xs, ys, zs, 'filled');



dist_max = 0;

for i = 1:length(x(:,1))
    for j = 1:length(x(1,:))
        dist = sqrt((xs -x(i,j))^2 + (ys -y(i,j))^2 + (zs -z(i,j))^2);
        
        if dist > dist_max
            dist_max = dist;
            row = i;
            column = j;
        end
    end
end

          


%% to divide beamwidth into different sections
%plot3([xs x(row,column)], [ys y(row,column)], [zs z(row,column)], 'LineWidth', 2);   % line connecting sonar and ellipsoid end
%theta = atan2(ys, xs);                    
division = 0.5*340* (n_sample / split_1 / 400E3);
dist_real = division * split_1;

%ellipsoid_dis = sqrt( (xs - x(row,column))^2 + (ys- y(row,column))^2 +(zs- z(row,column))^2 );
%divisions = ellipsoid_dis/split_1;

%m = [x(row,column) y(row,column) z(row,column)] - [xs ys zs];
%m = m * dist_real/dist_max;


% for i = 1 : split_1
%    
%    % get x and y alng sonar ellipsoid line for divisions
%    x_new = xs + i * 1/split_1 * m(1); %cos(pi+theta);
%    y_new = ys + i * 1/split_1 * m(2); %sin(pi+theta);
%    z_new = zs + i * 1/split_1 * m(3); %sin(pi+theta);
%    %plot3(x_new, y_new, z_new, 'Marker', 'o','MarkerSize' , 5, 'MarkerEdgeColor', [0 1 0],...
%           %   'MarkerFaceColor' , [0 1 0], 'color',  'g');   % green dot
%        r1 = sqrt((x_new-xs)^2 + (y_new-ys)^2 + (z_new-zs)^2)
% %          w_max = max(az(N*i/split_1,:));
% %          w_min = min(az(N*i/split_1,:));
% %          
% %         plot(xs+2*i* cos(w_min+pi/2:0.1:w_max),ys+ 2*i*sin(w_min+pi/2:0.1:w_max), 'LineWidth', 2);
% 
%          if i < 2
%             % for z_i = zs -1: 0.1: zs+1
%          %  plot3(xs+r1* cos(x0-pi/(2.5):0.1:x0+pi/(2.5)),ys+ r1*sin(x0-pi/(2.5):0.1:x0+pi/(2.5)), zs * ones(1,1+floor(2*pi/(2.5*0.1))), 'LineWidth', 2);
%            %  end
%              
%          elseif i > 1 && i < split_1
%                  
%   %  plot3(xs+r1* cos(x0-pi/(3*i/2):0.1:x0+pi/(3*i/2)),ys+ r1*sin(x0-pi/(3*i/2):0.1:x0+pi/(3*i/2)), zs * ones(1,1+floor(2*pi/((3*i/2)*0.1))), 'LineWidth', 2);
%          end
% %          
%   
%   % [xx, yy] = meshgrid(x_new-1: 0.1: x_new+1, y_new-1: 0.1: y_new+1)      
%   % zz = (-m(1) * xx - m(2) * yy + x_new * m(1) + y_new * m(2) + z_new * m(3)) / m(3);
%   % surf(xx,yy,zz);
%   % hSurface1 = surf(xx,yy,zz,'EdgeColor','none','FaceLighting','phong');
% %set(hSurface1,'FaceColor','r','FaceAlpha',0.5);
%    
%    
% end
end

