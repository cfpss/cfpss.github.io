function total_leaf = hazlewood_Paper_modified_2(newBR,IPP,etas)
%obj = readObj('/home/hassan/Documents/obj_test/Fir_OBJ/Fir.obj');   %use appropriate file
%obj = readObj('/home/hassan/Documents/tree2/Tree.obj');
obj_1 = readObj('/home/hassantanveer/obj_test/hazelnutbush/Hazelnut.obj');
%shading interp
%colormap(gray(256));
%lighting phong
%camproj('perspective');
%axis square
%axis off
%axis equal
%axis tight
%cameramenu

total_leaf = [];
face = [];
obj_1.f.v = obj_1.f.v';
obj_1.f.vn = obj_1.f.vn';
%set(gcf,'renderer','opengl')
%set(0,'DefaultFigureRenderer','painters');

 %%%%%%%%%%%%%%%%%%%%%%%%%%%% trunk %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:37700
        v1 = obj_1.f.v{i}(1);
        v2 = obj_1.f.v{i}(2);
        v3 = obj_1.f.v{i}(3);
        v4 = obj_1.f.v{i}(4);
       % if (v4>24498)
         %   continue;
       % end
        face = [face;[v1,v2,v3]];
        
        x1 = (obj_1.v(v1,1) + obj_1.v(v2,1) + obj_1.v(v3,1))/3;
        y1 = (obj_1.v(v1,2) + obj_1.v(v2,2) + obj_1.v(v3,2))/3;
        z1 = (obj_1.v(v1,3) + obj_1.v(v2,3) + obj_1.v(v3,3))/3;
        
        plot3(x1,z1,y1,' ');
        hold on;
    end
    trunk = face(1:37700,:);
    
    % IPP(2,1)
    %for ipp = 1:length(IPP(:,1))
    tmp_1 = (obj_1.v(:,1))/1.5 + IPP(4,1);
    tmp_2 = (obj_1.v(:,3))/1.5+ IPP(4,2);
    figure
    trisurf(trunk,tmp_1,tmp_2,obj_1.v(:,2),'FaceColor',[0.2, 0, 0]);
    grid off;
    hold on;


%%%%%%leaves%%%%%%%%%%%%%%%%%%%%
j =1;
leafs = [];
leafs_norm = [];
for L = 40589:42522
    
             %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 2.5*rand(1)/3,obj_1.v(L,3)/1.5+ IPP(4,2)+ 2.5*rand(1)/3,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%% %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 2.8*rand(1)/3,obj_1.v(L,3)/1.5+ IPP(4,2)+ 2.8*rand(1)/3,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%%%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 2*rand(1)/3,obj_1.v(L,3)/1.5+ IPP(4,2)+ 2*rand(1)/3,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
        %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 5*rand(1)/6,obj_1.v(L,3)/1.5+ IPP(4,2)+ 5*rand(1)/6,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
    %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 6*rand(1)/7,obj_1.v(L,3)/1.5+ IPP(4,2)+ 6*rand(1)/7,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
      %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 5*rand(1)/7,obj_1.v(L,3)/1.5+ IPP(4,2)+ 5*rand(1)/7,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
    %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 2*rand(1)/8,obj_1.v(L,3)/1.5+ IPP(4,2)+ 2*rand(1)/8,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
      %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 8*rand(1)/9,obj_1.v(L,3)/1.5+ IPP(4,2)+ 8*rand(1)/9,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
    
      %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 9*rand(1)/8,obj_1.v(L,3)/1.5+ IPP(4,2)+ 9*rand(1)/8,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
      %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 2*rand(1)/7,obj_1.v(L,3)/1.5+ IPP(4,2)+ 2*rand(1)/7,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
    
      %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + rand(1)/2,obj_1.v(L,3)/1.5+ IPP(4,2)+ rand(1)/2,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
    
             leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 3*rand(1)/4,obj_1.v(L,3)/1.5+ IPP(4,2)+ 3*rand(1)/4,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
     leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 2*rand(1)/5,obj_1.v(L,3)/1.5+ IPP(4,2)+ 2*rand(1)/5,obj_1.v(L,2)];
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
    leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 4*rand(1)/7,obj_1.v(L,3)/1.5+ IPP(4,2)+ 4*rand(1)/7,obj_1.v(L,2)];
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%ok%%%%
             leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 3*rand(1)/2,obj_1.v(L,3)/1.5+ IPP(4,2)+ 3*rand(1)/2,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            %%%%ok%%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1),obj_1.v(L,3)/1.5+ IPP(4,2),obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%
            leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + rand(1)/2,obj_1.v(L,3)/1.5+ IPP(4,2)+ rand(1)/2,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            %%%%ok%%%%%
                leafs(j,:) = [obj_1.v(L,1)/1.5+ IPP(4,1) + 3*rand(1)/2,obj_1.v(L,3)/1.5+ IPP(4,2)+ 3*rand(1)/2,obj_1.v(L,2)];
            
            leafs_norm (j,:) = obj_1.vn(L,:);
            total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
            j = j+1;
            
            end

% 
% 
% 
% %%%%%%%%%%%next tree%%%%%%%%%%%%%%%%%%
%

obj = readObj('/home/hassantanveer/47-mapletree/MapleTree.obj');
%obj = readObj('/home/hassan/Documents/47-mapletree/MapleTreeStem.obj');

obj.v(:,1) = obj.v(:,1);
obj.v(:,2) = obj.v(:,2);
obj.v(:,3) = obj.v(:,3);
ns = 2;

% for FIR Tree
obj.f.v = obj.f.v';
obj.f.vn = obj.f.vn';
%set(gcf,'renderer','opengl')
%set(0,'DefaultFigureRenderer','painters');


pre_eta = 0;
NewBR = newBR;

IPP = [IPP(1:3,:); IPP(5,:)];
%IPP = [IPP(1:4,:)];
hold on;
for ipp = 1:length(IPP(:,1))
    %leafs = [];
    newangle = [];
  
    lsys_branch = [];
    eta = etas(ipp);
    newBR = NewBR(2 * pre_eta + 1 : 2 * pre_eta + 2*eta, :);
    pre_eta = pre_eta + eta;
    face = [];
    
    %figure
    
    
    %for i=1:length(obj.f.v)
    
    %for i=1:length(obj.f.v)-10000
    
    rand_leaf = [];
    
    y = 1;
    
for lsys = 1:length(newBR)
    if newBR(lsys,4) < 1
        continue;
    end
    
    if (newBR(lsys,2) == 0 && newBR(lsys,3) == 0 && newBR(lsys,5) == 0 && newBR(lsys,6) == 0)
        lsys_branch(y,:) = newBR(lsys,:);
        y = y+1;
    end
end

kk = 1;
for p = 1:length(lsys_branch(:,1))
    
    for k = 1:length(newBR)
        if newBR(k,4) == lsys_branch(p,4)
            if (newBR(k,2) == 0 && newBR(k,3) == 0 && newBR(k,5) == 0 && newBR(k,6) == 0)
                continue
            else
                newangle(kk) = atan2(newBR(k,6),newBR(k,5));
                kk = kk+1;
                
                newangle(kk) = atan2(newBR(k,6),newBR(k,5)) + 3* 43.75*pi/180;
                kk = kk+1;
            end
        end
    end
    
end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% trunk %%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:100
    v1 = obj.f.v{i}(1);
    v2 = obj.f.v{i}(2);
    v3 = obj.f.v{i}(3);
    v4 = obj.f.v{i}(4);
    if (v4>24498)
        continue;
    end
    face = [face;[v1,v2,v3]];
    
    x1 = (obj.v(v1,1) + obj.v(v2,1) + obj.v(v3,1))/3;
    y1 = (obj.v(v1,2) + obj.v(v2,2) + obj.v(v3,2))/3;
    z1 = (obj.v(v1,3) + obj.v(v2,3) + obj.v(v3,3))/3;
    
    %plot3(x1,z1,y1,' ');
    hold on;
end
trunk = face(1:100,:);

    tmp_1 = (obj.v(:,1)*0.5) + IPP(ipp,1);
    tmp_2 = (obj.v(:,3)*0.5) + IPP(ipp,2);

    trisurf(trunk,tmp_1,tmp_2,obj.v(:,2),'FaceColor',[0.6350, 0.0780, 0.18400]);


%%%%%%%%%%%%%%%%%%%%%%%%%%% Branch %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i=3073:3863
face = [];
for i=193:312
    v1 = obj.f.v{i}(1);
    v2 = obj.f.v{i}(2);
    v3 = obj.f.v{i}(3);
    v4 = obj.f.v{i}(4);
    if (v4>24498)
        continue;
    end
    face = [face;[v1,v2,v3]];
    
    x1 = (obj.v(v1,1) + obj.v(v2,1) + obj.v(v3,1))/3;
    y1 = (obj.v(v1,2) + obj.v(v2,2) + obj.v(v3,2))/3;
    z1 = (obj.v(v1,3) + obj.v(v2,3) + obj.v(v3,3))/3;
    
    plot3(x1,z1,y1,' ');
    hold on;
end

tmpz = obj.v(:,2) -24;

i = 1;

previous_height = lsys_branch(1,4);
tmp_i = 1;
x1 = obj.v(:,1);
x2 = obj.v(:,3);

rotation = 0;
rand_rot = [];
count1 = 1;
kk = 1;
while i <= length(lsys_branch(:,1))
    rand_b = 1;
    count = 1;
    
    
    for j = 201:328
        tmpz(j) = tmpz(j) + (j-201)*5*rand_b*10E-03;
        
        if (j== 221 || j== 234 || j== 244 || j== 251 || j== 260 || j== 273 ...
                || j== 281 || j== 291 || j== 305 || j== 314 || j== 322)
            
            Sb_root{i,count} = tmpz(j);
            count = count + 1;
            
        end
        
    end
    
    rand_a = 1*pi/12;
    rand_rot = [rand_rot; rand_a];
    
    count1 = count1 + 1;
    
    rotation = newangle(kk)+rand_a;
    kk = kk+1;
        %for ipp = 1:length(IPP(:,1))
        tmp_1 = IPP(ipp,1);
        tmp_2 = IPP(ipp,2);
    
    trisurf(face,...
        (cos(rotation)*obj.v(:,1)+sin(rotation)*obj.v(:,3))*(1/(tmp_i*ns + 1)) + tmp_1,...
        (-sin(rotation)*obj.v(:,1)+cos(rotation)*obj.v(:,3))*(1/(tmp_i*ns + 1)) + tmp_2,...
       tmpz + lsys_branch(i,4)*3,...
        'FaceColor',[0.2,0,0 ]);
    
    rotation = newangle(kk)+rand_a;
    kk = kk+1;
    

    trisurf(face,...
        (cos(rotation)*obj.v(:,1)+sin(rotation)*obj.v(:,3))*(1/(tmp_i*ns + 1))+tmp_1,...
        (-sin(rotation)*obj.v(:,1)+cos(rotation)*obj.v(:,3))*(1/(tmp_i*ns + 1))+tmp_2,...
        tmpz+ 2 + lsys_branch(i,4)*3,...
        'FaceColor',[0.2,0,0 ]);
    
    
    i = i+1;
    
    if (i <= length(lsys_branch(:,1)))
        if (previous_height < lsys_branch(i,4))
            tmp_i = tmp_i + 1;
            previous_height = lsys_branch(i,4);
            
        else
            tmp_i = 1.5;
            previous_height = lsys_branch(i,4);
        end
    end
    
    tmpz = obj.v(:,2) - 24;
    
end


%%%%%%%%%%% Sub Branch%%%%%%%%%%%%%%
face = [];
for i=3073:3863
    v1 = obj.f.v{i}(1);
    v2 = obj.f.v{i}(2);
    v3 = obj.f.v{i}(3);
    v4 = obj.f.v{i}(4);
    if (v4>24498)
        continue;
    end
    face = [face;[v1,v2,v3]];
    
    
    x1 = (obj.v(v1,1) + obj.v(v2,1) + obj.v(v3,1))/3;
    y1 = (obj.v(v1,2) + obj.v(v2,2) + obj.v(v3,2))/3;
    z1 = (obj.v(v1,3) + obj.v(v2,3) + obj.v(v3,3))/3;
    
    plot3(x1,z1,y1,' ');
    hold on;
end


tmpz = obj.v(:,2) - 24;
tmpx = obj.v(:,1);
tmpy = obj.v(:,3);

kk = 1;
i = 1;

previous_height = lsys_branch(1,4);
tmp_i = 1;
x1 = obj.v(:,1);
x2 = obj.v(:,3);

rotation = 0;
count2 = 1;

sbl = [3353,3433,3513,3591,3674,3754,3833,3919,3990,4085];
while i <= length(lsys_branch(:,1))
    
    facetmp = [];
    count = 1;
    rand_SB = 1;
    
    
    skip = randi(length(sbl));
    
    if skip == length(sbl)
        tmp = [sbl(skip-1):sbl(end)];
    elseif skip == 1
        tmp = [sbl(skip):sbl(skip+2)-1];
    else
        tmp = [sbl(skip-1):sbl(skip+1)-1];
    end
    tmp = [tmp,[sbl(6+randi(4)):sbl(end)]];
    for pp = 1:length(face)
        sv = sum(ismember(face(pp,:),tmp));
        if ~sv
            facetmp = [facetmp;face(pp,:)];
        end
    end
    
    
    for j = 3273:4152
        
        if (j== 3273 || j== 3353 || j== 3433 || j== 3513 || j== 3591 || j== 3674 ...
                || j== 3754 || j== 3833 || j== 3919 || j== 3990 || j== 4085)
            
            tmpz(j) = Sb_root{i,count};
            count = count + 1;
            root = j;
            shf = tmpz(j);
        end
        tmpz(j) = shf + (j-root)*5*rand_SB*10E-03;
    end
    
    
    
    rand_a = rand_rot(count2);
    count2 = count2 + 1;
    
    
    rotation = newangle(kk)+rand_a;
    kk = kk+1;
    
     tmp_1 = IPP(ipp,1);
        tmp_2 = IPP(ipp,2);
    
     trisurf(facetmp,...
         (cos(rotation)*tmpx+sin(rotation)*tmpy)*(1/(tmp_i*ns + 1))+tmp_1,...
         (-sin(rotation)*tmpx+cos(rotation)*tmpy)*(1/(tmp_i*ns + 1))+tmp_2,...
         tmpz + lsys_branch(i,4)*3,...
         'FaceColor',[0.2,0,0 ]);     
    
    rotation = newangle(kk)+rand_a;
    kk = kk+1;
    
    trisurf(facetmp,...
        (cos(rotation)*tmpx+sin(rotation)*tmpy)*(1/(tmp_i*ns + 1))+tmp_1,...
        (-sin(rotation)*tmpx+cos(rotation)*tmpy)*(1/(tmp_i*ns + 1))+tmp_2,...
        tmpz +2+ lsys_branch(i,4)*3,...
        'FaceColor',[0.2,0,0 ]);
    
    
    i = i+1;
    
    if (i <= length(lsys_branch(:,1)))
        if (previous_height < lsys_branch(i,4))
            tmp_i = tmp_i + 1;
            previous_height = lsys_branch(i,4);
            
        else
            tmp_i = 1.5;
            previous_height = lsys_branch(i,4);
        end
    end
    tmpz = obj.v(:,2) - 24;
    tmpx = obj.v(:,1);
    tmpy = obj.v(:,3);
end



%%%%%%%%%%% Leaf %%%%%%%%%%%%%%%%
face = [];
j = 1;
kk = 1;

for i=18841:18956
    v1 = obj.f.v{i}(1);
    v2 = obj.f.v{i}(2);
    v3 = obj.f.v{i}(3);
    v4 = obj.f.v{i}(4);
    if (v4>24498)
        continue;
    end
    face = [face;[v1,v2,v3]];
    
    %x1 = (obj.v(v1,1) + obj.v(v2,1) + obj.v(v3,1))/3;
    x1 = obj.v(v1,1);
    z1 = obj.v(v1,3);
    y1 = obj.v(v1,2);
    
    x2 = obj.v(v2,1);
    z2 = obj.v(v2,3);
    y2 = obj.v(v2,2);
    
    x3 = obj.v(v3,1);
    z3 = obj.v(v3,3);
    y3 = obj.v(v3,2);
    
    %y1 = (obj.v(v1,2) + obj.v(v2,2) + obj.v(v3,2))/3;
    %z1 = (obj.v(v1,3) + obj.v(v2,3) + obj.v(v3,3))/3;
    
    plot3(x1,z1,y1,' ');
    hold on;
    
end

tmp_angle = 0;
tmpz = obj.v(:,2) - 24;
i = 1;

previous_height = lsys_branch(1,4);
tmp_i = 1;
x1 = obj.v(:,1);
x2 = obj.v(:,3);

rotation1 = 0;
k = 1;
count3 = 1;
while i <= length(lsys_branch(:,1))
    
    count = 1;
    
    rand_SB = 1;
    rand_leaf = [rand_leaf; rand_SB];
    
    
      rand_a = rand_rot(count3);

    count3 = count3 + 1;
    
        
    rotation1 = newangle(kk)+rand_a;
    kk = kk+1;
      rotation2 = newangle(kk)+rand_a;
    kk = kk+1;
    
    
    
    for L = 20793:21255
        if (L== 20793 || L== 20857 || L== 20893 || L== 20933 || L== 20973 || L== 21069 ...
                || L== 21117 || L== 21165 || L== 21185 || L== 21209 || L== 21253)
            
            tmpz(L) = Sb_root{i,count};
            count = count + 1;
            root = L;
            shf = tmpz(L);
        end
            tmpz(L) = shf +1+ (L-root)*9*rand_leaf(k)*10E-03;

    
    
     tmp_1 = IPP(ipp,1);
        tmp_2 = IPP(ipp,2); 
    

    
        
    leafs(j,:) = [(cos(rotation1)*obj.v(L,1)+sin(rotation1)*obj.v(L,3))*(1/(tmp_i*ns + 1))+tmp_1,...
        (-sin(rotation1)*obj.v(L,1)+cos(rotation1)*obj.v(L,3))*(1/(tmp_i*ns + 1))+tmp_2,...
        (tmpz(L)+lsys_branch(i,4)*3)];
    leafs_norm (j,:) = obj.vn(v4,:);
    total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
    j = j+1;
    

   
    leafs(j,:) = [(cos(rotation2)*obj.v(L,1)+sin(rotation2)*obj.v(L,3))*(1/(tmp_i*ns + 1))+tmp_1,...
        (-sin(rotation2)*obj.v(L,1)+cos(rotation2)*obj.v(L,3))*(1/(tmp_i*ns + 1))+tmp_2,...
        (tmpz(L)+2+lsys_branch(i,4)*3)];
    leafs_norm (j,:) = obj.vn(v4,:);
    total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
    j = j+1;
%     

    
    
    end    
    i = i +1;
    
    if (i <= length(lsys_branch(:,1)))
        if (previous_height < lsys_branch(i,4))
            tmp_i = tmp_i + 1;
            previous_height = lsys_branch(i,4);
            
        else
            tmp_i = 1.5;
            previous_height = lsys_branch(i,4);
        end
    end
    k = k + 1;
   
    
    tmpz = obj.v(:,2) -24;
end

end
end