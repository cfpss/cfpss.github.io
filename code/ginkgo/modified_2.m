function [total_leaf,xn,yn,radius] = modified_2(newBR)
%obj = readObj('/home/hassan/Documents/obj_test/Fir_OBJ/Fir.obj');   %use appropriate file
%obj = readObj('/home/hassan/Documents/tree2/Tree.obj');
obj = readObj('/home/hassan/Documents/47-mapletree/MapleTree.obj');
%obj = readObj('/home/hassan/Documents/47-mapletree/MapleTreeStem.obj');

obj.v(:,1) = obj.v(:,1);
obj.v(:,2) = obj.v(:,2);
obj.v(:,3) = obj.v(:,3);

shading interp
colormap(gray(256));
lighting phong
camproj('perspective');
axis square
axis off
axis equal
axis tight
cameramenu



 face = [];
 total_leaf = [];

figure

% for FIR Tree
obj.f.v = obj.f.v'
obj.f.vn = obj.f.vn'
%for i=1:length(obj.f.v)

%for i=1:length(obj.f.v)-10000

rand_leaf = [];
bnum = 8;
y = 1;

for lsys = 1:length(newBR)
if (newBR(lsys,2) == 0 && newBR(lsys,3) == 0 && newBR(lsys,5) == 0 && newBR(lsys,6) == 0)
    continue;
end

lsys_branch(y,:) = newBR(lsys,:);
y = y+1;
end

 for a = 1:length(lsys_branch)
%     if (lsys_branch(a,5) == 0)
% newangle(a) = -abs(atan2(lsys_branch(a,6),lsys_branch(a,5)));
%     else
         newangle(a) = atan2(lsys_branch(a,6),lsys_branch(a,5));
%     end
 end
 %for i=1:300  
 IPP = IPP_points();
 %IPP = IPP(1:2,:);
 
xn = 0;
yn = 0;

for i = 1:length(IPP(:,1))
xn = xn + IPP(i,1);
yn = yn + IPP(i,2);
end

xn = xn / length(IPP(:,1));
yn = yn / length(IPP(:,2));

d = 0;
for i = 1:length(IPP(:,1))
d1 = sqrt((xn-IPP(i,1))^2 +(yn-IPP(i,2))^2);

if d1 > d
    x1 = IPP(i,1);
    y1 = IPP(i,2);
d = d1;
end
end
radius = 10 + sqrt((x1 - xn)^2 + (y1 - yn)^2);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%% trunk %%%%%%%%%%%%%%%%%%%%%%%%%%%
 for i=1:192
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
trunk = face(1:192,:);

% IPP(2,1)
for ipp = 1:length(IPP(:,1))
tmp_1 = obj.v(:,1) + IPP(ipp,1);
tmp_2 = obj.v(:,3) + IPP(ipp,2);

    trisurf(trunk,tmp_1,tmp_2,obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);

end
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
%branch_F = face(193:312,:);

% trisurf(face,obj.v(:,1),obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);
% trisurf(face,-obj.v(:,1),-obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);

tmpz = obj.v(:,2);
% for i = 1:bnum
   
tmp_angle = 0;

i = 1;

previous_height = lsys_branch(1,4);
tmp_i = 1;
x1 = obj.v(:,1);
x2 = obj.v(:,3);

rotation = 0;

    while i <= length(lsys_branch)
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
   


if (i == 1)
tmp_angle = newangle(i);
end

rotation = rotation + newangle(i) - tmp_angle;

for ipp = 1:length(IPP(:,1))
tmp_1 = IPP(ipp,1);
tmp_2 = IPP(ipp,2);

trisurf(face,...
(cos(rotation)*obj.v(:,1)+sin(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_1,...
(-sin(rotation)*obj.v(:,1)+cos(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_2,...
tmpz + lsys_branch(i,4)*3,...
'FaceColor',[0.26,0.33,1.0 ]);

end


tmp_angle = newangle(i);
i = i+1;

rotation = rotation + newangle(i) - tmp_angle;

for ipp = 1:length(IPP(:,1))
tmp_1 = IPP(ipp,1);
tmp_2 = IPP(ipp,2);

trisurf(face,...
(cos(rotation)*obj.v(:,1)+sin(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_1,...
(-sin(rotation)*obj.v(:,1)+cos(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_2,...
tmpz + lsys_branch(i,4)*3,...
'FaceColor',[0.26,0.33,1.0 ]);

end

tmp_angle = newangle(i);

i = i+1;

if (i <= length(lsys_branch))
if (previous_height < lsys_branch(i,4))
tmp_i = tmp_i + 1;
previous_height = lsys_branch(i,4);

else
tmp_i = 1.5; 
previous_height = lsys_branch(i,4);
end
end

tmpz = obj.v(:,2);


   
%     trisurf(face,...    
% (cos(2*pi*i/bnum)*obj.v(:,1)+sin(2*pi*i/bnum)*obj.v(:,3))*(10-1.2*i)/10,...
% (-sin(2*pi*i/bnum)*obj.v(:,1)+cos(2*pi*i/bnum)*obj.v(:,3))*(10-1.2*i)/10,...
% tmpz+3*(i-1),...
% 'FaceColor',[0.26,0.33,1.0 ]);
%     trisurf(face,...    
% (cos(2*pi*i/bnum)*-obj.v(:,1)+sin(2*pi*i/bnum)*-obj.v(:,3))*(10-1.2*i)/10,...
% (-sin(2*pi*i/bnum)*-obj.v(:,1)+cos(2*pi*i/bnum)*-obj.v(:,3))*(10-1.2*i)/10,...
% tmpz+3*(i-1),...
% 'FaceColor',[0.26,0.33,1.0 ]);
% tmpz = obj.v(:,2);

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
%Sbranch_F = face(3073:3863,:);
%trisurf(face,obj.v(:,1),obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);
%trisurf(face,-obj.v(:,1),-obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);

tmpz = obj.v(:,2);
 
    
  tmp_angle = 0;

i = 1;

previous_height = lsys_branch(1,4);
tmp_i = 1;
x1 = obj.v(:,1);
x2 = obj.v(:,3);

rotation = 0;

    while i <= length(lsys_branch)
           count = 1;
    rand_SB = rand(1);    
    
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
  
if (i == 1)
tmp_angle = newangle(i);
end

rotation = rotation + newangle(i) - tmp_angle;


for ipp = 1:length(IPP(:,1))
tmp_1 = IPP(ipp,1);
tmp_2 = IPP(ipp,2);

trisurf(face,...
(cos(rotation)*obj.v(:,1)+sin(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_1,...
(-sin(rotation)*obj.v(:,1)+cos(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_2,...
tmpz + lsys_branch(i,4)*3,...
'FaceColor',[0.26,0.33,1.0 ]);

end

tmp_angle = newangle(i);
i = i+1;

rotation = rotation + newangle(i) - tmp_angle;

for ipp = 1:length(IPP(:,1))
tmp_1 = IPP(ipp,1);
tmp_2 = IPP(ipp,2);

trisurf(face,...
(cos(rotation)*obj.v(:,1)+sin(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_1,...
(-sin(rotation)*obj.v(:,1)+cos(rotation)*obj.v(:,3))*(10-2.2*tmp_i)/10 + tmp_2,...
tmpz + lsys_branch(i,4)*3,...
'FaceColor',[0.26,0.33,1.0 ]);
end


tmp_angle = newangle(i);

i = i+1;

if (i <= length(lsys_branch))
if (previous_height < lsys_branch(i,4))
tmp_i = tmp_i + 1;
previous_height = lsys_branch(i,4);

else
tmp_i = 1.5; 
previous_height = lsys_branch(i,4);
end
end
tmpz = obj.v(:,2);
    end



%%%%%%%%%%% Leaf %%%%%%%%%%%%%%%%
face = [];
j = 1;
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
% 
% leafs(j,:) = [x1,z1,y1];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
% 
% leafs(j,:) = [-x1,-z1,y1];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
% 
% 
% leafs(j,:) = [x2,z2,y2];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
% 
% leafs(j,:) = [-x2,-z2,y2];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
% 
% 
% leafs(j,:) = [x3,z3,y3];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
% 
% leafs(j,:) = [-x3,-z3,y3];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
end

%trisurf(face,obj.v(:,1),obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);
%trisurf(face,-obj.v(:,1),-obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);

% for k = 1:length(leafs)
%   leafs(j,:) = [-x1,-z1,y1];
% leafs_norm (j,:) = obj.vn(v4,:);
% total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
% j = j+1;
% end

  tmp_angle = 0;
tmpz = obj.v(:,2);
i = 1;

previous_height = lsys_branch(1,4);
tmp_i = 1;
x1 = obj.v(:,1);
x2 = obj.v(:,3);

  rotation1 = 0;
  k = 1;
  
while i <= length(lsys_branch)
    ii = i;
    i1 = i;
   
    count = 1;
    rand_SB = rand(1);
    rand_leaf = [rand_leaf; rand_SB];
   
    if (i == 1)
            tmp_angle = newangle(i);
    end
       
        rotation1 = rotation1 + newangle(i) - tmp_angle;
       
        tmp_angle = newangle(i);
        i = i+1;
        i2 = i;
       
        rotation2 = rotation1 + newangle(i) - tmp_angle;
     
        tmp_angle = newangle(i);
        i = i+1;
       
   
    for L = 20793:21255
        if (L== 20793 || L== 20857 || L== 20893 || L== 20933 || L== 20973 || L== 21069 ...
                || L== 21117 || L== 21165 || L== 21185 || L== 21209 || L== 21253)
           
            tmpz(L) = Sb_root{ii,count};
            count = count + 1;
            root = L;
            shf = tmpz(L);
        end
       
        tmpz(L) = shf + (L-root)*9*rand_leaf(k)*10E-03;
       
           
       for ipp = 1:length(IPP(:,1))
tmp_1 = IPP(ipp,1);
tmp_2 = IPP(ipp,2);
       
        leafs(j,:) = [((cos(rotation1)*obj.v(L,1)+sin(rotation1)*obj.v(L,3))*(10-2.2*tmp_i)/10) + tmp_1,...
            ((-sin(rotation1)*obj.v(L,1)+cos(rotation1)*obj.v(L,3))*(10-2.2*tmp_i)/10) + tmp_2,...
            (tmpz(L)+lsys_branch(i1,4)*3)];
        
            leafs_norm (j,:) = obj.vn(v4,:);
        total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
        j = j+1;
       
       end
       
 
       
       
    for ipp = 1:length(IPP(:,1))
tmp_1 = IPP(ipp,1);
tmp_2 = IPP(ipp,2); 
       
        leafs(j,:) = [((cos(rotation2)*obj.v(L,1)+sin(rotation2)*obj.v(L,3))*(10-2.2*tmp_i)/10) + tmp_1,...
            ((-sin(rotation2)*obj.v(L,1)+cos(rotation2)*obj.v(L,3))*(10-2.2*tmp_i)/10) + tmp_2,...
            (tmpz(L)+lsys_branch(i2,4)*3)];
        
     leafs_norm (j,:) = obj.vn(v4,:);
        total_leaf = [total_leaf;[leafs(j,:), leafs_norm(j,:)]];
        j = j+1;
    
    end
    
              
       
       
       
       
    end
            if (i <= length(lsys_branch))
            if (previous_height < lsys_branch(i2+1,4))
                tmp_i = tmp_i + 1;
                previous_height = lsys_branch(i2+1,4);
               
            else
                tmp_i = 1.5;
                previous_height = lsys_branch(i2+1,4);
            end
        end
      k = k + 1;
    rotation1 = rotation2;
      
  
%     trisurf(face,...    
% (cos(2*pi*i/bnum)*obj.v(:,1)+sin(2*pi*i/bnum)*obj.v(:,3))*(10-1.2*i)/10,...
% (-sin(2*pi*i/bnum)*obj.v(:,1)+cos(2*pi*i/bnum)*obj.v(:,3))*(10-1.2*i)/10,...
% tmpz+3*i,...
% 'FaceColor',[0.26,0.33,1.0 ]);
%     trisurf(face,...    
% (cos(2*pi*i/bnum)*-obj.v(:,1)+sin(2*pi*i/bnum)*-obj.v(:,3))*(10-1.2*i)/10,...
% (-sin(2*pi*i/bnum)*-obj.v(:,1)+cos(2*pi*i/bnum)*-obj.v(:,3))*(10-1.2*i)/10,...
% tmpz+3*i,...
% 'FaceColor',[0.26,0.33,1.0 ]);
tmpz = obj.v(:,2);
end

end