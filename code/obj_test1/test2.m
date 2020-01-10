clear all;
close all;

%[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/obj_test/Fir_OBJ/Fir.obj');
%[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/tree2/Tree.obj');
[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/47-mapletree/MapleTree.obj');
%[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/47-mapletree/MapleTreeStem.obj');

%trisurf(F,V(:,1),V(:,2),V(:,3),'FaceColor',[0.26,0.33,1.0 ]);
 %trisurf(F,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ])
 
% %plot(V(:,2),V(:,3))
% 
 %light('Position',[-1.0,-1.0,100.0],'Style','infinite');
 %lighting phong;
 %hold on
%a= rand(1);
rand_leaf = [];
bnum = 7;
 %for trunk
 for i=1:192

v1 = F(i,1);
v2 = F(i,2);
v3 = F(i,3);
x1 = (V(v1,1) + V(v2,1) + V(v3,1))/3;
y1 = (V(v1,2) + V(v2,2) + V(v3,2))/3;
z1 = (V(v1,3) + V(v2,3) + V(v3,3))/3;
leafs(i,:) = [x1,y1,z1];
i;
plot3(x1,z1,y1,' ');
hold on;
end
trunk = F(1:192,:);
trisurf(trunk,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
 
 
%V=V*roty(10);

% for sub
for i=193:312  
v1 = F(i,1);
v2 = F(i,2);
v3 = F(i,3);
x1 = (V(v1,1) + V(v2,1) + V(v3,1))/3;
y1 = (V(v1,2) + V(v2,2) + V(v3,2))/3;
z1 = (V(v1,3) + V(v2,3) + V(v3,3))/3;
leafs(i,:) = [x1,y1,z1];
i;
plot3(x1,z1,y1,' ');
hold on;
end

branch_F = F(193:312,:);
 trisurf(branch_F,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
trisurf(branch_F,-V(:,1),-V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);

tmpz = V(:,2);
for i = 1:bnum
    rand_b = rand(1); 
   count = 1;
   
   
   for j = 201:328
      tmpz(j) = tmpz(j) + (j-201)*5*rand_b*10E-03; 
      
      if (j== 221 || j== 234 || j== 244 || j== 251 || j== 260 || j== 273 ...
          || j== 281 || j== 291 || j== 305 || j== 314 || j== 322)
      
      Sb_root{i,count} = tmpz(j);
      count = count + 1;
      
      end
      
          
   end
  
    trisurf(branch_F,...    
(cos(2*pi*i/bnum)*V(:,1)+sin(2*pi*i/bnum)*V(:,3))*(10-1.2*i)/10,...
(-sin(2*pi*i/bnum)*V(:,1)+cos(2*pi*i/bnum)*V(:,3))*(10-1.2*i)/10,...
tmpz+3*i,...
'FaceColor',[0.26,1,0.33 ]);
    trisurf(branch_F,...    
(cos(2*pi*i/bnum)*-V(:,1)+sin(2*pi*i/bnum)*-V(:,3))*(10-1.2*i)/10,...
(-sin(2*pi*i/bnum)*-V(:,1)+cos(2*pi*i/bnum)*-V(:,3))*(10-1.2*i)/10,...
tmpz+3*i,...
'FaceColor',[0.26,1,0.33 ]);
tmpz = V(:,2);
end

% instead of 10-1.2*i we can 18-1.2*bnum-1.2*i

%hold on
%trisurf(branch_F,V(:,1)*-0.7,V(:,3)*0.7,V(:,2)+10*-a,'FaceColor',[0.26,1,0.33 ]);
%trisurf(branch_F,-V(:,1),V(:,3),V(:,2)+20*a,'FaceColor',[0.26,1,0.33 ]);


% 
for i=3073:3863
v1 = F(i,1);
v2 = F(i,2);
v3 = F(i,3);
x1 = (V(v1,1) + V(v2,1) + V(v3,1))/3;
y1 = (V(v1,2) + V(v2,2) + V(v3,2))/3;
z1 = (V(v1,3) + V(v2,3) + V(v3,3))/3;
leafs(i,:) = [x1,y1,z1];
i;
plot3(x1,z1,y1,' ');
hold on;
end

Sbranch_F = F(3073:3863,:);
 trisurf(Sbranch_F,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
 trisurf(Sbranch_F,-V(:,1),-V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
tmpz = V(:,2);
for i = 1:bnum
    count = 1;
    rand_SB = rand(1);
    rand_leaf = [rand_leaf; rand_SB]
   for j = 3273:4152
                 if (j== 3273 || j== 3353 || j== 3433 || j== 3513 || j== 3591 || j== 3674 ...
          || j== 3754 || j== 3833 || j== 3919 || j== 3990 || j== 4085)
      
       tmpz(j) = Sb_root{i,count};
      count = count + 1;
      root = j;
      shf = tmpz(j);
                 end     
%                if ((j>=3349 && j<= 3352) || (j>=3895 && j<= 3988))
%                continue;    
%                    end
               
                 
      tmpz(j) = shf + (j-root)*5*rand_SB*10E-03; 
   end
  
    trisurf(Sbranch_F,...    
(cos(2*pi*i/bnum)*V(:,1)+sin(2*pi*i/bnum)*V(:,3))*(10-1.2*i)/10,...
(-sin(2*pi*i/bnum)*V(:,1)+cos(2*pi*i/bnum)*V(:,3))*(10-1.2*i)/10,...
tmpz+3*i,...
'FaceColor',[0.26,1,0.33 ]);
 trisurf(Sbranch_F,...    
(cos(2*pi*i/bnum)*-V(:,1)+sin(2*pi*i/bnum)*-V(:,3))*(10-1.2*i)/10,...
(-sin(2*pi*i/bnum)*-V(:,1)+cos(2*pi*i/bnum)*-V(:,3))*(10-1.2*i)/10,...
tmpz+3*i,...
'FaceColor',[0.26,1,0.33 ]);

tmpz = V(:,2);
end
%  
%hold on


%trisurf(Sbranch_F,V(:,1)*-0.7,V(:,3)*0.7,V(:,2)+10*-a,'FaceColor',[0.26,1,0.33 ]);
%trisurf(Sbranch_F,-V(:,1),V(:,3),V(:,2)+20*a,'FaceColor',[0.26,1,0.33 ]);


%for leaf
for i=18841:18956
v1 = F(i,1);
v2 = F(i,2);
v3 = F(i,3);
x1 = (V(v1,1) + V(v2,1) + V(v3,1))/3;
y1 = (V(v1,2) + V(v2,2) + V(v3,2))/3;
z1 = (V(v1,3) + V(v2,3) + V(v3,3))/3;
leafs(i,:) = [x1,y1,z1];
i;
plot3(x1,z1,y1,' ');
hold on;
%drawnow
end

SBLeaf = F(18841:18956,:);
trisurf(SBLeaf,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
trisurf(SBLeaf,-V(:,1),-V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
tmpz = V(:,2);
for i = 1:bnum
    count = 1;
     for j = 20793:21255
                 if (j== 20793 || j== 20857 || j== 20893 || j== 20933 || j== 20973 || j== 21069 ...
          || j== 21117 || j== 21165 || j== 21185 || j== 21209 || j== 21253)
      
       tmpz(j) = Sb_root{i,count};
      count = count + 1;
      root = j;
      shf = tmpz(j);
                 end     
          
                 
      tmpz(j) = shf + (j-root)*9*rand_leaf(i)*10E-03; 
   end   
   
    trisurf(SBLeaf,...    
(cos(2*pi*i/bnum)*V(:,1)+sin(2*pi*i/bnum)*V(:,3))*(10-1.2*i)/10,...
(-sin(2*pi*i/bnum)*V(:,1)+cos(2*pi*i/bnum)*V(:,3))*(10-1.2*i)/10,...
tmpz+3*i,...
'FaceColor',[0.26,1,0.33 ]);
    trisurf(SBLeaf,...    
(cos(2*pi*i/bnum)*-V(:,1)+sin(2*pi*i/bnum)*-V(:,3))*(10-1.2*i)/10,...
(-sin(2*pi*i/bnum)*-V(:,1)+cos(2*pi*i/bnum)*-V(:,3))*(10-1.2*i)/10,...
tmpz+3*i,...
'FaceColor',[0.26,1,0.33 ]);
tmpz = V(:,2);
end



  


