clear all;
close all;

%[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/obj_test/Fir_OBJ/Fir.obj');
%[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/tree2/Tree.obj');
%[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/47-mapletree/MapleTree.obj');
[V,F] = read_vertices_and_faces_from_obj_file('/home/hassan/Documents/47-mapletree/MapleTreeStem.obj');

%trisurf(F,V(:,1),V(:,2),V(:,3),'FaceColor',[0.26,0.33,1.0 ]);
 %trisurf(F,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ])
 
% %plot(V(:,2),V(:,3))
% 
% light('Position',[-1.0,-1.0,100.0],'Style','infinite');
% lighting phong;
 %hold on
a= rand(1);

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
 
 
V=V*roty(10);

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
%trisurf(branch_F,V(:,1),V(:,3),V(:,2)+10*a,'FaceColor',[0.26,1,0.33 ]);
%trisurf(branch_F,-V(:,1),V(:,3),V(:,2)+20*a,'FaceColor',[0.26,1,0.33 ]);


% 
for i=3073:3140
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
Sbranch_F = F(3073:3140,:);
trisurf(Sbranch_F,V(:,1),V(:,3),V(:,2),'FaceColor',[0.26,1,0.33 ]);
%trisurf(Sbranch_F,V(:,1)*0.2,V(:,3)*0.2,V(:,2)+10*a,'FaceColor',[0.26,1,0.33 ]);
%trisurf(Sbranch_F,-V(:,1)*0.2,V(:,3)*0.2,V(:,2)+20*a,'FaceColor',[0.26,1,0.33 ]);