function total_leaf = modified_1()
%obj = readObj('/home/hassan/Documents/obj_test/Fir_OBJ/Fir.obj');   %use appropriate file


%obj = readObj('/home/hassan/Documents/tree2/Tree.obj');
obj = readObj('/home/hassan/Documents/47-mapletree/MapleTree.obj');

%obj = readObj('/home/hassan/Documents/47-mapletree/MapleTreeStem.obj');



obj.v(:,1) = obj.v(:,1);
obj.v(:,2) = obj.v(:,2);
obj.v(:,3) = obj.v(:,3);

%patch('Faces', obj.f.v', 'Vertices', obj.v);
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


% for FIR Tree
obj.f.v = obj.f.v'
obj.f.vn = obj.f.vn'
%for i=1:length(obj.f.v)

%for i=1:length(obj.f.v)-10000
 for i=1:300  
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

   
%leafs(i,:) = [x1,y1,z1];
%leafs_norm (i,:) = obj.vn(v4,:);

%total_leaf = [total_leaf;[leafs(i,:), leafs_norm(i,:)]];
plot3(x1,z1,y1,' ');
hold on;
end

trisurf(face,obj.v(:,1),obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);

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

   
%leafs(i,:) = [x1,y1,z1];
%leafs_norm (i,:) = obj.vn(v4,:);

%total_leaf = [total_leaf;[leafs(i,:), leafs_norm(i,:)]];
plot3(x1,z1,y1,' ');
hold on;
end

trisurf(face,obj.v(:,1),obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);

for i=18841:18956  
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

   
leafs(i,:) = [x1,y1,z1];
leafs_norm (i,:) = obj.vn(v4,:);

total_leaf = [total_leaf;[leafs(i,:), leafs_norm(i,:)]];
%plot3(x1,z1,y1,' ');
hold on;
end

%trisurf(face,obj.v(:,1),obj.v(:,3),obj.v(:,2),'FaceColor',[0.26,0.33,1.0 ]);
end