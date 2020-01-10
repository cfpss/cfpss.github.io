function SecondBR=zhu_branch_rotate_plot(BR, Raa, Rab, Rb)
% This function take the Branch parameters, rotate each 90% and 
% plot the branching structure. 

RotM = rotz(90);
N = size(BR,1); % number of branches.    
xlim0=[min([BR(:,3);BR(:,5)]),max([BR(:,3);BR(:,5)])];
ylim0=[min([BR(:,4);BR(:,6)]),max([BR(:,4);BR(:,6)])];
zlim0=[min([BR(:,5);BR(:,7)]),max([BR(:,5);BR(:,7)])];
SecondBR = zeros(size(BR,1), size(BR,2));
% plot the first branch. 
endrotated = RotM*BR(1, 2:4)'; endrotated(3) = endrotated(3) + Raa*Rab*Rb; 
tiprotated = RotM*BR(1, 5:7)'; tiprotated(3) = tiprotated(3) + Raa*Rab*Rb;
cylinder2P(BR(1, 8), 10, BR(1, 2:4), BR(1, 5:7));
cylinder2P(BR(1, 8), 10, endrotated', tiprotated');   
SecondBR(1,:) = [BR(1,1) endrotated' tiprotated' BR(1,8)];
xlim(xlim0);ylim(ylim0);zlim(zlim0);
%shading flat; camlight;
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
hold on
for n = 2:N
    hold on
    endrotated = RotM*BR(n, 2:4)';
    endrotated(3) = endrotated(3) + Raa*Rab*Rb; 
    tiprotated = RotM*BR(n, 5:7)';
    tiprotated(3) = tiprotated(3) + Raa*Rab*Rb;
    cylinder2P(BR(n, 8), 10, BR(n, 2:4), BR(n, 5:7));
    cylinder2P(BR(n, 8), 10, endrotated', tiprotated');   
    cylinder2P(BR(n, 8), 10, BR(n, 2:4), BR(n, 5:7));    
    SecondBR(n,:) = [BR(n,1) endrotated' tiprotated' BR(n,8)];
end
view(-145,20);axis equal;
camlight; 
axis off
