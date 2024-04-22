
clc; clear; close;

obj = Mesh();

obj.numnp = 20;
obj.numel = 12;

obj.nen = 4;
obj.ndm = 2;

 R = 0.5;
 s = R / sqrt(2.);

% Nodes
obj.x(1,:)  = [-5, -5];
obj.x(2,:)  = [-3, -5];
obj.x(3,:)  = [ 3, -5];
obj.x(4,:)  = [ 5, -5];
obj.x(5,:)  = [-5, -3];
obj.x(6,:)  = [-3, -3];
obj.x(7,:)  = [ 3, -3];
obj.x(8,:)  = [ 5, -3];
obj.x(9,:)  = [-5,  3];
obj.x(10,:) = [-3,  3];
obj.x(11,:) = [ 3,  3];
obj.x(12,:) = [ 5,  3];
obj.x(13,:) = [-5,  5];
obj.x(14,:) = [-3,  5];
obj.x(15,:) = [ 3,  5];
obj.x(16,:) = [ 5,  5];
obj.x(17,:) = [-s, -s];
obj.x(18,:) = [ s, -s];
obj.x(19,:) = [-s,  s];
obj.x(20,:) = [ s,  s];

% Elements
obj.ix(1,:)  = [ 1,  2,  6,  5];
obj.ix(2,:)  = [ 2,  3,  7,  6];
obj.ix(3,:)  = [ 3,  4,  8,  7];
obj.ix(4,:)  = [ 5,  6, 10,  9];
obj.ix(5,:)  = [ 6, 17, 19, 10];
obj.ix(6,:)  = [ 6,  7, 18, 17];
obj.ix(7,:)  = [ 7, 11, 20, 18];
obj.ix(8,:)  = [ 7,  8, 12, 11];
obj.ix(9,:)  = [19, 20, 11, 10]; 
obj.ix(10,:) = [ 9, 10, 14, 13];
obj.ix(11,:) = [10, 11, 15, 14];
obj.ix(12,:) = [11, 12, 16, 15];

nodes = [17, 18, 19, 20];
obj.setCoordinateSystem(nodes, 'POLAR', [ 0,  0]);
obj.split(  1,  2, 4);
obj.split(  2,  3, 8);
obj.split( 14, 15, 8);
obj.split(  3,  4, 4);

obj.split(  1,  5, 4);
obj.split(  5,  9, 8);
obj.split(  8, 12, 8);
obj.split(  9, 13, 4);

obj.split(  6, 17, 8);

obj.smoothInternal(50);
obj.plotMesh();
