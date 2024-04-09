 
clc; clear; close;

obj = Mesh();

obj.nen = 4;
obj.ndm = 2;

s2 = sqrt(2.);

% Nodes
obj.x(1,:)  = [-1.0    1.0];
obj.x(2,:)  = [ 0.0    1.0];
obj.x(3,:)  = [ 0.7    0.7];
obj.x(4,:)  = [ 1.0    0.0];
obj.x(5,:)  = [ 1.0   -1.0];
obj.x(6,:)  = [ 2.0   -1.0];
obj.x(7,:)  = [ 2.0    0.0];
obj.x(8,:)  = [ s2     s2 ];
obj.x(9,:)  = [ 0.0    2.0];
obj.x(10,:) = [-1.0    2.0];

% Elements
obj.ix(1,:) = [1     2     9    10];
obj.ix(2,:) = [2     3     8     9];
obj.ix(3,:) = [3     4     7     8];
obj.ix(4,:) = [4     5     6     7];

obj.numel = length(obj.ix);
obj.numnp = length(obj.x);

nodes = [9, 8, 7, 2, 3, 4];
obj.setCoordinateSystem(nodes, 'POLAR', [ 0,  0]);

obj.split( 1,   2,  25);
obj.split( 2,   3,  25);
obj.split( 3,   4,  25);
obj.split( 4,   5,  25);
obj.split( 1,  10,  20);

obj.smoothInternal(20);

obj.plotMesh();
