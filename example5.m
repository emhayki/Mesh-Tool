 
clc; clear; close;

obj = Mesh();

obj.nen = 4;
obj.ndm = 2;

s2 = sqrt(2.) + 0.025;
s3 = 0.6 + s2;

% Nodes
obj.x(1,:)  = [ -.5    1.3 ];
obj.x(2,:)  = [  0.5   1.3 ];
obj.x(3,:)  = [  1.3   0.64];
obj.x(4,:)  = [  1.5   0.0 ];
obj.x(5,:)  = [  2.5   0.0 ];
obj.x(6,:)  = [  s3     s2 ];
obj.x(7,:)  = [  0.5   2.4];
obj.x(8,:)  = [ -0.5   2.4];
obj.x(9,:)  = [ -1.3   0.64];
obj.x(10,:) = [ -1.5   0.0 ];
obj.x(11,:) = [ -2.5   0.0 ];
obj.x(12,:) = [ -s3     s2 ];

% Elements
obj.ix(1,:) = [1     2    7   8];
obj.ix(2,:) = [2     3    6   7];
obj.ix(3,:) = [3     4    5   6];
obj.ix(4,:) = [12    9    1   8];
obj.ix(5,:) = [11   10    9  12];

obj.numel = length(obj.ix);
obj.numnp = length(obj.x);

nodes = [1  2  3  4  5  6  7  8  9    10    11    12];
obj.setCoordinateSystem(nodes, 'POLAR', [ 0,  0]);

obj.split( 1,   2,  14);
obj.split( 9,   1,  20);
obj.split( 2,   3,  20);
obj.split( 9,  10,  20);
obj.split( 3,   4,  20);
obj.split( 11, 10,  15);

obj.smoothInternal(20);
obj.plotMesh();
