
clc; clear; close;

obj = Mesh();

obj.nen = 4;
obj.ndm = 2;
 
% Nodes
obj.x(1,:)  = [0, 0];
obj.x(2,:)  = [1, 0];
obj.x(3,:)  = [2, 0];
obj.x(4,:)  = [2, 1];
obj.x(5,:)  = [1, 1];
obj.x(6,:)  = [0, 1];
obj.x(7,:)  = [0, 2];
obj.x(8,:)  = [1, 2];

% Elements
obj.ix(1,:)  = [ 1,  2,  5,  6];
obj.ix(2,:)  = [ 2,  3,  4,  5];
obj.ix(3,:)  = [ 6,  5,  8,  7];

obj.numnp = length(obj.x);
obj.numel = length(obj.ix);

obj.split(  1,  2, 4);
obj.split(  2,  3, 4);

obj.split(  1,  6, 4);
obj.split(  6,  7, 4);

obj.smoothInternal(50);
obj.plotMesh();
