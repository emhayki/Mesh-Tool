
classdef Mesh < handle

    properties
        nen             %  Number of element nodes
        ndm             %  Number of dimensions
        numnp           %  Number of nodal points
        numel           %  Number of elements
        ix              %  Element connectivity matrix
        x               %  Node coordinates matrix
        coorSystem      %  Coordinate system of the mesh
    end

    methods
% ================
        function split(obj, n1, n2, s)

            if s < 2
                error('Subdivision factor must be >= 2.');
            end

           Numnp = obj.numnp;
            n2el = generateNodesToElements(obj);
            n2se = cell(Numnp, 1);
            
         spltEdg = [];
       queueElem = {};
            
            i = 1;
            while true
                if ismember(n2, n2el{n1}{i})
                    break;
                else
                    i = i + 1;
                end
                if i > length(n2el{n1})
                    error('split: invalid (n1,n2)');
                end
            end
            
            queueElem{end + 1} = n2el{n1}{i};
            queueNode = [n1 n2];
            
            while ~isempty(queueElem)
                [queueElem, queueNode, spltEdg, n2se, n2el, obj] = ...
                 ElementSplit(s, queueElem(1), queueNode, spltEdg, n2se, n2el, obj);
            end
            
            i = 1;
            j = 1;
            p = [];
            ln2el = size(n2el,1);

            while ~isempty(n2el)
                currentElement = n2el{i}{1};
                         found = false;

                for elemIndex = 1:size(p,1)
                    if p(elemIndex,:) == currentElement
                        found = true;
                        break;
                    end
                end

                if ~found
                    p(j,:) = currentElement;
                    j = j + 1;
                end

                n2el{i}{1} = [];
                i = i + 1;

                if i > ln2el
                    for kk = 1:i-1
                        n2el{kk} = n2el{kk}(~cellfun('isempty', n2el{kk}));
                    end
                 n2el = n2el(~cellfun('isempty', n2el));
                ln2el = size(n2el,1);
                    i = 1;
                end
            end
            
            obj.ix = unique(p, 'rows', 'first');
            X = zeros(length(n2el), obj.ndm);
            X(1:size(obj.x,1), :) = obj.x;
            obj.numel = size(obj.ix,1);
            
            if isempty(obj.coorSystem)
                if obj.ndm == 2
                    obj = setCoordinateSystem1(obj, 'CARTESIAN2D', X);
                elseif obj.ndm == 3
                    obj = setCoordinateSystem1(obj, 'CARTESIAN3D', X);
                else
                    error('split: invalid ndm!');
                end
            else
                for i = 1:Numnp
                    obj.coorSystem(i) = getFromCartesian(obj, obj.coorSystem(i), X(i, :));
                end
            end

            
            dfact = 1. / s;
            Numnp = Numnp + 1;
            for i = 1:3:length(spltEdg)
                i1 = min(spltEdg(i), spltEdg(i+1));
                i2 = max(spltEdg(i), spltEdg(i+1));
                fact1 = 1.;
                fact2 = 0.;
                for k = 1:s-1
                    fact1 = fact1 - dfact;
                    fact2 = fact2 + dfact;
                    obj.coorSystem(Numnp) = obj.interpolate(obj.coorSystem(i1), obj.coorSystem(i2), fact1, fact2);
                    X(Numnp,:) = obj.transformToCartesian(obj.coorSystem(Numnp));
                    Numnp = Numnp + 1;
                end
            end
            
            obj.x = X;
            obj.numnp = size(obj.x,1);
            obj.ix = findElementOrder(obj, p, obj.x);
            obj.numel = size(obj.ix,1);
            
            fprintf('The element split has been successful.\n');
            fprintf('   new numel = %d\n', obj.numel);
            fprintf('   new numnp = %d\n\n', obj.numnp);
        end
% ================
        function n2e = generateNodesToElements(obj)
            n2e = cell(obj.numnp, 1);
            for i = 1:size(obj.ix,1)
                 elem = obj.ix(i,:);
                for j = 1:obj.nen
                    nodeIndex = elem(j);
                    if isempty(n2e{nodeIndex})
                        n2e{nodeIndex} = {elem};
                    else
                        n2e{nodeIndex}{end+1} = elem;
                    end
                end 
            end 
        end
% ================
        function [queueElem, queueNode, spltEdg, n2se, n2el, obj] = ElementSplit(s, queueElem, queueNode, spltEdg, n2se, n2el, obj)
       
             n1 = queueNode(1);
             n2 = queueNode(2);
              m = [n1, n2];
            ixv = [queueElem{1} queueElem{1}(1)];

            for i = 1:5
                if isempty(setxor(m, ixv(1:2)))
                    break; 
                end
                ixv = circshift(ixv, -1);
                ixv(end) = ixv(1);
            end
                 
            n1 = ixv(1);
            n2 = ixv(2);
            n3 = ixv(3); 
            n4 = ixv(4); 
           n12 = zeros(s + 1,1);
           n43 = zeros(s + 1,1);

            n12(1) = n1;
            n43(1) = n4;
    
            n12(s + 1) = n2;
            n43(s + 1) = n3;
        
        [k, n2se, n2el, spltEdg, queueElem, queueNode] = splitEdge(obj, n1, n2, s, n2se, n2el, spltEdg, queueElem, queueNode);
            for j = 1:s-1
                n12(j + 1) = abs(k + j - 1);
            end
        
        [k, n2se, n2el, spltEdg, queueElem, queueNode] = splitEdge(obj, n4, n3, s, n2se, n2el, spltEdg, queueElem, queueNode);
        
        
            for j = 1:s-1
                n43(j + 1) = abs(k + j - 1);
            end
        
            Ix(1) = n12(1);
            Ix(2) = n12(2);
            Ix(3) = n43(2);
            Ix(4) = n43(1);
            
            bix = queueElem{1};
        
            for i = 1:length(n2el)
                for j = 1:length(n2el{i})
                    if(n2el{i}{j} == bix)
                        n2el{i}{j} = Ix;
                    end
                end
            end
        
        
        n2el = updateNodesToElements1(obj, n2el,Ix);
        
             for i = 2:s
                Ix(1) = n12(i);
                Ix(2) = n12(i + 1);
                Ix(3) = n43(i + 1);
                Ix(4) = n43(i);
                 n2el = updateNodesToElements1(obj, n2el, Ix);
             end
      
            n2el = updateNodesToElements2(obj, n2el, n3);
            n2el = updateNodesToElements2(obj, n2el, n4);
           
            queueElem(1) = [];
            queueNode(1) = [];
            queueNode(1) = [];
        
        end
% ================
        function [k, n2se, n2el, spltEdg, queueElem, queueNode] = splitEdge(~, n1, n2, s, n2se, n2el, spltEdg, queueElem, queueNode)

            nmn = min(n1, n2);
            nmx = max(n1, n2);
        
            for i = 1:length(n2se{nmn})
                if(spltEdg(n2se{nmn}(i) + 2) == nmx)
                    k = spltEdg(n2se{nmn}(i) + 3);
                    if(n1 >= n2); k = 2 - k - s; end
                    return;
                end
            end
        
            k = length(n2se) + 1;
            n2se{k + s - 2, 1} = [];
            n2se{nmn} = [n2se{nmn} length(spltEdg)];
            spltEdg = [spltEdg nmn nmx k];
        
            for i = 1:length(n2el{n1})
                elemIndex = n2el{n1}{i}; 
                 contains = any(cellfun(@(x) isequal(x, elemIndex), queueElem));
                if ~contains
                    if ismember(n2, elemIndex) 
                        queueElem{end + 1} = elemIndex;
                        queueNode = [queueNode, n1, n2]; 
                    end
                end
            end

            if n1 >= n2
                k = 2 - k - s;
            end

        end
% ================
        function n2e = updateNodesToElements1(~, n2e, ix)
            for i = 1:length(ix)
                k = ix(i);
                if (k > length(n2e))
                    for j = (length(n2e) + 1):k
                    n2e{j,1} = {};
                    end
                end
                n2e{k} = [n2e{k} ix];
            end
        end
% ================
        function n2e = updateNodesToElements2(~, n2e, n)
            i = 1;
            while i <= length(n2e{n})     
                if ~any(n2e{n}{i} == n) 
                    n2e{n}(i) = []; 
                else
                    i = i + 1;
                end
            end
        end
% ================
        function obj = setCoordinateSystem1(obj, t, x)
                     Numnp = obj.numnp;
            obj.coorSystem = repmat(struct('type', '', 'xc', [], 'coor', []), Numnp, 1); 
            for i = 1:Numnp
                obj.coorSystem(i).type = t;
                obj.coorSystem(i) = getFromCartesian(obj, obj.coorSystem(i), x(i, :));
            end
        end
% ================
        function elem = getFromCartesian(~, elem, x)
            switch elem.type
                case 'CARTESIAN3D'
                    elem.coor = x(1:3); 
        
                case 'CARTESIAN2D'
                    elem.coor = x(1:2); 
        
                case 'POLAR'
                    dx = x(1) - elem.xc(1);
                    dy = x(2) - elem.xc(2);
                    elem.coor = [atan2(dy, dx), sqrt(dx^2 + dy^2)];
        
                case 'SPHERICAL'
                    dx = x(1) - elem.xc(1);
                    dy = x(2) - elem.xc(2);
                    dz = x(3) - elem.xc(3);
                    elem.coor = [sqrt(dx^2 + dy^2 + dz^2), atan2(dy, dx), atan2(sqrt(dx^2 + dy^2), dz)];
            end
        end
% ================
        function cs = interpolate(obj, c1, c2, f1, f2)

            if strcmp(c1.type, c2.type)
                
                coor1 = [c1.coor(1), c1.coor(2)];
                coor2 = [c2.coor(1), c2.coor(2)];
            
                type = c1.type;
            
                if strcmp(type, 'SPHERICAL')
                    error('SPHERICAL interpolation is not implemented.');
                end
            
                if strcmp(type, 'POLAR') || strcmp(type, 'CYLINDRICAL')
                   
                    if coor2(1) + pi + 1e-4 < coor1(1)
                       coor2(1) = coor2(1) + 2*pi;
                        
                    elseif coor1(1) + pi + 1e-4 < coor2(1)
                           coor1(1) = coor1(1) + 2*pi;
            
                    end
                end
            
                coor = f1 * coor1 + f2 * coor2;
                xc = c1.xc;
            
            else
            
                x1 = obj.transformToCartesian(c1);
                x2 = obj.transformToCartesian(c2);
                
                coor = zeros(1, 2);
                for i = 1:2
                    coor(i) = f1 * x1(i) + f2 * x2(i);
                end
                
                if strcmp(c1.type, 'POLAR') || strcmp(c1.type, 'CARTESIAN2D') || strcmp(c2.type, 'POLAR') || strcmp(c2.type, 'CARTESIAN2D')
                    type = 'CARTESIAN2D';
                else
                    type = 'CARTESIAN3D';
                end
                
            end
            
            cs.type = type;
            cs.coor = coor;
            cs.xc = [ 0,  0]; 
        end
% ================
        function x = transformToCartesian(~, csElement)

            switch csElement.type
                case 'CARTESIAN3D'
                    x(1) = csElement.coor(1);
                    x(2) = csElement.coor(2);
                    x(3) = csElement.coor(3);
        
                case 'CARTESIAN2D'
                    x(1) = csElement.coor(1);
                    x(2) = csElement.coor(2);
        
                case 'CYLINDRICAL'
                    x(1) = csElement.xc(1) + cos(csElement.coor(1)) * csElement.coor(2);
                    x(2) = csElement.xc(2) + sin(csElement.coor(1)) * csElement.coor(2);
                    x(3) = csElement.coor(3);
        
                case 'POLAR'
                    x(1) = csElement.xc(1) + cos(csElement.coor(1)) * csElement.coor(2);
                    x(2) = csElement.xc(2) + sin(csElement.coor(1)) * csElement.coor(2);
  
                case 'SPHERICAL'
                    x(1) = csElement.xc(1) + cos(csElement.coor(2)) * sin(csElement.coor(3)) * csElement.coor(1);
                    x(2) = csElement.xc(2) + sin(csElement.coor(2)) * sin(csElement.coor(3)) * csElement.coor(1);
                    x(3) = csElement.xc(3) + cos(csElement.coor(3)) * csElement.coor(1);
            end
        end
% ================
        function orderedElements = findElementOrder(~, elements, nodeCoords)

              numElements = size(elements, 1);
                centroids = zeros(numElements, 2);
                
                for i = 1:numElements
                    elementNodes = elements(i, :);
                    xMean = mean(nodeCoords(elementNodes, 1));
                    yMean = mean(nodeCoords(elementNodes, 2));
                    centroids(i, :) = [xMean, yMean];
                end
                
                     [~, sortOrder] = sortrows(centroids, [2, 1]);
                elementsSortedByRow = elements(sortOrder, :);
               centroidsSortedByRow = centroids(sortOrder, :);
                
                uniqueY = unique(centroidsSortedByRow(:, 2));
                orderedElements = [];
                
                for i = 1:length(uniqueY)
                       rowIndices = find(centroidsSortedByRow(:, 2) == uniqueY(i));
                    elementsInRow = elementsSortedByRow(rowIndices, :);
                    [~, rowSortOrder] = sort(centroidsSortedByRow(rowIndices, 1));
                    if mod(i, 2) == 0
                        rowSortOrder = rowSortOrder(end:-1:1);
                    end
                    orderedElements = [orderedElements; elementsInRow(rowSortOrder, :)];
                end
        end
% ================
        function obj = smoothInternal(obj, iter)
               Numnp = obj.numnp;
                   X = obj.x;
                 n2n = generateNodesToNodes(obj);
            
                if (obj.ndm == 2)
                    bndNodeToNode = ComputeBoundaryData2D(obj);
                elseif(obj.ndm == 3)
                   error('smoothInternal: 3D not implemented');
                else
                   error('smoothInternal: 2D and 3D only!');
                end
                    
                for it = 1:iter
                    for i = 1:Numnp
                        if isempty(bndNodeToNode{i}(:))
                            xi = X(n2n{i}(1), :);
                            for j = 2:length(n2n{i})
                               xj = X(n2n{i}(j), :);
                               xi = xi + xj;
                            end
                             fact = 1. / length(n2n{i});
                               xi = xi * fact;
                          X(i, :) = xi;
                        end
                    end
                end 
                obj.x = X;
        end
% ================
        function n2n = generateNodesToNodes(obj)
                  Ix = obj.ix;
            numNodes = max(Ix(:));
            n2n = cell(numNodes, 1);
            for elem = 1:size(Ix, 1)
                for n = 1:size(Ix, 2)
                    currentNode = Ix(elem, n);
                    n2n{currentNode} = unique([n2n{currentNode}, Ix(elem, Ix(elem, :) ~= currentNode)]);
                end
            end
        end
% ================
        function bndNodeToNode = ComputeBoundaryData2D(obj)
                Numnp = obj.numnp;
                Numel = obj.numel;
                  tmp = cell(1, Numnp);
            for elemIndex = 1:Numel
                edge = getMeshEdges(obj, obj.ix(elemIndex, :), obj.nen);

               for i = 1:2:length(edge)
                   n1 = edge(i);
                   n2 = edge(i + 1);
                    if ~isempty(tmp{n2}) && any(tmp{n2} == n1)
                        tmp{n2}(tmp{n2} == n1) = [];
                    else
                        tmp{n1}(end+1) = n2;
                    end
                end
            end

            bndNodeToNode = cell(size(tmp));
            for i = 1:Numnp
                for j = 1:length(tmp{i})
                     bndNodeToNode{i} = [ bndNodeToNode{i}, tmp{i}(j)];
                     bndNodeToNode{tmp{i}(j)} = [ bndNodeToNode{tmp{i}(j)}, i];
                end
            end
        end
% ================
        function edge = getMeshEdges(~, ix, nen)
            switch nen
                case 3
                    edge = [ix(1), ix(2), ix(2), ix(3), ix(3), ix(1)];
                case 4  
                    edge = [ix(1), ix(2), ix(2), ix(3), ix(3), ix(4), ix(4), ix(1)];
                otherwise
                    error('getMeshEdges: Element Type not implemented');
            end
        end
% ================
        function obj = setCoordinateSystem(obj, node, t, xx)
            Numnp = obj.numnp; 

                if (isempty(obj.coorSystem) || length(obj.coorSystem) < Numnp)
                    obj.coorSystem = repmat(struct('type', '', 'xc', xx, 'coor', []), Numnp, 1); 

                    if obj.ndm == 2
                        for i = 1:Numnp
                            obj.coorSystem(i).type = 'CARTESIAN2D';
                        end

                    elseif obj.ndm == 3
                        for i = 1:Numnp
                            obj.coorSystem(i).type = 'CARTESIAN3D';
                        end
                    end
            
                    for i = 1:Numnp
                        obj.coorSystem(i) = getFromCartesian(obj, obj.coorSystem(i),obj.x(i,:)); 
                    end
                end
            
            
                for i = 1:length(node)
                    if ~isempty(xx)
                        for j = 1:obj.ndm
                            obj.coorSystem(node(i)).xc(j) = xx(j);
                        end
                    end
                    obj.coorSystem(node(i)) = getFromCartesian(obj, obj.coorSystem(i),obj.x(node(i),:));
                    obj.coorSystem(node(i)).type = t;
                end
        end
% ================
        function plotMesh(obj, varargin)

            showNodes    = false;
            showElements = false;
        
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'nodes')
                    showNodes = true;
                elseif strcmp(varargin{i}, 'elements')
                    showElements = true;
                end
            end
        
            patch('Vertices',  obj.x, 'Faces',  obj.ix, 'FaceColor', 'white', 'EdgeColor', 'black');
            hold on;
        
            if showNodes
                 plot(obj.x(:,1), obj.x(:,2), 'ro');
                for i = 1:size(obj.x, 1)
                    text(obj.x(i,1), obj.x(i,2), sprintf('%d', i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color', 'red', 'Clipping', 'on');
                end
            end
       
            if showElements
                for i = 1:size(obj.ix, 1)
                    centroid = mean(obj.x(obj.ix(i,:), :), 1);
                    text(centroid(1), centroid(2), sprintf('%d', i), ...
                         'VerticalAlignment', 'cap', 'HorizontalAlignment', 'center', 'Color', 'blue', 'Clipping', 'on');
                end
            end
        
               hold off;
               axis equal;

            xRange = range(obj.x(:,1));
            yRange = range(obj.x(:,2));
            maxRange = max(xRange, yRange);
            margin = maxRange * 0.05; 
            xlim([min(obj.x(:,1)) - margin, max(obj.x(:,1)) + margin]);
            ylim([min(obj.x(:,2)) - margin, max(obj.x(:,2)) + margin]);
        end
% ================
    end

end
