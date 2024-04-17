MeshTool is a MATLAB class library designed for advanced mesh manipulation and analysis. It enables users to dynamically modify meshes, such as splitting elements and interpolating node coordinates, making it a valuable tool for computational geometry and finite element analysis.

## Features

- **Dynamic Mesh Modification:** Supports operations like element splitting, allowing for localised mesh refinement.
- **Flexible Coordinate Systems:** Supports multiple coordinate systems (Cartesian, Polar, Spherical).
- **Visualisation:** Includes built-in methods for mesh visualisation, with options to display nodes and element numbering.

## Installation

To use MeshTool, clone this repository or download the `Mesh.m` file into your MATLAB working directory:

```bash
git clone https://github.com/emhayki/MeshTool.git
```

## Usage
To use the tool, follow these steps:

1. **Create a Mesh Instance:** Instantiate the **Mesh** object. This will be the primary object for your mesh manipulations.
```
obj = Mesh();
```


2. **Define Mesh Properties:** Set the number of element nodes (**nen**) and dimensions (**ndm**) for the mesh. In this example, a 2D mesh composed of quadrilateral elements is being created.
```
obj.nen = 4;
obj.ndm = 2;
```

3. **Specify Node Coordinates:** Define the coordinates of each node in the mesh. 
```
obj.x(1,:)  = [0, 0];
obj.x(2,:)  = [1, 0];
obj.x(3,:)  = [2, 0];
obj.x(4,:)  = [2, 1];
obj.x(5,:)  = [1, 1];
obj.x(6,:)  = [0, 1];
obj.x(7,:)  = [0, 2];
obj.x(8,:)  = [1, 2];
```

4. **Define Elements by Nodes:** Construct elements by specifying the indices of the nodes that form each element. 
```
obj.ix(1,:)  = [ 1,  2,  5,  6];
obj.ix(2,:)  = [ 2,  3,  4,  5];
obj.ix(3,:)  = [ 6,  5,  8,  7];
```

5. **Update Node and Element Counts:** Update the total number of nodes and elements.
```
obj.numnp = length(obj.x);
obj.numel = length(obj.ix);
```

6. **Refine the Mesh:** Use the **split** method to refine specific elements by dividing them into smaller elements. 
```
obj.split(1, 2, 4); % Split between nodes 1 and 2 into 4 segments
obj.split(2, 3, 4); % Split between nodes 2 and 3 into 4 segments
obj.split(1, 6, 4); % Split between nodes 1 and 6 into 4 segments
obj.split(6, 7, 4); % Split between nodes 6 and 7 into 4 segments
```


7. **Smooth the Mesh:** Optionally, smooth the mesh to improve the quality of the mesh after splitting. This is achieved by averaging node positions over a specified number of iterations.
```
obj.smoothInternal(50);
```

8. **Visualise the Mesh:** Display the refined mesh using the plotMesh method. To view node numbers, use `obj.plotMesh('nodes');`, to see element numbers, use `obj.plotMesh('elements');`, or to display both, use `obj.plotMesh('nodes', 'elements');`.

![Mesh_Refinement](https://github.com/emhayki/MeshTool/assets/135982304/6c2481ee-0a90-4c55-99a1-0bcc4037dc24)

The runnable code for this example is provided in `example1.m`. 

## Examples
The following table showcases various meshes that can be generated using MeshTool.

|[Example 3](/example3.m) | [Example 4](/example4.m)  | [Example 5](/example5.m)  |
| ------------ | -------------------- | -------------------- |
|![example3](https://github.com/emhayki/MeshTool/assets/135982304/21bbb120-f848-4509-98f2-77b318a0b5a7)|![example4](https://github.com/emhayki/MeshTool/assets/135982304/f86fc5a7-e3e6-45d3-980f-cd9cb2647c4d)|![example5](https://github.com/emhayki/MeshTool/assets/135982304/425d64dc-ff37-4378-9e1c-676108530e97)|



## Contributing
Contributions are welcome! Please feel free to submit pull requests with bug fixes, improvements, or new features.

