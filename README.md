# Foil Solver

My potential solver started as an Airfoil solver. At least what I had planned for my course Numerical Methods.  It turned out to be a full potential solver for a defined field with a given object. I began by creating the field, then defined objects, and ideas just exploded in my mind with ideas afterwards that I couldn't complete the ideas. I need to admit, that it is not a full-scale potential solver yet. There are no sources or sinks, pressure solver etc. Even if my project lacks these tools, I am sure that they can be easily implemented. So, in this course project, you can create objects in a field, and you can arrange boundary conditions Neumann or Dirichlet for any potential or you can directly go for an airfoil flow tool to see foil overflow. I also build a C-type meshing tool but I haven't tried to use its Jacobian yet. 

## Getting Started

To be honest, at the moment the potential solver fully depends on square mesh, and the solver visits every cell. I was a bit late to write full-scale 2D TDMA, I was committed to doing it I hit a deadline. It is quite slow at the moment unfortunately but TDMA is always a good friend for us in this kind of situation. My other big plan was to write a suitable GUI for fancy-looking tools. But again deadline was strict.

## The Capabilities 


![Screenshot_1](https://github.com/MehmetHakanSari/PotentialSolver/assets/112701635/3f06afe8-5ec2-4183-aa60-9dd939493daa)


