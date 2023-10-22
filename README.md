# Foil Solver

My potential solver started as an Airfoil solver. At least what I had planned for my course Numerical Methods.  It turned out to be a full 2D potential solver for a defined field with a given object. I began by creating the field, then defined objects, and ideas just exploded in my mind with ideas afterwards that I couldn't complete the ideas. I need to admit, that it is not a full-scale potential solver yet. There are no sources or sinks, pressure solver etc. Even if my project lacks these tools, I am sure that they can be easily implemented. So, in this course project, you can create objects in a field, and you can arrange boundary conditions Neumann or Dirichlet for any potential or you can directly go for an airfoil flow tool to see foil overflow. I also build a C-type meshing tool but I haven't tried to use its Jacobian yet. 

## Getting Started

To be honest, at the moment the potential solver fully depends on square mesh, and the solver visits every cell. I was a bit late to write full-scale 2D TDMA, I was committed to doing it I hit a deadline. It is quite slow at the moment unfortunately but TDMA is always a good friend for us in this kind of situation. My other big plan was to write a suitable GUI for fancy-looking tools. But again deadline was strict.

## The Capabilities 

We start by defining the computational space. For this case, it consists of simply two spheres and one rectangle:  

![Ekran görüntüsü 2023-10-22 152224](https://github.com/MehmetHakanSari/PotentialSolver/assets/112701635/c0ed1a32-b152-40a4-84a0-f31a68c9ebb8)

We can define boundary conditions in both Dirichlet and Neumann. The North wall is the Neumann BB whereas others are given as Dirichlet. Then we solve the potential inside the domain. The BC of the objects are regarded as walls. Then a streamplot can be plotted by ready streamplot functions provided by Mathplotlib. 

![Ekran görüntüsü 2023-10-22 152235](https://github.com/MehmetHakanSari/PotentialSolver/assets/112701635/d8b2333f-dd87-4e7e-a87f-8d4c20f2d3fc)

The ultimate goal is to solve for an airfoil. The airfoil is created via the formulations provided on the web. For the stream solver, we need to solve the grid first for grid transformation. The transformed grid is given such in the figure. First, interpolation between the specified outer domain and the object of interest is employed to speed up the grid solver. 

![Ekran görüntüsü 2023-10-22 151921](https://github.com/MehmetHakanSari/PotentialSolver/assets/112701635/e626b9de-41e8-4133-9431-18d333cd47b4)

Afterwards, the stream is solved for the grid. The below figures compare the flow field around an airfoil with two methods. The stream function and potential function are interrelated to each other by velocity fields. Stream solutions do not require additional steps to show the stream around an object as streamplot is basically doing the same calculations. But from the potential solver, one is required to take derivatives in order to obtain the velocity field and then use streamplot to show the streamline:

![Ekran görüntüsü 2023-10-22 152047](https://github.com/MehmetHakanSari/PotentialSolver/assets/112701635/b062577a-2499-415d-b6b4-458a48888beb)
