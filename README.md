# geomath
### A framework for simulation and computer graphics

## Brief
geomath is a little framework that implements most of the common math features needed 
when dealing with numerical simulation or computer graphics. 
It focuses mainly on algebra and geometry and provide a simple API while showing high performance 
thanks to many optimization that are allowed only in a 2D-3D-4D context.

## Features
- Stack allocated Matrix and Vectors with common algebra and operators
- All common 3D and 4D transforms (perspective, rotations, ...)
- Coordinates manipulation (polar, cylindrical, spherical, ...)
- Point's kinematics

## Usage
The framework is yet not documented, it still in development. However it won't be to difficult to use it if you are familiar with numpy or Matlab.
However if you need examples take a look at :
- [this example project](https://github.com/samiBendou/nbodies)
- [this library](https://github.com/samiBendou/dynamics)

### Why another maths framework
nalgebra and cgmath do the job well, however if you prefer use a more concise syntax this framework provides methods
(that can be chained) and generators to compute very quickly a considerable range of transform matrix, vectors, ...

It's really made to make maths easier when coding, the goal is to reach the feeling of Matlab.

Furthermore, the implementation is the result of the previous optimizations found during the two other math framework that
I published. 

It runs very fast and provide an API that shows clearly where computation is expensive so you can adopt the best patterns
to reach the desired performance.