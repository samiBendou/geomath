
# geomath
[![crates](https://img.shields.io/crates/v/geomath?style=flat-square)](https://crates.io/crates/geomath)
[![documentation](https://img.shields.io/static/v1?label=docs.rs&message=documentation&color=blue&style=flat-square)](https://docs.rs/geomath)
[![license](https://img.shields.io/crates/l/geomath?style=flat-square)](https://github.com/samiBendou/geomath/blob/master/LICENSE)

### Simulation, Graphics, Geometry

## Brief
A framework that aims to provide a complete 2D-3D-4D maths toolbox for the Rust language.

It's general purpose and exposes a vast and simple API while showing high performance 
thanks to many optimization that are allowed only in a 2D-3D-4D context.

## Features
- Stack allocated matrices and vectors with common algebra and transforms
- All common 2D, 3D and 4D transforms (rigid body, rotations, ...)
- Coordinates manipulation (polar, cylindrical, spherical, ...)
- Point's kinematics and trajectory representation
- Documentation and examples

## Usage

The [documentation](https://docs.rs/geomath) contains a short introductions for each module to easily get started with the framework.

However if you need examples take a look at :
- [this example project](https://github.com/samiBendou/nbodies)
- [this library](https://github.com/samiBendou/dynamics)

### Why another maths framework

Instead of providing tools only for linear algebra or computer graphics, this framework is made to cover all theses usages at
the same time.

It's really made to make maths easier when coding, the goal is to reach the feeling of Matlab and to provide a general purpose
 2D, 3D and 4D maths toolbox.

Furthermore, the implementation is the result of the previous optimizations found during the two other math framework that
I published on npm [meca3](https://www.npmjs.com/package/meca3) and [space3](https://www.npmjs.com/package/space3). 

It runs very fast and provide an API that shows clearly where computation is expensive so you can adopt the best patterns
to reach the desired performance.

## Contribution

The framework still largely incomplete, more features are scheduled to be implemented,
some optimizations can still be perform and more unit testing is required.

Any suggestion or issues are welcomed, be free to contribute the geomath project and improve the Rust math experience.


