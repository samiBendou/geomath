//!
//! geomath is a general purpose maths framework that aims to provide efficient real-time tools for the
//! following domains:
//!
//! * Linear algebra
//! * Computational Geometry
//! * Computer Graphics and Vision
//! * Physics and Kinematics simulation
//! * Numerical simulation
//!
//! It relies on a vast API that consists on various traits each specifying a math property.
//! It provides most of the common abstraction used in the above domains :
//! * 2D, 3D, 4D matrices, trajectories and points
//! * 2D, 3D, 4D, 6D vectors
//!
//! The implementation benefits from many optimizations allowed when working with small dimensions :
//! * Straight forward code
//! * Hand-optimized operations
//! * Fully stack allocated
//!
//! However, it still uncomplete and other features and improvements are planned,
//! take a look at the [GitHub developement branch](https://github.com/samiBendou/geomath)  for more details.
//!
//! # Get started
//! If you already familiar with linear algebra libraries, you can easily get started since the
//! API is really straight forward,
//! basic vector space operations are implemented as operators and do not need any traits.
//! Other operations can be imported with the prelude.
//!
//! ## Basic principles
//! Some simple conventions have been choosed either to gain clarity or performance :
//! * Vectors are columns matrices when operating with a matrix
//! * Vectors cartesian coordinates are accessible using `vector.i` syntax with `i` the name of the coordinate
//! * Matrices components can be accessed using `matrix.ij` where `i` is the row index and `j` the column index
//! * Rows and columns are indexed with `x`, `y`, `z`, ... starting with `x`
//! * Vectors of size N can be left multiplied by matrices of size N+1
//!
//! **Note :** When left-multiplying a vector of size N with a matrix of size N+1 the vector is seen as an homogeneous
//! vector (ie. vector of size N) with added component at 1.
//! It's a way to perfrom efficient multiplications with transformation matrices.
//!
//! ### Example
//! ```
//! use geomath::matrix;
//! let m = matrix::consts::EYE_2;
//! assert_eq!(m.xx, 1.);
//! assert_eq!(m.xy, 0.);
//! assert_eq!(m.yx, 0.);
//! assert_eq!(m.yy, 1.);
//! ```
//!
//! ### Multiply assign
//! The multiply assign operator between matrices and vectors reverses the order of the operands
//! which means `u *= m` computes in fact `u = m * u`. This allows faster assign multiply implementation.
//!
//! ## The prelude
//! As said above, the prelude contains the traits specifications, it is split in a root and 2 modules :
//! ```
//! use geomath::prelude::*; // most of the traits
//! use geomath::prelude::coordinates::*; // coordinates systems and transformations
//! use geomath::prelude::transforms::*; // matrix transforms (translations, rotations, ...)
//! ```
//!
//! ## Constants
//! Instead of providing static generators for unit vectors and identity matrices, these
//! are hard coded as constants :
//! ```
//! use geomath::vector;
//! let u = vector::consts::EX_3; // unit vector in the positive x (1, 0, 0)
//! let v = vector::consts::N_EZ_3; // unit vector in the negative z direction (0, 0, -1)
//! let w = vector::consts::ZEROS_3; // zero vector of 3D space (0, 0, 0)
//! ```
//! **Note :** All the modules have a `consts` submodule.
//!
//! However if you are more familiar with the generators, you can still do :
//! ```
//! use geomath::prelude::Initializer;
//! use geomath::vector::Vector3;
//! let u = Vector3::zeros(); // zero vector of 3D space (0, 0, 0)
//! let v = Vector3::ones(); // one vector of 3D space (1, 1, 1)
//! ```
//!
//! ## Short constructors
//! To construct a vector a matrix, you can use short constructors :
//! ```
//! use geomath::vector::vec3;
//! use geomath::matrix::mat2;
//! let u = vec3(2., 3., 5.);
//! let m = mat2(1., 2., 3., 4.); // ordered by rows, from left to right
//! ```
//!
//! ## Metric space features
//! The metric space fetuares for vectors and matrices can be either computed using methods or operators :
//! ```
//! use geomath::vector;
//! use geomath::prelude::*;
//! let u = vector::consts::EX_3;
//! let v = vector::consts::N_EZ_3;
//! assert_eq!(u % v, u.distance(&v));
//! assert_eq!(u | v, u.dot(&v));
//! assert_eq!(!u, u.magnitude());
//! ```

/// Vector structures and documentation
pub mod vector;
/// Matrix structures and documentation
pub mod matrix;
/// Point structures and documentation
pub mod point;
/// Trajectory structures and documentation
pub mod trajectory;
/// Prelude traits and documentation
pub mod prelude;
/// Macros mentioned bellow
pub mod macros;

