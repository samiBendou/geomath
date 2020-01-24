//!
//! The prelude specifies most of the features provided by the objects of the framework.
//! This goes from math properties to array conversion and vectors concatenation.
//!

/// Construct objects either filled with ones or zeros
pub trait Initializer {
    /// Object filled with zeros
    fn zeros() -> Self;
    /// Object filled with ones
    fn ones() -> Self;
}

/// Fill an object with a value
pub trait Reset<T> {
    /// Fills an object with zeros
    fn reset0(&mut self) -> &mut Self;
    /// Fills a object with ones
    fn reset1(&mut self) -> &mut Self;

    /// Fills a object with a given value
    fn reset(&mut self, val: &T) -> &mut Self;
}

/// Set an object from an array and convert an object onto an array
pub trait Array<T> {
    /// Convert the object to array
    fn to_array(&self) -> T;
    /// Set the object from array
    fn set_array(&mut self, arr: &T) -> &mut Self;
}

/// Set an object from a vector and convert an object onto a vector
pub trait Vector<T> {
    /// Convert the object to vector
    fn to_vector(&self) -> T;
    /// Set the object from vector
    fn set_vector(&mut self, arr: &T) -> &mut Self;
}

/// Split an object onto two part that can be concatenated or accessed separately
///
/// It aims to be the upper and lower parts of a vector or the upper diagonal and lower diagonal parts of a matrix.
pub trait Split<T> {
    /// Get the two parts of the split object
    fn split(&self) -> [T; 2];
    /// Concatenates two part to construct an object
    fn concat(lhs: &T, rhs: &T) -> Self;
    /// Get the first part of the object
    fn upper(&self) -> T;
    /// Get the second part of the object
    fn lower(&self) -> T;
    /// Set the first part of the object
    fn set_upper(&mut self, val: &T) -> &mut Self;
    /// Set the second part of the object
    fn set_lower(&mut self, val: &T) -> &mut Self;
}

/// Operations between objects of a metric space
///
/// A metric space is considered as a space where it exists a norm based on dot product.
/// It can be the metric operations between two vectors or matrices such as the distance or the magnitude.
///
/// **Note :** For matrices and vectors you can also use operators : `%` for distance, `|` for dot product, `!` for magnitude
pub trait Metric where
Self: Copy {
    /// Dot product of the two objects
    fn dot(&self, other: &Self) -> f64;
    /// Squared distance between the two objects
    fn distance2(&self, other: &Self) -> f64;
    /// Distance between the two objects
    fn distance(&self, other: &Self) -> f64;
    /// Squared magnitude of an object
    fn magnitude2(&self) -> f64;
    /// Magnitude of an object
    fn magnitude(&self) -> f64;
    /// Get the normalized vector, ie. vector with same direction and magnitude 1
    fn normalized(&self) -> Self {
        let mut ret = *self;
        ret.set_normalized();
        ret
    }
    /// Normalizes the vector, ie. sets magnitude to 1 without changing direction
    fn set_normalized(&mut self) -> &mut Self;
}
/// Operations related to the angle of two vectors
///
/// Do not use theses features with zero vector since divisions with the magnitude are performed.
pub trait Angle where
Self: Metric {
    /// Cosine of the angle between the two vectors, non-oriented
    fn cos(&self, rhs: &Self) -> f64 {
        self.dot(rhs) / (self.magnitude() * rhs.magnitude())
    }
    /// Sine of the angle between the two vectors, non-oriented
    fn sin(&self, rhs: &Self) -> f64 {
        self.area(rhs) / (self.magnitude() * rhs.magnitude())
    }
    /// Area of the parallelepiped formed by two vectors, non-oriented
    fn area(&self, rhs: &Self) -> f64;
    /// Angle between two vectors, non-oriented
    fn angle(&self, rhs: &Self) -> f64 {
        self.cos(rhs).acos()
    }
}

/// Stable cross product of two vectors
pub trait Cross where
Self: Copy + Clone {
    /// Get the cross product between two vectors
    fn cross(&self, rhs: &Self) -> Self {
        let mut ret = *self;
        ret.set_cross(rhs);
        ret
    }
    /// Set the cross product between two vectors
    fn set_cross(&mut self, rhs: &Self) -> &mut Self;
}

/// Interpolations between two objects
pub trait Interpolation {
    /// Get linear interpolation
    fn lerp(&self, other: &Self, s: f64) -> Self where
        Self: Copy + Clone {
        let mut ret = *self;
        ret.set_lerp(other, s);
        ret
    }
    /// Get cubic Hermite's interpolation, ie. with two tangent values
    fn herp(&mut self, other: &Self, other1: &Self, other2: &Self, s: f64) -> Self where
        Self: Copy + Clone {
        let mut ret = *self;
        ret.set_herp(other, other1, other2, s);
        ret
    }

    /// Get cubic Bezier's interpolation, ie. with two control points
    fn berp(&mut self, other: &Self, other1: &Self, other2: &Self, s: f64) -> Self where
        Self: Copy + Clone {
        let mut ret = *self;
        ret.set_berp(other, other1, other2, s);
        ret
    }

    /// Set the linear interpolation
    fn set_lerp(&mut self, other: &Self, s: f64) -> &mut Self;
    /// Set the Hermite's interpolation
    fn set_herp(&mut self, other: &Self, other1: &Self, other2: &Self, s: f64) -> &mut Self;
    /// Set the Bezier's interpolation
    fn set_berp(&mut self, other: &Self, other1: &Self, other2: &Self, s: f64) -> &mut Self;
}

/// Access the matrix as an array of rows ordered from left to right
pub trait Rows<T> where
    Self: std::marker::Sized + Initializer {
    /// Get a matrix from rows array
    fn from_rows(rows: &T) -> Self {
        let mut ret = Self::zeros();
        ret.set_rows(rows);
        ret
    }
    /// Get the rows of a matrix
    fn rows(&self) -> T;

    /// Set the rows of a matrix
    fn set_rows(&mut self, rows: &T) -> &mut Self;
}

/// Square matrices algebra
pub trait Algebra<T> where
    Self: std::marker::Sized + Copy + Clone,
    T: std::marker::Sized + Copy + Clone {

    /// Get the determinant of the matrix
    fn determinant(&self) -> f64;

    /// Get the inverse matrix
    fn inverse(&self) -> Self {
        let mut ret = *self;
        ret.set_inverse();
        ret
    }

    /// Get the transposed matrix
    fn transposed(&self) -> Self {
        let mut ret = *self;
        ret.set_transposed();
        ret
    }

    /// Get the adjugate matrix
    fn adjugate(&self) -> Self {
        let mut ret = *self;
        ret.set_adjugate();
        ret
    }

    /// Set inverse matrix
    fn set_inverse(&mut self) -> &mut Self;

    /// Set transposed matrix
    fn set_transposed(&mut self) -> &mut Self;

    /// Set adjugate matrix
    fn set_adjugate(&mut self) -> &mut Self;
}

/// Coordinates accessors
///
/// Get and Set coordinates of 2D and 3D vectors.
///
/// ## Conventions
/// A coherent naming convention between polar, cylindrical and spherical coordinates has been established.
///
/// ### Coordinates naming
/// ![coordinates](https://github.com/samiBendou/geomath/raw/dev/media/coordinates_diagram.png)
///
/// **3D coordinates representation in geomath**
///
/// The four basic coordinates systems are represented
/// * (x, y, z) the cartesian coordinates
/// * (rho, phi, z) the polar/cylindrical coordinates
/// * (r, phi, theta) the spherical coordinates
///
/// The traits are implemented such that there is no code repetition between the systems that have coordinates
/// in common. For example, since the angle phi is common to polar and spherical coordinates
/// systems it will be only defined in the `Polar` trait.
///
/// **Notes :**
///
/// * The value of phi remains the same in both polar/cylindrical and spherical coordinates systems
/// * The value of r in the schema above is denoted `radius` in the code
/// * The z coordinate is common to cartesian and cylindrical coordinates systems
///
/// ### Local basis
/// You can generate vectors of the local basis of each coordinate system, the convention is described
/// in the schemas bellow.
///
/// ![coordinates](https://github.com/samiBendou/geomath/raw/dev/media/local_%20basis_diagram.png)
///
pub mod coordinates {

    /// Polar coordinates accessors
    pub trait Polar {
        /// Get a vector from polar coordinates
        fn from_polar(rho: f64, phi: f64) -> Self;

        /// Set a vector from polar coordinates
        fn set_polar(&mut self, rho: f64, phi: f64) -> &mut Self;

        /// Get a unit polar radial vector from angle
        fn unit_rho(phi: f64) -> Self;

        /// Get a unit tangent/prograde vector from angle
        fn unit_phi(phi: f64) -> Self;

        /// Get the rho coordinate
        fn rho(&self) -> f64;

        /// Get the phi coordinate
        fn phi(&self) -> f64;

        /// Set the rho coordinate
        fn set_rho(&mut self, rho: f64) -> &mut Self;

        /// Set the phi coordinate
        fn set_phi(&mut self, phi: f64) -> &mut Self;
    }

    /// Cylindrical coordinates accessors
    pub trait Cylindrical {
        /// Get a vector from cylindrical coordinates
        fn from_cylindrical(rho: f64, phi: f64, z: f64) -> Self;

        /// Set a vector from cylindrical coordinates
        fn set_cylindrical(&mut self, rho: f64, phi: f64, z: f64) -> &mut Self;
    }

    /// Spherical coordinates accessors
    pub trait Spherical {
        /// Get a vector from spherical coordinates
        fn from_spherical(radius: f64, phi: f64, theta: f64) -> Self;
        /// Set a vector from spherical coordinates
        fn set_spherical(&mut self, radius: f64, phi: f64, theta: f64) -> &mut Self;
        /// Get a unit spherical radial vector from angles
        fn unit_radius(phi: f64, theta: f64) -> Self;
        /// Get a unit tangent/normal vector from angles
        fn unit_theta(phi: f64, theta: f64) -> Self;
        /// Get the theta coordinate
        fn theta(&self) -> f64;
        /// Set the theta coordinates
        fn set_theta(&mut self, theta: f64) -> &mut Self;
    }

    /// Homogeneous coordinates transforms
    pub trait Homogeneous<T> {
        /// Get a inhomogeneous vector from an homogeneous vector
        ///
        /// The inhomogeneous vector has the last component removed and the others
        /// divided by the value of the last component of the homogeneous vector.
        fn from_homogeneous(vector: &T) -> Self;

        /// Get a vector transformed to an homogeneous vector
        ///
        /// It is the same vector but with a component with the value 1 added.
        fn to_homogeneous(&self) -> T;
    }
}

/// Transform matrix
///
/// Set and generate transform matrices.
/// The transform matrices can represent common 2D and 3D geometrical operations such as translation, rotation, ...
pub mod transforms {
    use crate::prelude::Initializer;
    use crate::vector::Vector3;

    /// Translation matrix
    pub trait Translation<T> where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        /// Get a translation matrix from translation vector
        fn from_translation(vector: &T) -> Self {
            let mut ret = Self::zeros();
            ret.set_translation(vector);
            ret
        }

        /// Set a translation matrix from translation vector
        fn set_translation(&mut self, vector: &T) -> &mut Self;
    }

    /// Rigid body matrix
    ///
    /// A rigid body transform is the combination of a rotation and a translation.
    pub trait Rigid<U, T> where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        /// Get a rigid body matrix from rotation matrix and translation vector
        fn from_rigid(rotation: &U, vector: &T) -> Self {
            let mut ret = Self::zeros();
            ret.set_rigid(rotation, vector);
            ret
        }

        /// Set a rigid body matrix from the given rotation matrix and translation vector
        fn set_rigid(&mut self, rotation: &U, vector: &T) -> &mut Self;
    }

    /// Similarity matrix
    ///
    /// A similarity is the combination of a scaling, a rotation and a translation.
    pub trait Similarity<U, T> where
        Self: std::marker::Sized + Copy + Clone + Initializer {
        /// Get a similarity matrix from scale factor, rotation matrix and translation vector
        fn from_similarity(scale: f64, rotation: &U, vector: &T) -> Self {
            let mut ret = Self::zeros();
            ret.set_similarity(scale, rotation, vector);
            ret
        }

        /// Set a similarity matrix from scale factor, rotation matrix and translation vector
        fn set_similarity(&mut self, scale: f64, rotation: &U, vector: &T) -> &mut Self;
    }

    /// 2D rotations matrix
    pub trait Rotation2 where
        Self: std::marker::Sized + Copy + Clone + Initializer {

        /// Get a rotation matrix from angle
        fn from_rotation(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation(angle);
            ret
        }

        /// Set a rotation matrix from angle
        fn set_rotation(&mut self, angle: f64) -> &mut Self;
    }

    /// 3D rotations matrix
    pub trait Rotation3 where
        Self: std::marker::Sized + Copy + Clone + Initializer {

        /// Get a rotation matrix from axis and angle
        fn from_rotation(angle: f64, axis: &Vector3) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation(angle, axis);
            ret
        }

        /// Get a rotation matrix around x-axis from given angle
        fn from_rotation_x(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation_x(angle);
            ret
        }

        /// Get a rotation matrix around y-axis from given angle
        fn from_rotation_y(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation_y(angle);
            ret
        }

        /// Get a rotation matrix around z-axis from given angle
        fn from_rotation_z(angle: f64) -> Self {
            let mut ret = Self::zeros();
            ret.set_rotation_z(angle);
            ret
        }

        /// Set a rotation matrix from axis and angle
        fn set_rotation(&mut self, angle: f64, axis: &Vector3) -> &mut Self;

        /// Set a rotation matrix around x-axis from given angle
        fn set_rotation_x(&mut self, angle: f64) -> &mut Self;

        /// Set a rotation matrix around y-axis from given angle
        fn set_rotation_y(&mut self, angle: f64) -> &mut Self;

        /// Set a rotation matrix around z-axis from given angle
        fn set_rotation_z(&mut self, angle: f64) -> &mut Self;
    }
}