//!
//! A point encapsulates position, speed and trajectory.
//!
//! It can be seen as a moving point where we store both the kinematic state of the point
//! and it's trajectory.
//!
//! ## Purpose
//! The main goal of the abstraction is to provide a handler for solving ODE of second order
//! such as the 2nd law of Newton.
//! It also provide some affine geometry features.
//!
//! ## Conventions
//! * The trajectory needs to be updated by the calling code
//! * A vector 2N is seen as an upper part containing the position and a lower part containing a speed
//! * Operations between points changes only position and speed
//! * Operations between points of size N and vectors of size N changes only position
//! * Operations between points of size N and vectors of size 2N changes position and speed
//! * Metric space operations between points are perform between their positions
//! * Equality between points checks at both positions and speed
//!
//! **Note :** The trajectory is never updated, you have to call `update_trajectory` method to
//! add the current position value to the trajectory.
//!
use std::fmt;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

use crate::prelude::*;
use crate::trajectory;
use crate::trajectory::Trajectory;
use crate::vector::{self, *};

/// 2D point
pub type Point2 = Point<Vector2>;
/// 3D point
pub type Point3 = Point<Vector3>;
/// 4D point
pub type Point4 = Point<Vector4>;

/// constant points
pub mod consts {
    use crate::vector;
    use super::*;

    /// 2D zero point
    pub const ZEROS_2: Point2 = Point2 {
        position: vector::consts::ZEROS_2,
        speed: vector::consts::ZEROS_2,
        trajectory: trajectory::consts::ZEROS_2,
    };

    /// 3D zero point
    pub const ZEROS_3: Point3 = Point3 {
        position: vector::consts::ZEROS_3,
        speed: vector::consts::ZEROS_3,
        trajectory: trajectory::consts::ZEROS_3,
    };

    /// 4D zero point
    pub const ZEROS_4: Point4 = Point4 {
        position: vector::consts::ZEROS_4,
        speed: vector::consts::ZEROS_4,
        trajectory: trajectory::consts::ZEROS_4,
    };
}

/// Generic structure and documentation for points
#[derive(Copy, Clone)]
pub struct Point<T> {
    pub position: T,
    pub speed: T,
    pub trajectory: Trajectory<T>,
}

impl From<Vector2> for Point2 {
    fn from(vector: Vector2) -> Self {
        Point::new(vector, vector::consts::ZEROS_2)
    }
}

impl From<Vector3> for Point3 {
    fn from(vector: Vector3) -> Self {
        Point::new(vector, vector::consts::ZEROS_3)
    }
}

impl From<Vector4> for Point2 {
    fn from(vector: Vector4) -> Self {
        Point::new(vector.upper(), vector.lower())
    }
}

impl From<Vector6> for Point3 {
    fn from(vector: Vector6) -> Self {
        Point3::new(vector.upper(), vector.lower())
    }
}

impl<T> Point<T> where
    T: Copy + Clone + AddAssign<T> + SubAssign<T> {
    /// Construct a point with given position and speed vectors
    //noinspection RsTypeCheck
    #[inline]
    pub fn new(position: T, speed: T) -> Point<T> {
        Point {
            position,
            speed,
            trajectory: Trajectory::from(position),
        }
    }

    /// Updates trajectory
    ///
    /// The value of `position` is pushed into the trajectory.
    #[inline]
    pub fn update_trajectory(&mut self) -> &mut Self {
        self.trajectory.push(&self.position);
        self
    }

    /// Change the origin
    ///
    /// This method also changes the origin of the trajectory.
    #[inline]
    pub fn reset_origin(&mut self, origin: &Self, old_origin: &Self) -> &mut Self {
        self.position += old_origin.position;
        self.speed += old_origin.speed;
        self.trajectory += old_origin.trajectory;
        self.position -= origin.position;
        self.speed -= origin.speed;
        self.trajectory -= origin.trajectory;
        self
    }
}

impl<T> Initializer for Point<T> where
    T: Initializer + Copy + Clone + AddAssign<T> + SubAssign<T> {
    #[inline]
    fn zeros() -> Self {
        Point::new(T::ones(), T::ones())
    }

    #[inline]
    fn ones() -> Self {
        Point::new(T::ones(), T::ones())
    }
}

impl<T> Reset<Point<T>> for Point<T> where
    T: Reset<T> + Copy + Clone + AddAssign<T> + SubAssign<T> {
    fn reset0(&mut self) -> &mut Self {
        self.position.reset0();
        self.speed.reset0();
        self
    }

    fn reset1(&mut self) -> &mut Self {
        self.position.reset1();
        self.speed.reset1();
        self
    }

    fn reset(&mut self, val: &Point<T>) -> &mut Self {
        self.position = val.position;
        self.speed = val.speed;
        self
    }
}

impl<T> Debug for Point<T> where
    T: Debug + Copy + Clone + AddAssign<T> + SubAssign<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(
            f,
            "position: {:?}\nspeed: {:?}\n",
            self.position,
            self.speed,
        )
    }
}

impl Array<[f64; 4]> for Point2 {
    #[inline]
    fn to_array(&self) -> [f64; 4] {
        [self.position.x, self.position.y, self.speed.x, self.speed.y]
    }

    #[inline]
    fn set_array(&mut self, array: &[f64; 4]) -> &mut Self {
        self.position.x = array[0];
        self.position.y = array[1];
        self.speed.x = array[2];
        self.speed.y = array[3];
        self
    }
}

impl Array<[f64; 6]> for Point3 {
    #[inline]
    fn to_array(&self) -> [f64; 6] {
        [self.position.x, self.position.y, self.position.z, self.speed.x, self.speed.y, self.speed.z]
    }

    #[inline]
    fn set_array(&mut self, array: &[f64; 6]) -> &mut Self {
        self.position.x = array[0];
        self.position.y = array[1];
        self.position.z = array[2];
        self.speed.x = array[3];
        self.speed.y = array[4];
        self.speed.z = array[5];
        self
    }
}

impl Vector<Vector4> for Point2 {
    #[inline]
    fn to_vector(&self) -> Vector4 {
        Vector4::concat(&self.position, &self.speed)
    }

    #[inline]
    fn set_vector(&mut self, vector: &Vector4) -> &mut Self {
        self.position.x = vector.x;
        self.position.y = vector.y;
        self.speed.x = vector.z;
        self.speed.y = vector.w;
        self
    }
}

impl Vector<Vector6> for Point3 {
    #[inline]
    fn to_vector(&self) -> Vector6 {
        Vector6::concat(&self.position, &self.speed)
    }

    #[inline]
    fn set_vector(&mut self, vector: &Vector6) -> &mut Self {
        self.position.x = vector.x;
        self.position.y = vector.y;
        self.position.z = vector.z;
        self.speed.x = vector.u;
        self.speed.y = vector.v;
        self.speed.z = vector.w;
        self
    }
}

impl<T> Metric for Point<T> where
    T: Metric + Copy + Clone + AddAssign<T> + SubAssign<T> {
    #[inline]
    fn dot(&self, other: &Self) -> f64 {
        self.position.dot(&other.position)
    }

    #[inline]
    fn distance2(&self, other: &Self) -> f64 {
        self.position.distance2(&other.position)
    }

    #[inline]
    fn distance(&self, other: &Self) -> f64 {
        self.position.distance(&other.position)
    }

    #[inline]
    fn magnitude2(&self) -> f64 {
        self.position.magnitude2()
    }

    #[inline]
    fn magnitude(&self) -> f64 {
        self.position.magnitude()
    }

    fn set_normalized(&mut self) -> &mut Self {
        self.position.set_normalized();
        self
    }
}

impl<T> PartialEq for Point<T> where
    T: PartialEq<T> + Copy + Clone + AddAssign<T> + SubAssign<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.position == other.position && self.speed == other.speed
    }

    #[inline]
    fn ne(&self, other: &Self) -> bool {
        self.position != other.position || self.speed != other.speed
    }
}

impl<T> Add<T> for Point<T> where
    T: AddAssign<T> + Copy + Clone + SubAssign<T> {
    type Output = Point<T>;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl<T> Add<Point<T>> for Point<T> where
    T: AddAssign<T> + Copy + Clone + SubAssign<T> {
    type Output = Point<T>;

    #[inline]
    fn add(self, rhs: Point<T>) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl Add<Vector4> for Point2 {
    type Output = Point2;

    #[inline]
    fn add(self, rhs: Vector4) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl Add<Vector6> for Point3 {
    type Output = Point3;

    #[inline]
    fn add(self, rhs: Vector6) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl<T> AddAssign<T> for Point<T> where
    T: AddAssign<T> + Copy + Clone + SubAssign<T> {
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        self.position += rhs;
    }
}

impl<T> AddAssign<Point<T>> for Point<T> where
    T: AddAssign<T> + Copy + Clone + SubAssign<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Point<T>) {
        self.position += rhs.position;
        self.speed += rhs.speed;
    }
}

impl AddAssign<Vector4> for Point2 {
    #[inline]
    fn add_assign(&mut self, rhs: Vector4) {
        self.position.x += rhs.x;
        self.position.y += rhs.y;
        self.speed.x += rhs.z;
        self.speed.y += rhs.w;
    }
}

impl AddAssign<Vector6> for Point3 {
    #[inline]
    fn add_assign(&mut self, rhs: Vector6) {
        self.position.x += rhs.x;
        self.position.y += rhs.y;
        self.position.z += rhs.z;
        self.speed.x += rhs.u;
        self.speed.y += rhs.v;
        self.speed.z += rhs.w;
    }
}

impl<T> Sub<T> for Point<T> where
    T: SubAssign<T> + Copy + Clone + AddAssign<T> {
    type Output = Point<T>;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl<T> Sub<Point<T>> for Point<T> where
    T: SubAssign<T> + Copy + Clone + AddAssign<T> {
    type Output = Point<T>;

    #[inline]
    fn sub(self, rhs: Point<T>) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl Sub<Vector4> for Point2 {
    type Output = Point2;

    #[inline]
    fn sub(self, rhs: Vector4) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl Sub<Vector6> for Point3 {
    type Output = Point3;

    #[inline]
    fn sub(self, rhs: Vector6) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl<T> SubAssign<T> for Point<T> where
    T: SubAssign<T> + Copy + Clone + AddAssign<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: T) {
        self.position -= rhs;
    }
}

impl<T> SubAssign<Point<T>> for Point<T> where
    T: SubAssign<T> + Copy + Clone + AddAssign<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Point<T>) {
        self.position -= rhs.position;
        self.speed -= rhs.speed;
    }
}

impl SubAssign<Vector4> for Point2 {
    #[inline]
    fn sub_assign(&mut self, rhs: Vector4) {
        self.position.x -= rhs.x;
        self.position.y -= rhs.y;
        self.speed.x -= rhs.z;
        self.speed.y -= rhs.w;
    }
}

impl SubAssign<Vector6> for Point3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Vector6) {
        self.position.x -= rhs.x;
        self.position.y -= rhs.y;
        self.position.z -= rhs.z;
        self.speed.x -= rhs.u;
        self.speed.y -= rhs.v;
        self.speed.z -= rhs.w;
    }
}

impl<T> Mul<f64> for Point<T> where
    T: MulAssign<f64> + Copy + Clone + AddAssign<T> + SubAssign<T> {
    type Output = Point<T>;

    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        let mut output = self;
        output *= rhs;
        output
    }
}

impl<T> MulAssign<f64> for Point<T> where
    T: MulAssign<f64> + Copy + Clone + AddAssign<T> + SubAssign<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        self.position *= rhs;
        self.speed *= rhs;
    }
}

impl<T> Div<f64> for Point<T> where
    T: DivAssign<f64> + Copy + Clone + AddAssign<T> + SubAssign<T> {
    type Output = Point<T>;

    #[inline]
    fn div(self, rhs: f64) -> Self::Output {
        let mut output = self;
        output /= rhs;
        output
    }
}

impl<T> DivAssign<f64> for Point<T> where
    T: DivAssign<f64> + Copy + Clone + AddAssign<T> + SubAssign<T> {
    #[inline]
    fn div_assign(&mut self, rhs: f64) {
        self.position /= rhs;
        self.speed /= rhs;
    }
}

