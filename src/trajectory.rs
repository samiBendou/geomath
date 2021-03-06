//!
//! A trajectory is a circular buffer containing vectors. It can be seen as the trajectory
//! of a moving point, in this case, the vectors are the position vectors of the point.
//!
//! **Note :** Use the method `push` to add a point to a trajectory.
//! ## Purpose
//! This abstraction aims represent only the last points of a complete trajectory.
//! It's main purpose is to allow drawing a trajectory consistently without allocating
//! new space each time a point is added.
//!
//! ## Conventions
//! * All the trajectories have the same fixed size
//! * Vector space operations can be performed between trajectories
//! * Translation can be performed with a vector using `+` and `-` operators
//! * Scaling can be performed using `*` and `/` operators
//! * When adding a point to a trajectory, the oldest point is erased
//!

use std::fmt;
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign};

use crate::prelude::{Initializer, Reset};
use crate::vector::{Vector2, Vector3, Vector4};
use crate::trajectory::consts::*;

/// 2D trajectory
pub type Trajectory2 = Trajectory<Vector2>;

/// 3D trajectory
pub type Trajectory3 = Trajectory<Vector3>;

/// 4D trajectory
pub type Trajectory4 = Trajectory<Vector4>;


/// constant trajectories
pub mod consts {
    use super::*;
    use crate::vector;


    /// Size of all the trajectories
    pub const TRAJECTORY_SIZE: usize = 256;

    /// 2D zero trajectory
    pub const ZEROS_2: Trajectory2 = Trajectory2 {
        positions: [vector::consts::ZEROS_2; TRAJECTORY_SIZE],
        index: 0,
    };

    /// 3D zero trajectory
    pub const ZEROS_3: Trajectory3 = Trajectory3 {
        positions: [vector::consts::ZEROS_3; TRAJECTORY_SIZE],
        index: 0,
    };

    /// 4D zero trajectory
    pub const ZEROS_4: Trajectory4 = Trajectory4 {
        positions: [vector::consts::ZEROS_4; TRAJECTORY_SIZE],
        index: 0,
    };
}

/// Generic structure and documentation for trajectories
#[derive(Copy, Clone)]
pub struct Trajectory<T> {
    positions: [T; TRAJECTORY_SIZE],
    index: usize,
}

impl<T> From<T> for Trajectory<T> where
    T: Copy + Clone
{
    fn from(position: T) -> Self {
        Trajectory::new([position; TRAJECTORY_SIZE], 0)
    }
}

impl<T> From<[T; TRAJECTORY_SIZE]> for Trajectory<T> where
    T: Copy + Clone {
    fn from(positions: [T; TRAJECTORY_SIZE]) -> Self {
        Trajectory::new(positions, 0)
    }
}

impl<T> Trajectory<T> where
    T: Copy + Clone {

    /// Construct a trajectory with given positions and index of the oldest position in the array
    #[inline]
    pub fn new(positions: [T; TRAJECTORY_SIZE], index: usize) -> Trajectory<T> {
        Trajectory { positions, index }
    }

    /// Add a vector to the trajectory
    #[inline]
    pub fn push(&mut self, position: &T) {
        self.positions[self.index] = *position;
        self.index = self.index_offset(1);
    }

    /// Get the positions vectors of the trajectory as array
    #[inline]
    pub fn positions(&self) -> &[T; TRAJECTORY_SIZE] { &self.positions }

    /// Get the position vector at index. The 0 index is the oldest position vector
    #[inline]
    pub fn position(&self, i: usize) -> &T {
        &self.positions[self.index_offset(i)]
    }

    /// Same as `position` but with a mutable reference
    #[inline]
    pub fn position_mut(&mut self, i: usize) -> &mut T {
        &mut self.positions[self.index_offset(i)]
    }

    /// Get the lastly added position vector
    pub fn last(&self) -> &T {
        &self.positions[self.index_offset(TRAJECTORY_SIZE - 1)]
    }

    /// Same as `last` but with a mutable reference
    pub fn last_mut(&mut self) -> &mut T {
        &mut self.positions[self.index_offset(TRAJECTORY_SIZE - 1)]
    }

    #[inline]
    fn index_offset(&self, i: usize) -> usize {
        (i + self.index) % TRAJECTORY_SIZE
    }
}

impl<T> Initializer for Trajectory<T> where
    T: Initializer + Copy + Clone {
    #[inline]
    fn zeros() -> Self {
        Trajectory::new([T::zeros(); TRAJECTORY_SIZE], 0)
    }

    #[inline]
    fn ones() -> Self {
        Trajectory::new([T::ones(); TRAJECTORY_SIZE], 0)
    }
}

impl<T> Reset<T> for Trajectory<T> where
    T: Reset<T> + Copy + Clone {
    #[inline]
    fn reset0(&mut self) -> &mut Self {
        for position in self.positions.iter_mut() {
            position.reset0();
        }
        self
    }

    #[inline]
    fn reset1(&mut self) -> &mut Self {
        for position in self.positions.iter_mut() {
            position.reset1();
        }
        self
    }

    #[inline]
    fn reset(&mut self, val: &T) -> &mut Self {
        for pos in self.positions.iter_mut() {
            pos.reset(val);
        }
        self
    }
}

impl<T> Debug for Trajectory<T> where
    T: Debug + Copy + Clone {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut buffer = String::new();
        for i in TRAJECTORY_SIZE - 32..TRAJECTORY_SIZE {
            buffer += format!("\n{:?}", self.positions[self.index_offset(i)]).as_str();
        }
        write!(f, "{}", buffer)
    }
}

impl<T> Add<Trajectory<T>> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;
    #[inline]
    fn add(self, rhs: Trajectory<T>) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl<T> Add<T> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;
    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        let mut ret = self;
        ret += rhs;
        ret
    }
}

impl<T> AddAssign<Trajectory<T>> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    #[inline]
    fn add_assign(&mut self, rhs: Trajectory<T>) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[self.index_offset(i)] += rhs.positions[rhs.index_offset(i)];
        }
    }
}

impl<T> AddAssign<T> for Trajectory<T> where
    T: AddAssign<T> + Copy + Clone {
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] += rhs;
        }
    }
}

impl<T> Sub<Trajectory<T>> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;
    #[inline]
    fn sub(self, rhs: Trajectory<T>) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl<T> Sub<T> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    type Output = Trajectory<T>;
    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        let mut ret = self;
        ret -= rhs;
        ret
    }
}

impl<T> SubAssign<Trajectory<T>> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    #[inline]
    fn sub_assign(&mut self, rhs: Trajectory<T>) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[self.index_offset(i)] -= rhs.positions[rhs.index_offset(i)];
        }
    }
}

impl<T> SubAssign<T> for Trajectory<T> where
    T: SubAssign<T> + Copy + Clone {
    #[inline]
    fn sub_assign(&mut self, rhs: T) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] -= rhs;
        }
    }
}

impl<T> Mul<f64> for Trajectory<T> where
    T: MulAssign<f64> + Copy + Clone {
    type Output = Trajectory<T>;
    #[inline]
    fn mul(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl<T> MulAssign<f64> for Trajectory<T> where
    T: MulAssign<f64> + Copy + Clone {
    #[inline]
    fn mul_assign(&mut self, rhs: f64) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] *= rhs;
        }
    }
}

impl<T> Div<f64> for Trajectory<T> where
    T: DivAssign<f64> + Copy + Clone {
    type Output = Trajectory<T>;
    #[inline]
    fn div(self, rhs: f64) -> Self::Output {
        let mut ret = self;
        ret /= rhs;
        ret
    }
}

impl<T> DivAssign<f64> for Trajectory<T> where
    T: DivAssign<f64> + Copy + Clone {
    #[inline]
    fn div_assign(&mut self, rhs: f64) {
        for i in 0..TRAJECTORY_SIZE {
            self.positions[i] /= rhs;
        }
    }
}

impl<T> Index<usize> for Trajectory<T> where
    T: Copy + Clone {
    type Output = T;
    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        &self.positions[self.index_offset(index)]
    }
}

impl<T> IndexMut<usize> for Trajectory<T> where
    T: Copy + Clone {
    #[inline]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.positions[self.index_offset(index)]
    }
}