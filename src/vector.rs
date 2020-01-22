use std::fmt::{Debug, Error, Formatter};
use std::ops::{
    Add, AddAssign,
    BitOr, Div,
    DivAssign,
    Mul,
    MulAssign,
    Neg, Not,
    Rem, Sub,
    SubAssign,
};
use crate::common::*;
use crate::common::coordinates::Polar;
use crate::matrix::{Matrix2, Matrix3, Matrix4};
use crate::point::Point2;
use crate::{impl_vector, impl_debug_vector};

pub mod consts {
    use super::*;

    pub const EX_2: Vector2 = Vector2 { x: 1., y: 0. };
    pub const N_EX_2: Vector2 = Vector2 { x: -1., y: 0. };
    pub const EY_2: Vector2 = Vector2 { x: 0., y: 1. };
    pub const N_EY_2: Vector2 = Vector2 { x: 0., y: -1. };
    pub const ZEROS_2: Vector2 = Vector2 { x: 0., y: 0. };
    pub const ONES_2: Vector2 = Vector2 { x: 1., y: 1. };

    pub const EX_3: Vector3 = Vector3 { x: 1., y: 0., z: 0. };
    pub const N_EX_3: Vector3 = Vector3 { x: -1., y: 0., z: 0. };
    pub const EY_3: Vector3 = Vector3 { x: 0., y: 1., z: 0. };
    pub const N_EY_3: Vector3 = Vector3 { x: 0., y: -1., z: 0. };
    pub const EZ_3: Vector3 = Vector3 { x: 0., y: 0., z: 1. };
    pub const N_EZ_3: Vector3 = Vector3 { x: 0., y: 0., z: -1. };
    pub const ZEROS_3: Vector3 = Vector3 { x: 0., y: 0., z: 0. };
    pub const ONES_3: Vector3 = Vector3 { x: 1., y: 1., z: 1. };

    pub const EX_4: Vector4 = Vector4 { x: 1., y: 0., z: 0., w: 0. };
    pub const N_EX_4: Vector4 = Vector4 { x: -1., y: 0., z: 0., w: 0. };
    pub const EY_4: Vector4 = Vector4 { x: 0., y: 1., z: 0., w: 0. };
    pub const N_EY_4: Vector4 = Vector4 { x: 0., y: -1., z: 0., w: 0. };
    pub const EZ_4: Vector4 = Vector4 { x: 0., y: 0., z: 1., w: 0. };
    pub const N_EZ_4: Vector4 = Vector4 { x: 0., y: 0., z: -1., w: 0. };
    pub const EW_4: Vector4 = Vector4 { x: 0., y: 0., z: 0., w: 1. };
    pub const N_EW_4: Vector4 = Vector4 { x: 0., y: 0., z: -0., w: -1. };
    pub const ZEROS_4: Vector4 = Vector4 { x: 0., y: 0., z: 0., w: 0. };
    pub const ONES_4: Vector4 = Vector4 { x: 1., y: 1., z: 1., w: 1. };
}

#[derive(Copy, Clone)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64,
}

#[derive(Copy, Clone)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Copy, Clone)]
pub struct Vector4 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

#[derive(Copy, Clone)]
pub struct Vector6 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub u: f64,
    pub v: f64,
    pub w: f64,
}

impl From<[f64; 2]> for Vector2 {
    fn from(array: [f64; 2]) -> Self {
        Vector2::new(array[0], array[1])
    }
}

impl From<[f64; 3]> for Vector3 {
    fn from(array: [f64; 3]) -> Self {
        Vector3::new(array[0], array[1], array[2])
    }
}

impl From<[f64; 4]> for Vector4 {
    fn from(array: [f64; 4]) -> Self {
        Vector4::new(array[0], array[1], array[2], array[3])
    }
}

impl From<[f64; 6]> for Vector6 {
    fn from(array: [f64; 6]) -> Self {
        Vector6::new(array[0], array[1], array[2], array[3], array[4], array[5])
    }
}

impl From<Point2> for Vector4 {
    fn from(point: Point2) -> Self {
        Vector4::new(point.position.x, point.position.y, point.speed.x, point.speed.y)
    }
}

impl_vector!(Vector2 {x, y}, 2);
impl_vector!(Vector3 {x, y, z}, 3);
impl_vector!(Vector4 {x, y, z, w}, 4);
impl_vector!(Vector6 {x, y, z, u, v, w}, 6);

impl_debug_vector!(Vector2 {x, y});
impl_debug_vector!(Vector3 {x, y, z});
impl_debug_vector!(Vector4 {x, y, z, w});
impl_debug_vector!(Vector6 {x, y, z, u, v, w});

impl MulAssign<Matrix2> for Vector2 {
    fn mul_assign(&mut self, rhs: Matrix2) {
        let x = self.x;
        let y = self.y;
        self.x = rhs.xx * x + rhs.xy * y;
        self.y = rhs.yx * x + rhs.yy * y;
    }
}

impl MulAssign<Matrix3> for Vector3 {
    fn mul_assign(&mut self, rhs: Matrix3) {
        let x = self.x;
        let y = self.y;
        let z = self.z;
        self.x = rhs.xx * x + rhs.xy * y + rhs.xz * z;
        self.y = rhs.yx * x + rhs.yy * y + rhs.yz * z;
        self.z = rhs.zx * x + rhs.zy * y + rhs.zz * z;
    }
}

impl MulAssign<Matrix3> for Vector2 {
    fn mul_assign(&mut self, rhs: Matrix3) {
        let x = self.x;
        let y = self.y;
        self.x = rhs.xx * x + rhs.xy * y + rhs.xz;
        self.y = rhs.yx * x + rhs.yy * y + rhs.yz;
    }
}

impl MulAssign<Matrix4> for Vector4 {
    fn mul_assign(&mut self, rhs: Matrix4) {
        let x = self.x;
        let y = self.y;
        let z = self.z;
        let w = self.w;
        self.x = rhs.xx * x + rhs.xy * y + rhs.xz * z + rhs.xw * w;
        self.y = rhs.yx * x + rhs.yy * y + rhs.yz * z + rhs.yw * w;
        self.z = rhs.zx * x + rhs.zy * y + rhs.zz * z + rhs.zw * w;
        self.w = rhs.wx * x + rhs.wy * y + rhs.wz * z + rhs.ww * w;
    }
}

impl MulAssign<Matrix4> for Vector3 {
    fn mul_assign(&mut self, rhs: Matrix4) {
        let x = self.x;
        let y = self.y;
        let z = self.z;
        self.x = rhs.xx * x + rhs.xy * y + rhs.xz * z + rhs.xw;
        self.y = rhs.yx * x + rhs.yy * y + rhs.yz * z + rhs.yw;
        self.z = rhs.zx * x + rhs.zy * y + rhs.zz * z + rhs.zw;
    }
}

impl Angle for Vector2 {
    fn area(&self, rhs: &Self) -> f64 {
        self.x * rhs.y - self.y * rhs.x
    }
}

impl Angle for Vector3 {
    fn area(&self, rhs: &Self) -> f64 {
        self.cross(rhs).magnitude()
    }
}

impl Cross for Vector3 {
    fn set_cross(&mut self, rhs: &Self) -> &mut Self {
        let x = self.x;
        let y = self.y;
        let z = self.z;
        self.x = y * rhs.z - z * rhs.y;
        self.y = z * rhs.x - x * rhs.z;
        self.z = x * rhs.y - y * rhs.x;
        self
    }
}

impl coordinates::Polar for Vector2 {
    #[inline]
    fn from_polar(rho: f64, phi: f64) -> Self {
        Vector2 { x: rho * phi.cos(), y: rho * phi.sin() }
    }

    #[inline]
    fn set_polar(&mut self, rho: f64, phi: f64) -> &mut Self {
        self.x = rho * phi.cos();
        self.y = rho * phi.sin();
        self
    }

    #[inline]
    fn unit_rho(ang: f64) -> Self {
        Vector2 { x: ang.cos(), y: ang.sin() }
    }

    #[inline]
    fn unit_phi(ang: f64) -> Self {
        Vector2 { x: -ang.sin(), y: ang.cos() }
    }

    #[inline]
    fn rho(&self) -> f64 {
        self.magnitude()
    }

    #[inline]
    fn phi(&self) -> f64 {
        (self.y).atan2(self.x)
    }

    #[inline]
    fn set_rho(&mut self, rho: f64) -> &mut Self {
        self.set_polar(rho, self.phi())
    }

    #[inline]
    fn set_phi(&mut self, phi: f64) -> &mut Self {
        self.set_polar(self.rho(), phi)
    }
}

impl coordinates::Polar for Vector3 {
    #[inline]
    fn from_polar(rho: f64, phi: f64) -> Self {
        Vector3 { x: rho * phi.cos(), y: rho * phi.sin(), z: 0. }
    }

    #[inline]
    fn set_polar(&mut self, rho: f64, phi: f64) -> &mut Self {
        self.x = rho * phi.cos();
        self.y = rho * phi.sin();
        self
    }

    #[inline]
    fn unit_rho(phi: f64) -> Self {
        Vector3 { x: phi.cos(), y: phi.sin(), z: 0. }
    }

    #[inline]
    fn unit_phi(phi: f64) -> Self {
        Vector3 { x: -phi.sin(), y: phi.cos(), z: 0. }
    }

    #[inline]
    fn rho(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    #[inline]
    fn phi(&self) -> f64 {
        (self.y).atan2(self.x)
    }

    #[inline]
    fn set_rho(&mut self, rho: f64) -> &mut Self {
        self.set_polar(rho, self.phi())
    }

    #[inline]
    fn set_phi(&mut self, phi: f64) -> &mut Self {
        self.set_polar(self.rho(), phi)
    }
}

impl coordinates::Cylindrical for Vector3 {
    #[inline]
    fn from_cylindrical(rho: f64, phi: f64, z: f64) -> Self {
        Vector3 { x: rho * phi.cos(), y: rho * phi.sin(), z }
    }

    #[inline]
    fn set_cylindrical(&mut self, rho: f64, phi: f64, z: f64) -> &mut Self {
        self.x = rho * phi.cos();
        self.y = rho * phi.sin();
        self.z = z;
        self
    }
}

impl coordinates::Spherical for Vector3 {
    #[inline]
    fn from_spherical(radius: f64, phi: f64, theta: f64) -> Self {
        let s = theta.sin();
        Vector3 {
            x: radius * s * phi.cos(),
            y: radius * s * phi.sin(),
            z: radius * theta.cos(),
        }
    }

    #[inline]
    fn set_spherical(&mut self, radius: f64, phi: f64, theta: f64) -> &mut Self {
        let s = theta.sin();
        self.x = radius * s * phi.cos();
        self.y = radius * s * phi.sin();
        self.z = radius * theta.cos();
        self
    }

    #[inline]
    fn unit_radius(phi: f64, theta: f64) -> Self {
        let s = theta.sin();
        Vector3 {
            x: s * phi.cos(),
            y: s * phi.sin(),
            z: theta.cos(),
        }
    }

    #[inline]
    fn unit_theta(phi: f64, theta: f64) -> Self {
        let c = theta.cos();
        Vector3 {
            x: c * phi.cos(),
            y: c * phi.sin(),
            z: -theta.sin(),
        }
    }

    #[inline]
    fn theta(&self) -> f64 {
        (self.x * self.x + self.y * self.y).atan2(self.z)
    }

    #[inline]
    fn set_theta(&mut self, theta: f64) -> &mut Self {
        self.set_spherical(self.magnitude(), self.phi(), theta)
    }
}

impl coordinates::Homogeneous<Vector3> for Vector2 {
    fn from_homogeneous(vector: &Vector3) -> Self {
        Vector2::new(vector.x / vector.z, vector.y / vector.z)
    }

    fn to_homogeneous(&self) -> Vector3 {
        Vector3::new(self.x, self.y, 1.)
    }
}

impl coordinates::Homogeneous<Vector4> for Vector3 {
    fn from_homogeneous(vector: &Vector4) -> Self {
        Vector3::new(vector.x / vector.w, vector.y / vector.w, vector.z / vector.w)
    }

    fn to_homogeneous(&self) -> Vector4 {
        Vector4::new(self.x, self.y, self.z, 1.)
    }
}

impl transforms::Rotation3 for Vector3 {
    fn set_rotation(&mut self, angle: f64, axis: &Vector3) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();
        let k = 1. - c;
        let x = self.x;
        let y = self.y;
        let z = self.z;
        let ux = axis.x;
        let uy = axis.y;
        let uz = axis.z;

        let k_uxy = k * ux * uy;
        let k_uxz = k * ux * uz;
        let k_uyz = k * uy * uz;

        self.x = (k * ux * ux + c) * x + (k_uxy - uz * s) * y + (k_uxz + uy * s) * z;
        self.y = (k_uxy + uz * s) * x + (k * uy * uy + c) * y + (k_uyz - ux * s) * z;
        self.z = (k_uxz - uy * s) * x + (k_uyz + ux * s) * y + (k * uz * uz + c) * z;
        self
    }

    fn set_rotation_x(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();
        let y = self.y;
        let z = self.z;

        self.y = y * c - z * s;
        self.z = y * s + z * c;
        self
    }

    fn set_rotation_y(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();
        let x = self.x;
        let z = self.z;

        self.x = x * c + z * s;
        self.z = z * c - x * s;
        self
    }

    fn set_rotation_z(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();
        let x = self.x;
        let y = self.y;

        self.x = x * c - y * s;
        self.y = x * s + y * c;
        self
    }
}

impl Array<[f64; 2]> for Vector2 {
    #[inline]
    fn to_array(&self) -> [f64; 2] {
        [self.x, self.y]
    }

    #[inline]
    fn set_array(&mut self, array: &[f64; 2]) -> &mut Self {
        self.x = array[0];
        self.y = array[1];
        self
    }
}

impl Array<[f64; 3]> for Vector3 {
    #[inline]
    fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    #[inline]
    fn set_array(&mut self, array: &[f64; 3]) -> &mut Self {
        self.x = array[0];
        self.y = array[1];
        self.z = array[2];
        self
    }
}

impl Array<[f64; 4]> for Vector4 {
    #[inline]
    fn to_array(&self) -> [f64; 4] {
        [self.x, self.y, self.z, self.w]
    }

    #[inline]
    fn set_array(&mut self, array: &[f64; 4]) -> &mut Self {
        self.x = array[0];
        self.y = array[1];
        self.z = array[2];
        self.w = array[3];
        self
    }
}

impl Array<[f64; 6]> for Vector6 {
    #[inline]
    fn to_array(&self) -> [f64; 6] {
        [self.x, self.y, self.z, self.u, self.v, self.w]
    }

    #[inline]
    fn set_array(&mut self, array: &[f64; 6]) -> &mut Self {
        self.x = array[0];
        self.y = array[1];
        self.z = array[2];
        self.u = array[3];
        self.v = array[4];
        self.w = array[5];
        self
    }
}

impl Split<Vector2> for Vector4 {
    fn split(&self) -> [Vector2; 2] {
        [self.upper(), self.lower()]
    }

    fn concat(lhs: &Vector2, rhs: &Vector2) -> Self {
        Vector4::new(lhs.x, lhs.y, rhs.x, rhs.y)
    }

    #[inline]
    fn upper(&self) -> Vector2 {
        Vector2::new(self.x, self.y)
    }

    #[inline]
    fn lower(&self) -> Vector2 {
        Vector2::new(self.z, self.w)
    }

    #[inline]
    fn set_upper(&mut self, vector: &Vector2) -> &mut Self {
        self.x = vector.x;
        self.y = vector.y;
        self
    }

    #[inline]
    fn set_lower(&mut self, vector: &Vector2) -> &mut Self {
        self.z = vector.x;
        self.w = vector.y;
        self
    }
}

impl Split<Vector3> for Vector6 {
    fn split(&self) -> [Vector3; 2] {
        [self.upper(), self.lower()]
    }

    fn concat(lhs: &Vector3, rhs: &Vector3) -> Self {
        Vector6::new(lhs.x, lhs.y, lhs.z, rhs.x, rhs.y, rhs.z)
    }

    #[inline]
    fn upper(&self) -> Vector3 {
        Vector3::new(self.x, self.y, self.z)
    }

    #[inline]
    fn lower(&self) -> Vector3 {
        Vector3::new(self.u, self.v, self.w)
    }

    #[inline]
    fn set_upper(&mut self, vector: &Vector3) -> &mut Self {
        self.x = vector.x;
        self.y = vector.y;
        self.z = vector.z;
        self
    }

    #[inline]
    fn set_lower(&mut self, vector: &Vector3) -> &mut Self {
        self.u = vector.x;
        self.v = vector.y;
        self.w = vector.z;
        self
    }
}

#[cfg(test)]
mod tests {
    mod vector3 {
        use crate::assert_near;
        use crate::common::*;
        use crate::common::transforms::Rotation3;
        use crate::vector;
        use crate::common::coordinates::*;
        use crate::vector::Vector3;

        #[test]
        fn new() {
            let u = Vector3::new(1., 2., 3.);
            assert_eq!(u.x, 1.);
            assert_eq!(u.y, 2.);
            assert_eq!(u.z, 3.);
        }

        #[test]
        fn magnitude() {
            let sqrt_2 = std::f64::consts::SQRT_2;
            let zero = vector::consts::ZEROS_3;
            let u = Vector3::new(1., 1., 0.);
            assert_eq!(u.magnitude2(), 2.);
            assert_eq!(u.magnitude(), sqrt_2);
            assert_eq!(u.distance2(&zero), 2.);
            assert_eq!(u.distance(&zero), sqrt_2);
        }

        #[test]
        fn partial_eq() {
            let u = Vector3::new(-4., 0., 1.);
            let v = Vector3::new(-2., 0., 1.);
            assert_eq!(u, u);
            assert_ne!(u, v);
        }

        #[test]
        fn polar_coordinates() {
            let u = Vector3::ones();
            assert_eq!(u.rho(), std::f64::consts::SQRT_2);
            assert_eq!(u.phi(), std::f64::consts::FRAC_PI_4);
        }

        #[test]
        fn normalized() {
            let mut u = Vector3::new(1., 1., 0.);
            let tol = 10. * std::f64::EPSILON;
            let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
            u.set_normalized();
            assert_near!(u.magnitude2(), 1f64, tol);
            assert_near!(u.x, inv_sqrt2,  tol);
            assert_near!(u.y, inv_sqrt2,  tol);
            assert_near!(u.z, 0.,  tol);
        }

        #[test]
        fn fmt() {
            let u = Vector3::new(1., 2., 3.);
            let formatted = format!("{:?}", u);
            assert_eq!(formatted.as_str(), "( 1.000  2.000  3.000 )");
        }

        #[test]
        fn distance() {
            let u = Vector3::new(1., 1., 0.);
            let v = vector::consts::ZEROS_3;
            assert_eq!(u.distance2(&v), 2f64);
        }

        #[test]
        fn arithmetic() {
            let mut u = Vector3::new(-4., 1., 1.);
            let v = Vector3::new(3., 2., -1.);

            assert_eq!(u + v, Vector3::new(-1., 3., 0.));
            assert_eq!(u - v, Vector3::new(-7., -1., 2.));
            assert_eq!(u * 2., Vector3::new(-8., 2., 2.));
            assert_eq!(u / 4., Vector3::new(-1., 0.25, 0.25));

            u += v;
            assert_eq!(u, Vector3::new(-1., 3., 0.));
        }

        #[test]
        fn angles() {
            let angle = std::f64::consts::FRAC_PI_4;
            let u = vector::consts::EX_3;
            let mut v = u;

            assert_eq!(u.angle(&v), 0.);
            v.set_rotation_z(angle);
            assert_near!(u.angle(&v), angle, std::f64::EPSILON);
            assert_near!(v.angle(&u), angle, std::f64::EPSILON);
            v.set_rotation_z(angle);
            assert_near!(u.angle(&v), 2. * angle, 10. * std::f64::EPSILON);
            v.set_rotation_z(angle);
            assert_near!(u.angle(&v), 3. * angle, 10. * std::f64::EPSILON);
            v.set_rotation_z(angle);
            assert_eq!(u.angle(&v), 4. * angle);
        }

        #[test]
        fn rotations_xyz() {
            let angle = std::f64::consts::FRAC_PI_2;
            let mut u = vector::consts::EX_3;
            let mut v = vector::consts::EY_3;
            let mut w = vector::consts::EZ_3;

            u.set_rotation_z(angle);
            assert_near!(u.distance2(&vector::consts::EY_3), 0., std::f64::EPSILON);

            v.set_rotation_x(angle);
            assert_near!(v.distance2(&vector::consts::EZ_3), 0., std::f64::EPSILON);

            w.set_rotation_y(angle);
            assert_near!(w.distance2(&vector::consts::EX_3), 0., std::f64::EPSILON);
        }

        #[test]
        fn rotations() {
            let angle = std::f64::consts::FRAC_PI_2;
            let mut u = vector::consts::EX_3;
            let mut v = vector::consts::EY_3;
            let mut w = vector::consts::EZ_3;

            let mut axis = vector::consts::EZ_3;
            u.set_rotation(angle, &axis);
            assert_near!(u.distance2(&vector::consts::EY_3), 0., std::f64::EPSILON);

            axis = vector::consts::EX_3;
            v.set_rotation(angle, &axis);
            assert_near!(v.distance2(&vector::consts::EZ_3), 0., std::f64::EPSILON);

            axis = vector::consts::EY_3;
            w.set_rotation(angle, &axis);
            assert_near!(w.distance2(&vector::consts::EX_3), 0., std::f64::EPSILON);
        }
    }
}