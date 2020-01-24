//!
//! A matrix is represented as a struct containing a `f64` field for each component.
//!

use std::fmt::{Debug, Error, Formatter};
use std::ops::*;
use crate::prelude::*;
use crate::vector::*;
use crate::prelude::transforms;
use crate::{impl_vector, impl_debug_matrix};

/// constant matrices
pub mod consts {
    use super::*;

    /// 2x2 zero matrix
    pub const ZEROS_2: Matrix2 = Matrix2 {xx: 0., xy: 0., yx: 0., yy: 0.};
    /// 2x2 one matrix
    pub const ONES_2: Matrix2 = Matrix2 {xx: 1., xy: 1., yx: 1., yy: 1.};
    /// 2x2 identity matrix
    pub const EYE_2: Matrix2 = Matrix2 {xx: 1., xy: 0., yx: 0., yy: 1.};

    /// 3x3 zero matrix
    pub const ZEROS_3: Matrix3 = Matrix3 {xx: 0., xy: 0., xz: 0., yx: 0., yy: 0., yz: 0., zx: 0., zy: 0., zz: 0.};
    /// 3x3 one matrix
    pub const ONES_3: Matrix3 = Matrix3 {xx: 1., xy: 1., xz: 1., yx: 1., yy: 1., yz: 0., zx: 1., zy: 1., zz: 1.};
    /// 3x3 identity matrix
    pub const EYE_3: Matrix3 = Matrix3 {xx: 1., xy: 0., xz: 0., yx: 0., yy: 1., yz: 0., zx: 0., zy: 0., zz: 1.};

    /// 4x4 zero matrix
    pub const ZEROS_4: Matrix4 = Matrix4 {xx: 0., xy: 0., xz: 0., xw: 0., yx: 0., yy: 0., yz: 0., yw: 0., zx: 0., zy: 0., zz: 0., zw: 0., wx: 0., wy: 0., wz: 0., ww: 0.};
    /// 4x4 one matrix
    pub const ONES_4: Matrix4 = Matrix4 {xx: 1., xy: 1., xz: 1., xw: 1., yx: 1., yy: 1., yz: 1., yw: 1., zx: 1., zy: 1., zz: 1., zw: 1., wx: 1., wy: 1., wz: 1., ww: 1.};
    /// 4x4 identity matrix
    pub const EYE_4: Matrix4 = Matrix4 {xx: 1., xy: 0., xz: 0., xw: 0., yx: 0., yy: 1., yz: 0., yw: 0., zx: 0., zy: 0., zz: 1., zw: 0., wx: 0., wy: 0., wz: 0., ww: 1.};
}

/// 2x2 matrix
#[derive(Copy, Clone)]
pub struct Matrix2 {
    pub xx: f64,
    pub xy: f64,
    pub yx: f64,
    pub yy: f64,
}

/// 3x3 matrix
#[derive(Copy, Clone)]
pub struct Matrix3 {
    pub xx: f64,
    pub xy: f64,
    pub xz: f64,
    pub yx: f64,
    pub yy: f64,
    pub yz: f64,
    pub zx: f64,
    pub zy: f64,
    pub zz: f64,
}

/// 4x4 matrix
#[derive(Copy, Clone)]
pub struct Matrix4 {
    pub xx: f64,
    pub xy: f64,
    pub xz: f64,
    pub xw: f64,
    pub yx: f64,
    pub yy: f64,
    pub yz: f64,
    pub yw: f64,
    pub zx: f64,
    pub zy: f64,
    pub zz: f64,
    pub zw: f64,
    pub wx: f64,
    pub wy: f64,
    pub wz: f64,
    pub ww: f64,
}

impl_vector!(Matrix2 {xx, xy, yx, yy}, 4, mat2);
impl_vector!(Matrix3 {xx, xy, xz, yx, yy, yz, zx, zy, zz}, 9, mat3);
impl_vector!(Matrix4 {xx, xy, xz, xw, yx, yy, yz, yw, zx, zy, zz, zw, wx, wy, wz, ww}, 16, mat4);

impl_debug_matrix!(Matrix2);
impl_debug_matrix!(Matrix3);
impl_debug_matrix!(Matrix4);

impl Rows<[Vector2; 2]> for Matrix2 {
    fn rows(&self) -> [Vector2; 2] {
        [
            vec2(self.xx, self.xy),
            vec2(self.yx, self.yy)
        ]
    }

    fn set_rows(&mut self, rows: &[Vector2; 2]) -> &mut Self {
        self.xx = rows[0].x;
        self.xy = rows[0].y;

        self.yx = rows[1].x;
        self.yy = rows[1].y;
        self
    }
}

impl Rows<[Vector3; 3]> for Matrix3 {
    fn rows(&self) -> [Vector3; 3] {
        [
            vec3(self.xx, self.xy, self.xz),
            vec3(self.yx, self.yy, self.yz),
            vec3(self.zx, self.zy, self.zz),
        ]
    }

    fn set_rows(&mut self, rows: &[Vector3; 3]) -> &mut Self {
        self.xx = rows[0].x;
        self.xy = rows[0].y;
        self.xz = rows[0].z;

        self.yx = rows[1].x;
        self.yy = rows[1].y;
        self.yz = rows[1].z;

        self.zx = rows[2].x;
        self.zy = rows[2].y;
        self.zz = rows[2].z;
        self
    }
}

impl Rows<[Vector4; 4]> for Matrix4 {
    fn rows(&self) -> [Vector4; 4] {
        [
            vec4(self.xx, self.xy, self.xz, self.xw),
            vec4(self.yx, self.yy, self.yz, self.yw),
            vec4(self.zx, self.zy, self.zz, self.zw),
            vec4(self.wx, self.wy, self.wz, self.ww),
        ]
    }

    fn set_rows(&mut self, rows: &[Vector4; 4]) -> &mut Self {
        self.xx = rows[0].x;
        self.xy = rows[0].y;
        self.xz = rows[0].z;
        self.xw = rows[0].w;

        self.yx = rows[1].x;
        self.yy = rows[1].y;
        self.yz = rows[1].z;
        self.yw = rows[1].w;

        self.zx = rows[2].x;
        self.zy = rows[2].y;
        self.zz = rows[2].z;
        self.zw = rows[2].w;

        self.wx = rows[3].x;
        self.wy = rows[3].y;
        self.wz = rows[3].z;
        self.ww = rows[3].w;
        self
    }
}

impl Mul<Matrix2> for Matrix2 {
    type Output = Matrix2;

    fn mul(self, rhs: Matrix2) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl Mul<Vector2> for Matrix2 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl MulAssign<Matrix2> for Matrix2 {
    fn mul_assign(&mut self, rhs: Matrix2) {
        let xx = self.xx;
        let yx = self.yx;
        let xy = self.xy;
        let yy = self.yy;

        self.xx = rhs.xx * xx + rhs.yx * xy;
        self.yx = rhs.xx * yx + rhs.yx * yy;
        self.xy = rhs.xy * xx + rhs.yy * xy;
        self.yy = rhs.xy * yx + rhs.yy * yy;
    }
}

impl Mul<Matrix3> for Matrix3 {
    type Output = Matrix3;

    fn mul(self, rhs: Matrix3) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl Mul<Vector3> for Matrix3 {
    type Output = Vector3;

    fn mul(self, rhs: Vector3) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl Mul<Vector2> for Matrix3 {
    type Output = Vector2;

    fn mul(self, rhs: Vector2) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl MulAssign<Matrix3> for Matrix3 {
    fn mul_assign(&mut self, rhs: Matrix3) {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;

        self.xx = rhs.xx * xx + rhs.yx * xy + rhs.zx * xz;
        self.yx = rhs.xx * yx + rhs.yx * yy + rhs.zx * yz;
        self.zx = rhs.xx * zx + rhs.yx * zy + rhs.zx * zz;

        self.xy = rhs.xy * xx + rhs.yy * xy + rhs.zy * xz;
        self.yy = rhs.xy * yx + rhs.yy * yy + rhs.zy * yz;
        self.zy = rhs.xy * zx + rhs.yy * zy + rhs.zy * zz;

        self.xz = rhs.xz * xx + rhs.yz * xy + rhs.zz * xz;
        self.yz = rhs.xz * yx + rhs.yz * yy + rhs.zz * yz;
        self.zz = rhs.xz * zx + rhs.yz * zy + rhs.zz * zz;
    }
}

impl Mul<Matrix4> for Matrix4 {
    type Output = Matrix4;

    fn mul(self, rhs: Matrix4) -> Self::Output {
        let mut ret = self;
        ret *= rhs;
        ret
    }
}

impl Mul<Vector4> for Matrix4 {
    type Output = Vector4;

    fn mul(self, rhs: Vector4) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl Mul<Vector3> for Matrix4 {
    type Output = Vector3;

    fn mul(self, rhs: Vector3) -> Self::Output {
        let mut ret = rhs;
        ret *= self;
        ret
    }
}

impl MulAssign<Matrix4> for Matrix4 {
    fn mul_assign(&mut self, rhs: Matrix4) {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let wx = self.wx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let wy = self.wy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;
        let wz = self.wz;
        let xw = self.xw;
        let yw = self.yw;
        let zw = self.zw;
        let ww = self.ww;

        self.xx = rhs.xx * xx + rhs.yx * xy + rhs.zx * xz + rhs.wx * xw;
        self.yx = rhs.xx * yx + rhs.yx * yy + rhs.zx * yz + rhs.wx * yw;
        self.zx = rhs.xx * zx + rhs.yx * zy + rhs.zx * zz + rhs.wx * zw;
        self.wx = rhs.xx * wx + rhs.yx * wy + rhs.zx * wz + rhs.wx * ww;

        self.xy = rhs.xy * xx + rhs.yy * xy + rhs.zy * xz + rhs.wy * xw;
        self.yy = rhs.xy * yx + rhs.yy * yy + rhs.zy * yz + rhs.wy * yw;
        self.zy = rhs.xy * zx + rhs.yy * zy + rhs.zy * zz + rhs.wy * zw;
        self.wy = rhs.xy * wx + rhs.yy * wy + rhs.zy * wz + rhs.wy * ww;

        self.xz = rhs.xz * xx + rhs.yz * xy + rhs.zz * xz + rhs.wz * xw;
        self.yz = rhs.xz * yx + rhs.yz * yy + rhs.zz * yz + rhs.wz * yw;
        self.zz = rhs.xz * zx + rhs.yz * zy + rhs.zz * zz + rhs.wz * zw;
        self.wz = rhs.xz * wx + rhs.yz * wy + rhs.zz * wz + rhs.wz * ww;

        self.xw = rhs.xw * xx + rhs.yw * xy + rhs.zw * xz + rhs.ww * xw;
        self.yw = rhs.xw * yx + rhs.yw * yy + rhs.zw * yz + rhs.ww * yw;
        self.zw = rhs.xw * zx + rhs.yw * zy + rhs.zw * zz + rhs.ww * zw;
        self.ww = rhs.xw * wx + rhs.yw * wy + rhs.zw * wz + rhs.ww * ww;
    }
}

impl Algebra<Matrix2> for Matrix2 {
    fn determinant(&self) -> f64 {
        let xx = self.xx;
        let yx = self.yx;
        let xy = self.xy;
        let yy = self.yy;

        xx * yy - xy * yx
    }

    fn set_inverse(&mut self) -> &mut Self {
        let xx = self.xx;
        let yx = self.yx;
        let xy = self.xy;
        let yy = self.yy;

        let mut det = xx * yy - xy * yx;

        if det == 0. {
            return self;
        }

        det = 1. / det;

        self.xx = yy * det;
        self.xy = -xy * det;
        self.yx = -yx * det;
        self.yy = xx * det;

        self
    }

    fn set_transposed(&mut self) -> &mut Self {
        let yx = self.yx;

        self.yx = self.xy;
        self.xy = yx;

        self
    }

    fn set_adjugate(&mut self) -> &mut Self {
        let xx = self.xx;

        self.xx = self.yy;
        self.yx = -self.yx;
        self.xy = -self.xy;
        self.yy = xx;

        self
    }
}

impl Algebra<Matrix3> for Matrix3 {
    fn determinant(&self) -> f64 {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;

        let dyx = zz * yy - zy * yz;
        let dyy = -zz * xy + zy * xz;
        let dyz = yz * xy - yy * xz;

        xx * dyx + yx * dyy + zx * dyz
    }

    fn set_inverse(&mut self) -> &mut Self {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;

        let dyx = zz * yy - zy * yz;
        let dyy = -zz * xy + zy * xz;
        let dyz = yz * xy - yy * xz;
        let mut det = xx * dyx + yx * dyy + zx * dyz;

        if det == 0. {
            return self;
        }

        det = 1. / det;
        self.xx = dyx * det;
        self.yx = (-zz * yx + zx * yz) * det;
        self.zx = (zy * yx - zx * yy) * det;
        self.xy = dyy * det;
        self.yy = (zz * xx - zx * xz) * det;
        self.zy = (-zy * xx + zx * xy) * det;
        self.zx = dyz * det;
        self.zy = (-yz * xx + yx * xz) * det;
        self.zz = (yy * xx - yx * xy) * det;

        self
    }

    fn set_transposed(&mut self) -> &mut Self {
        let yx = self.yx;
        let zx = self.zx;
        let yz = self.yz;

        self.yx = self.xy;
        self.xy = yx;
        self.zx = self.xz;
        self.xz = zx;
        self.yz = self.zy;
        self.zy = yz;

        self
    }

    fn set_adjugate(&mut self) -> &mut Self {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;

        self.xx = yy * zz - zy * yz;
        self.yx = zx * yz - yx * zz;
        self.zx = yx * zy - zx * yy;
        self.xy = zy * xz - xy * zz;
        self.yy = xx * zz - zx * xz;
        self.zy = zx * xy - xx * zy;
        self.xz = xy * yz - yy * xz;
        self.yz = yx * xz - xx * yz;
        self.zz = xx * yy - yx * xy;

        self
    }
}

impl Algebra<Matrix4> for Matrix4 {
    fn determinant(&self) -> f64 {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let wx = self.wx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let wy = self.wy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;
        let wz = self.wz;
        let xw = self.xw;
        let yw = self.yw;
        let zw = self.zw;
        let ww = self.ww;

        let d00 = xx * yy - yx * xy;
        let d01 = xx * zy - zx * xy;
        let d02 = xx * wy - wx * xy;
        let d03 = yx * zy - zx * yy;
        let d04 = yx * wy - wx * yy;
        let d05 = zx * wy - wx * zy;
        let d06 = xz * yw - yz * xw;
        let d07 = xz * zw - zz * xw;
        let d08 = xz * ww - wz * xw;
        let d09 = yz * zw - zz * yw;
        let d10 = yz * ww - wz * yw;
        let d11 = zz * ww - wz * zw;

        d00 * d11 - d01 * d10 + d02 * d09 + d03 * d08 - d04 * d07 + d05 * d06
    }

    fn set_inverse(&mut self) -> &mut Self {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let wx = self.wx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let wy = self.wy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;
        let wz = self.wz;
        let xw = self.xw;
        let yw = self.yw;
        let zw = self.zw;
        let ww = self.ww;

        let d00 = xx * yy - yx * xy;
        let d01 = xx * zy - zx * xy;
        let d02 = xx * wy - wx * xy;
        let d03 = yx * zy - zx * yy;
        let d04 = yx * wy - wx * yy;
        let d05 = zx * wy - wx * zy;
        let d06 = xz * yw - yz * xw;
        let d07 = xz * zw - zz * xw;
        let d08 = xz * ww - wz * xw;
        let d09 = yz * zw - zz * yw;
        let d10 = yz * ww - wz * yw;
        let d11 = zz * ww - wz * zw;

        let mut det = d00 * d11 - d01 * d10 + d02 * d09 + d03 * d08 - d04 * d07 + d05 * d06;

        if det == 0. {
            return self;
        }

        det = 1.0 / det;

        self.xx = (yy * d11 - zy * d10 + wy * d09) * det;
        self.yx = (zx * d10 - yx * d11 - wx * d09) * det;
        self.zx = (yw * d05 - zw * d04 + ww * d03) * det;
        self.wx = (zz * d04 - yz * d05 - wz * d03) * det;
        self.xy = (zy * d08 - xy * d11 - wy * d07) * det;
        self.yy = (xx * d11 - zx * d08 + wx * d07) * det;
        self.zy = (zw * d02 - xw * d05 - ww * d01) * det;
        self.wy = (xz * d05 - zz * d02 + wz * d01) * det;
        self.xz = (xy * d10 - yy * d08 + wy * d06) * det;
        self.yz = (yx * d08 - xx * d10 - wx * d06) * det;
        self.zz = (xw * d04 - yw * d02 + ww * d00) * det;
        self.wz = (yz * d02 - xz * d04 - wz * d00) * det;
        self.xw = (yy * d07 - xy * d09 - zy * d06) * det;
        self.yw = (xx * d09 - yx * d07 + zx * d06) * det;
        self.zw = (yw * d01 - xw * d03 - zw * d00) * det;
        self.ww = (xz * d03 - yz * d01 + zz * d00) * det;

        self
    }

    fn set_transposed(&mut self) -> &mut Self {
        let yx = self.yx;
        let zx = self.zx;
        let wx = self.wx;
        let zy = self.zy;
        let wy = self.wy;
        let wz = self.wz;

        self.yx = self.xy;
        self.zx = self.xz;
        self.wx = self.xw;
        self.xy = yx;
        self.zy = self.yz;
        self.wy = self.yw;
        self.xz = zx;
        self.yz = zy;
        self.wz = self.zw;
        self.xw = wx;
        self.yw = wy;
        self.zw = wz;

        self
    }

    fn set_adjugate(&mut self) -> &mut Self {
        let xx = self.xx;
        let yx = self.yx;
        let zx = self.zx;
        let wx = self.wx;
        let xy = self.xy;
        let yy = self.yy;
        let zy = self.zy;
        let wy = self.wy;
        let xz = self.xz;
        let yz = self.yz;
        let zz = self.zz;
        let wz = self.wz;
        let xw = self.xw;
        let yw = self.yw;
        let zw = self.zw;
        let ww = self.ww;

        self.xx = yy * (zz * ww - wz * zw) - yz * (zy * ww - wy * zw) + yw * (zy * wz - wy * zz);
        self.yx = -(yx * (zz * ww - wz * zw) - yz * (zx * ww - wx * zw) + yw * (zx * wz - wx * zz));
        self.zx = yx * (zy * ww - wy * zw) - yy * (zx * ww - wx * zw) + yw * (zx * wy - wx * zy);
        self.wx = -(yx * (zy * wz - wy * zz) - yy * (zx * wz - wx * zz) + yz * (zx * wy - wx * zy));
        self.xy = -(xy * (zz * ww - wz * zw) - xz * (zy * ww - wy * zw) + xw * (zy * wz - wy * zz));
        self.yy = xx * (zz * ww - wz * zw) - xz * (zx * ww - wx * zw) + xw * (zx * wz - wx * zz);
        self.zy = -(xx * (zy * ww - wy * zw) - xy * (zx * ww - wx * zw) + xw * (zx * wy - wx * zy));
        self.wy = xx * (zy * wz - wy * zz) - xy * (zx * wz - wx * zz) + xz * (zx * wy - wx * zy);
        self.xz = xy * (yz * ww - wz * yw) - xz * (yy * ww - wy * yw) + xw * (yy * wz - wy * yz);
        self.yz = -(xx * (yz * ww - wz * yw) - xz * (yx * ww - wx * yw) + xw * (yx * wz - wx * yz));
        self.zz = xx * (yy * ww - wy * yw) - xy * (yx * ww - wx * yw) + xw * (yx * wy - wx * yy);
        self.wz = -(xx * (yy * wz - wy * yz) - xy * (yx * wz - wx * yz) + xz * (yx * wy - wx * yy));
        self.xw = -(xy * (yz * zw - zz * yw) - xz * (yy * zw - zy * yw) + xw * (yy * zz - zy * yz));
        self.yw = xx * (yz * zw - zz * yw) - xz * (yx * zw - zx * yw) + xw * (yx * zz - zx * yz);
        self.zw = -(xx * (yy * zw - zy * yw) - xy * (yx * zw - zx * yw) + xw * (yx * zy - zx * yy));
        self.ww = xx * (yy * zz - zy * yz) - xy * (yx * zz - zx * yz) + xz * (yx * zy - zx * yy);

        self
    }
}

impl transforms::Translation<Vector2> for Matrix3 {
    fn set_translation(&mut self, vector: &Vector2) -> &mut Self {
        self.xx = 1.;
        self.xy = 0.;
        self.xz = vector.x;
        self.yx = 0.;
        self.yy = 1.;
        self.yz = vector.y;
        self.zx = 0.;
        self.zy = 0.;
        self.zz = 1.;

        self
    }
}

impl transforms::Translation<Vector3> for Matrix4 {
    fn set_translation(&mut self, vector: &Vector3) -> &mut Self {
        self.xx = 1.;
        self.xy = 0.;
        self.xz = 0.;
        self.xw = vector.x;
        self.yx = 0.;
        self.yy = 1.;
        self.yz = 0.;
        self.yw = vector.y;
        self.zx = 0.;
        self.zy = 0.;
        self.zz = 1.;
        self.zw = vector.z;
        self.wx = 0.;
        self.wy = 0.;
        self.wz = 0.;
        self.ww = 1.;

        self
    }
}

impl transforms::Rigid<Matrix2, Vector2> for Matrix3 {
    fn set_rigid(&mut self, rotation: &Matrix2, vector: &Vector2) -> &mut Self {
        self.xx = rotation.xx;
        self.xy = rotation.xy;
        self.xz = vector.x;
        self.yx = rotation.yx;
        self.yy = rotation.yy;
        self.yz = vector.y;
        self.zx = 0.;
        self.zy = 0.;
        self.zz = 1.;

        self
    }
}

impl transforms::Rigid<Matrix3, Vector3> for Matrix4 {
    fn set_rigid(&mut self, rotation: &Matrix3, vector: &Vector3) -> &mut Self {
        self.xx = rotation.xx;
        self.xy = rotation.xy;
        self.xz = rotation.xz;
        self.xw = vector.x;
        self.yx = rotation.yx;
        self.yy = rotation.yy;
        self.yz = rotation.yz;
        self.yw = vector.y;
        self.zx = rotation.zx;
        self.zy = rotation.zy;
        self.zz = rotation.zz;
        self.zw = vector.z;
        self.wx = 0.;
        self.wy = 0.;
        self.wz = 0.;
        self.ww = 1.;

        self
    }
}

impl transforms::Similarity<Matrix2, Vector2> for Matrix3 {
    fn set_similarity(&mut self, scale: f64, rotation: &Matrix2, vector: &Vector2) -> &mut Self {
        self.xx = scale * rotation.xx;
        self.xy = scale * rotation.xy;
        self.xz = vector.x;
        self.yx = scale * rotation.yx;
        self.yy = scale * rotation.yy;
        self.yz = vector.y;
        self.zx = 0.;
        self.zy = 0.;
        self.zz = 1.;

        self
    }
}

impl transforms::Similarity<Matrix3, Vector3> for Matrix4 {
    fn set_similarity(&mut self, scale: f64, rotation: &Matrix3, vector: &Vector3) -> &mut Self {
        self.xx = scale * rotation.xx;
        self.xy = scale * rotation.xy;
        self.xz = scale * rotation.xz;
        self.xw = vector.x;
        self.yx = scale * rotation.yx;
        self.yy = scale * rotation.yy;
        self.yz = scale * rotation.yz;
        self.yw = vector.y;
        self.zx = scale * rotation.zx;
        self.zy = scale * rotation.zy;
        self.zz = scale * rotation.zz;
        self.zw = vector.z;
        self.wx = 0.;
        self.wy = 0.;
        self.wz = 0.;
        self.ww = 1.;

        self
    }
}

impl transforms::Rotation2 for Matrix2 {
    fn set_rotation(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();

        self.xx = c;
        self.xy = -s;
        self.yx = s;
        self.yy = c;

        self
    }
}

impl transforms::Rotation3 for Matrix3 {
    fn set_rotation(&mut self, angle: f64, axis: &Vector3) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();
        let k = 1. - c;

        let ux = axis.x;
        let uy = axis.y;
        let uz = axis.z;

        let k_uxy = k * ux * uy;
        let k_uxz = k * ux * uz;
        let k_uyz = k * uy * uz;

        self.xx = k * ux * ux + c;
        self.xy = k_uxy - uz * s;
        self.xz = k_uxz + uy * s;
        self.yx = k_uxy + uz * s;
        self.yy = k * uy * uy + c;
        self.yz = k_uyz - ux * s;
        self.zx = k_uxz - uy * s;
        self.zy = k_uyz + ux * s;
        self.zz = k * uz * uz + c;

        self
    }

    fn set_rotation_x(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();

        self.xx = 1.;
        self.xy = 0.;
        self.xz = 0.;
        self.yx = 0.;
        self.yy = c;
        self.yz = -s;
        self.zx = 0.;
        self.zy = s;
        self.zz = c;

        self
    }

    fn set_rotation_y(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();

        self.xx = c;
        self.xy = 0.;
        self.xz = s;
        self.yx = 0.;
        self.yy = 1.;
        self.yz = 0.;
        self.zx = -s;
        self.zy = 0.;
        self.zz = c;

        self
    }

    fn set_rotation_z(&mut self, angle: f64) -> &mut Self {
        let c = angle.cos();
        let s = angle.sin();

        self.xx = c;
        self.xy = -s;
        self.xz = 0.;
        self.yx = s;
        self.yy = c;
        self.yz = 0.;
        self.zx = 0.;
        self.zy = 0.;
        self.zz = 1.;

        self
    }
}

#[cfg(test)]
mod tests {
    mod matrix3 {
        use crate::prelude::transforms::Rotation3;
        use crate::vector::{self};
        use crate::matrix::{self, *};

        #[test]
        fn arithmetic() {
            let a = matrix::consts::EYE_3 * 2.;
            let b = matrix::consts::ONES_3;
            let c = a * b;
            assert_eq!(c, b * 2.);
        }

        #[test]
        fn determinant() {
            let a = matrix::consts::EYE_3 * 2.;
            assert_eq!(a.determinant(), 8.);
        }

        #[test]
        fn transposed() {
            let a = mat3(1., 2., 3., 4., 5., 6., 7., 8., 9.);
            assert_eq!(a.transposed(), mat3(1., 4., 7., 2., 5., 8., 3., 6., 9.));
        }

        #[test]
        fn inverse() {
            let a = matrix::consts::EYE_3 * 2.;
            assert_eq!(a.inverse(), matrix::consts::EYE_3 * 0.5);
        }

        #[test]
        fn adjugate() {
            let a = matrix::consts::EYE_3 * 2.;
            assert_eq!(a.adjugate(), a.inverse() * a.determinant());
        }

        #[test]
        fn rotations() {
            let angle = std::f64::consts::FRAC_PI_8;
            let rot_x = Matrix3::from_rotation_x(angle);
            let rot_y = Matrix3::from_rotation_y(angle);
            let rot_z = Matrix3::from_rotation_z(angle);

            let mut axis = vector::consts::EX_3;
            let mut rot = Matrix3::from_rotation(angle, &axis);
            assert_eq!(rot, rot_x);

            axis = vector::consts::EY_3;
            rot = Matrix3::from_rotation(angle, &axis);
            assert_eq!(rot, rot_y);

            axis = vector::consts::EZ_3;
            rot = Matrix3::from_rotation(angle, &axis);
            assert_eq!(rot, rot_z);
        }
    }

    mod matrix4 {
        use crate::assert_near;
        use crate::prelude::*;
        use crate::prelude::coordinates::Homogeneous;
        use crate::prelude::transforms::{Rigid, Rotation3, Similarity, Translation};
        use crate::matrix::{self, *};
        use crate::vector::{self, *};

        #[test]
        fn arithmetic() {
            let a = matrix::consts::EYE_4 * 2.;
            let b = matrix::consts::ONES_4;
            let c = a * b;
            assert_eq!(c, b * 2.);
        }

        #[test]
        fn determinant() {
            let a = matrix::consts::EYE_4 * 2.;
            assert_eq!(a.determinant(), 16.);
        }

        #[test]
        fn transposed() {
            let a = mat4(1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.);
            assert_eq!(a.transposed(), mat4(1., 5., 9., 13., 2., 6., 10., 14., 3., 7., 11., 15., 4., 8., 12., 16.));
        }

        #[test]
        fn inverse() {
            let a = matrix::consts::EYE_4 * 2.;
            assert_eq!(a.inverse(), matrix::consts::EYE_4 * 0.5);
        }

        #[test]
        fn adjugate() {
            let a = matrix::consts::EYE_4 * 2.;
            assert_eq!(a.adjugate(), a.inverse() * a.determinant());
        }

        #[test]
        fn translations() {
            let unit_x = vector::consts::EX_3;
            let a = Matrix4::from_translation(&unit_x);
            let u = unit_x.to_homogeneous();
            let translated = a * u;
            assert_eq!(Vector3::from_homogeneous(&translated), (unit_x * 2.));
        }

        #[test]
        fn rigid() {
            let angle = std::f64::consts::FRAC_PI_2;
            let unit_x = vector::consts::EX_3;
            let rotation_z = Matrix3::from_rotation_z(angle);
            let a = Matrix4::from_rigid(&rotation_z, &unit_x);
            let u = unit_x.to_homogeneous();
            let moved = a * u;
            assert_eq!(Vector3::from_homogeneous(&moved), vec3(1., 1., 0.));
        }

        #[test]
        fn similarity() {
            let angle = std::f64::consts::FRAC_PI_2;
            let scale = 2.;
            let unit_x = vector::consts::EX_3;
            let rotation_z = Matrix3::from_rotation_z(angle);
            let a = Matrix4::from_similarity(scale, &rotation_z, &unit_x);
            let u = unit_x.to_homogeneous();
            let moved = a * u;
            assert_near!(Vector3::from_homogeneous(&moved).distance2(&vec3(1., 2., 0.)), 0.,  std::f64::EPSILON);
        }
    }
}