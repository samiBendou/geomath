#[macro_export]
macro_rules! assert_near {
    ($val: expr, $exp: expr, $tol: expr) => {
        assert!(
        ($val - $exp).abs() < $tol,
        "Approximation failed\nvalue: {}\nexpected: {}\n tolerance: {}",
        $val, $exp, $tol
        )
    }
}

#[macro_export]
macro_rules! impl_vector {
    ($VectorN:ident { $($field:ident),+ }, $n: expr) => {
        impl From<f64> for $VectorN {
            fn from(s: f64) -> Self {
                $VectorN { $($field: s),+ }
            }
        }

        impl $VectorN {
            #[inline]
            pub fn new($($field: f64),+) -> Self {
                $VectorN { $($field: $field),+ }
            }
        }

        impl Initializer for $VectorN {
            #[inline]
            fn zeros() -> Self {
                $VectorN { $($field: 0.),+ }
            }

            #[inline]
            fn ones() -> Self {
                $VectorN { $($field: 1.),+ }
            }
        }

        impl Reset<$VectorN> for $VectorN {
            #[inline]
            fn reset0(&mut self) -> &mut Self {
                $(self.$field = 0.;)+
                self
            }

            #[inline]
            fn reset1(&mut self) -> &mut Self {
                $(self.$field = 1.;)+
                self
            }

            #[inline]
            fn reset(&mut self, val: &$VectorN) -> &mut Self {
                $(self.$field = val.$field;)+
                self
            }
        }

        impl Metric for $VectorN {

            #[inline]
            fn dot(&self, rhs: &Self) -> f64 {
                let mut ret = 0.;
                $(ret += self.$field * rhs.$field;)+
                ret
            }

            #[inline]
            fn magnitude2(&self) -> f64 {
                let mut ret = 0.;
                $(ret += self.$field * self.$field;)+
                ret
            }

            #[inline]
            fn magnitude(&self) -> f64 {
                self.magnitude2().sqrt()
            }

            #[inline]
            fn distance2(&self, rhs: &Self) -> f64 {
                let mut ret = 0.;
                let mut distance;
                $(
                    distance = self.$field - rhs.$field;
                    ret += distance * distance;
                )+
                ret
            }

            #[inline]
            fn distance(&self, rhs: &Self) -> f64 {
                self.distance2(rhs).sqrt()
            }

            #[inline]
            fn set_normalized(&mut self) -> &mut Self {
                let magnitude = self.magnitude();
                $(self.$field /= magnitude;)+
                self
            }
        }

        impl BitOr<$VectorN> for $VectorN {
            type Output = f64;

            #[inline]
            fn bitor(self, rhs: $VectorN) -> Self::Output {
                self.dot(&rhs)
            }
        }

        impl Not for $VectorN {
            type Output = f64;

            #[inline]
            fn not(self) -> Self::Output {
                self.magnitude()
            }
        }

        impl Rem<$VectorN> for $VectorN {
            type Output = f64;

            #[inline]
            fn rem(self, rhs: Self) -> Self::Output {
                self.distance(&rhs)
            }
        }

        impl PartialEq for $VectorN {

            #[inline]
            fn eq(&self, other: &Self) -> bool {
                self.distance2(other) < std::f64::MIN_POSITIVE
            }

            #[inline]
            fn ne(&self, other: &Self) -> bool {
                self.distance2(other) >= std::f64::MIN_POSITIVE
            }
        }

        impl Add<$VectorN> for $VectorN {
            type Output = Self;

            #[inline]
            fn add(self, rhs: Self) -> Self::Output {
                $VectorN { $($field: self.$field + rhs.$field),+ }
            }
        }

        impl AddAssign<$VectorN> for $VectorN {

            #[inline]
            fn add_assign(&mut self, rhs: $VectorN) {
                $(self.$field += rhs.$field;)+
            }
        }

        impl Sub<$VectorN> for $VectorN {
            type Output = Self;

            #[inline]
            fn sub(self, rhs: Self) -> Self::Output {
                $VectorN { $($field: self.$field - rhs.$field),+ }
            }
        }

        impl SubAssign<$VectorN> for $VectorN {

            #[inline]
            fn sub_assign(&mut self, rhs: $VectorN) {
                $(self.$field -= rhs.$field;)+
            }
        }

        impl Neg for $VectorN {
            type Output = Self;

            #[inline]
            fn neg(self) -> Self::Output {
                $VectorN { $($field: -self.$field),+ }
            }
        }

        impl Mul<f64> for $VectorN {
            type Output = Self;

            #[inline]
            fn mul(self, rhs: f64) -> Self::Output {
                $VectorN { $($field: self.$field * rhs),+ }
            }
        }

        impl MulAssign<f64> for $VectorN {

            #[inline]
            fn mul_assign(&mut self, rhs: f64) {
                $(self.$field *= rhs;)+
            }
        }

        impl Div<f64> for $VectorN {
            type Output = Self;

            #[inline]
            fn div(self, rhs: f64) -> Self::Output {
                $VectorN { $($field: self.$field / rhs),+ }
            }
        }

        impl DivAssign<f64> for $VectorN {

            #[inline]
            fn div_assign(&mut self, rhs: f64) {
                $(self.$field /= rhs;)+
            }
        }

        impl Interpolation for $VectorN {
            fn set_lerp(&mut self, other: &Self, s: f64) -> &mut Self {
                $(self.$field += (other.$field - self.$field) * s;)+
                self
            }

            fn set_herp(&mut self, other: &Self, other1: &Self, other2: &Self, s: f64) -> &mut Self {
                let s2 = s * s;
                let t0 = s2 * (2. * s - 3.) + 1.;
                let t1 = s2 * (s - 2.) + s;
                let t2 = s2 * (s - 1.);
                let t3 = s2 * (3. - 2. * s);
                $(self.$field = self.$field * t0 + other1.$field * t1 + other2.$field * t2 + other.$field * t3;)+
                self
            }

            fn set_berp(&mut self, other: &Self, other1: &Self, other2: &Self, s: f64) -> &mut Self {
                let s2 = s * s;
                let inv = 1. - s;
                let inv2 = inv * inv;
                let t0 = inv2 * inv;
                let t1 = 3. * s * inv2;
                let t2 = 3. * s2 * inv;
                let t3 = s2 * s;
                $(self.$field = self.$field * t0 + other1.$field * t1 + other2.$field * t2 + other.$field * t3;)+
                self
            }
        }
    }
}

#[macro_export]
macro_rules! impl_debug_vector {
($VectorN:ident { $($field:ident),+ }) => {
        impl Debug for $VectorN {
            fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
                let mut buffer = String::from("(");
                if self.magnitude() > 1000. {
                    $(buffer += format!(" {:.3e} ", self.$field).as_str();)+
                } else {
                    $(buffer += format!(" {:.3} ", self.$field).as_str();)+
                }
                buffer += ")";
                write!(f, "{}", buffer)
            }
        }
    }
}

#[macro_export]
macro_rules! impl_debug_matrix {
($MatrixN:ident) => {
        impl Debug for $MatrixN {
            fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
                let rows = self.rows();
                let mut buffer = String::from("");
                for row in rows.iter() {
                    buffer += &format!("\n{:?}", row);
                }
                write!(f, "{}", buffer)
            }
        }
    }
}
