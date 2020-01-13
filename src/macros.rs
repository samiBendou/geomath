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
