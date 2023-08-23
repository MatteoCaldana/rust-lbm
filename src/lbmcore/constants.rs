pub const Q: usize = 9;
pub const CX: &'static [f64] = &[0., 1., -1., 0.,0., 1.,-1.,-1.,1.];
pub const CY: &'static [f64] = &[0.,0.,0.,1.,-1.,1.,-1.,1.,-1.];
pub const W: &'static [f64] = &[
    4. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 9.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
];
pub const NS: &'static [usize] = &[0usize, 2usize, 1usize, 4usize, 3usize, 6usize, 5usize, 8usize, 7usize];