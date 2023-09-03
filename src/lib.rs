pub mod lbm_core;
pub mod lbm_app;
pub mod render;
pub mod tests;
pub mod geometry;

// TODO: parallelize with rayon https://web.dev/webassembly-threads/
// TODO: check correctenss

// fn print_type_of<T>(_: &T) {
//     println!("{}", std::any::type_name::<T>())
// }

// fn meshgrid(x: &Array1<f64>, y: &Array1<f64>) -> (Array2<f64>, Array2<f64>) {
//     let mut xx = Array2::<f64>::ones((x.len(), y.len()));
//     let mut yy = Array2::<f64>::ones((x.len(), y.len()));
//     for i in 0..x.len() {
//         for j in 0..y.len() {
//             xx[[i, j]] = x[i];
//             yy[[i, j]] = y[j];
//         }
//     }
//     return (xx, yy);
// }
