pub mod lbmcore;
pub mod render;
pub mod tests;

// sources
// https://github.com/ndbaker1/bloe
// https://github.com/jviquerat/lbm/

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

// pub fn run_main_ndarray() {
//     use std::time::Instant;

//     let start_time = Instant::now();

//     let mut _lattice = get_cavity_lattice();
//     init_lattice(&mut _lattice);

//     for it in 0.._lattice.it_max {
//         if it % 100 == 0 {
//             println!("Iteration: {:>6}/{}", it, _lattice.it_max);
//         }
//         step_lattice(&mut _lattice, it);
//     }

//     let elapsed = start_time.elapsed();
//     println!("Elapsed: {:.2?}", elapsed);
// }
