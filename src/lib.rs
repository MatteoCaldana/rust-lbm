use ndarray::*;
pub mod lbmcore;

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

fn set_inlets(
    sigma: f64,
    u_lbm: f64,
    it: usize,
    u_top: &mut Array2<f64>,
    u_bot: &mut Array2<f64>,
    u_left: &mut Array2<f64>,
    u_right: &mut Array2<f64>,
) {
    let ret: f64 = 1.0 - f64::exp(-((it * it) as f64) / (2.0 * sigma * sigma));
    for i in 0..u_top.shape()[1] {
        u_top[[0, i]] = u_lbm * ret;
        u_bot[[0, i]] = 0.0;
        u_left[[1, i]] = 0.0;
        u_right[[1, i]] = 0.0;
    }
}

pub fn run_main_ndarray() {
    use std::time::Instant;
    use std::collections::HashMap;

    let start_time = Instant::now();

    let re_lbm: f64 = 100.0;
    let npts: usize = 100;
    let u_lbm: f64 = 0.2;
    let rho_lbm: f64 = 1.0;
    let t_max: f64 = 20.0;
    let x_min: f64 = 0.0;
    let x_max: f64 = 1.0;
    let y_min: f64 = 0.0;
    let y_max: f64 = 1.0;
    let c_s: f64 = 1.0 / f64::sqrt(3.0);
    let ny: usize = npts;
    let nu_lbm: f64 = u_lbm * (npts as f64) / re_lbm;
    let tau_lbm: f64 = 0.5 + nu_lbm / (c_s * c_s);
    let dt: f64 = re_lbm * nu_lbm / ((npts * npts) as f64);
    //let dx: f64 = (y_max - y_min) / (ny as f64);
    //let dy: f64 = dx;
    let nx: usize = f64::round((ny as f64) * (x_max - x_min) / (y_max - y_min)) as usize;
    let it_max: usize = f64::round(t_max / dt) as usize;
    let sigma: f64 = 10.0 * (nx as f64);

    let lx: usize = nx - 1;
    let ly: usize = ny - 1;
    const Q: usize = 9;

    let tau_p_lbm: f64 = tau_lbm;
    let lambda_trt: f64 = 1.0 / 4.0; // Best for stability
    let tau_m_lbm: f64 = lambda_trt / (tau_p_lbm - 0.5) + 0.5;
    let om_p_lbm: f64 = 1.0 / tau_p_lbm;
    let om_m_lbm: f64 = 1.0 / tau_m_lbm;
    //let om_lbm: f64 = 1.0 / tau_lbm;

    // # D2Q9 Velocities
    let c = arr2(&[
        [0., 0.],
        [1., 0.],
        [-1., 0.],
        [0., 1.],
        [0., -1.],
        [1., 1.],
        [-1., -1.],
        [-1., 1.],
        [1., -1.],
    ]);

    // # Weights
    let w = arr1(&[
        4. / 9.,
        1. / 9.,
        1. / 9.,
        1. / 9.,
        1. / 9.,
        1. / 36.,
        1. / 36.,
        1. / 36.,
        1. / 36.,
    ]);

    // # Array for bounce-back
    let ns = arr1(&[0usize, 2, 1, 4, 3, 6, 5, 8, 7]);

    // # Density arrays
    let mut g = Array3::<f64>::zeros((Q, nx, ny));
    let mut g_eq = Array3::<f64>::zeros((Q, nx, ny));
    let mut g_up = Array3::<f64>::zeros((Q, nx, ny));

    // # Boundary conditions
    let mut u_left = Array2::<f64>::zeros((2, ny));
    let mut u_right = Array2::<f64>::zeros((2, ny));
    let mut u_top = Array2::<f64>::zeros((2, nx));
    let mut u_bot = Array2::<f64>::zeros((2, nx));

    // # Physical fields
    let mut rho = Array2::<f64>::ones((nx, ny));
    let mut u = Array3::<f64>::zeros((2, nx, ny));

    // Initialize and compute first equilibrium
    set_inlets(
        sigma,
        u_lbm,
        0,
        &mut u_top,
        &mut u_bot,
        &mut u_left,
        &mut u_right,
    );
    rho *= rho_lbm;
    lbmcore::lattice::compute_equilibrium(&u, &c, &w, &rho, &mut g_eq);
    g.assign(&g_eq);


    let mut profiler = HashMap::new();
    
    for it in 0..it_max {
        if it % 100 == 0 {
            println!("Iteration: {:>6}/{}", it, it_max);
        }
        // 1. Set inlets
        let t0 = Instant::now();
        set_inlets(
            sigma,
            u_lbm,
            0,
            &mut u_top,
            &mut u_bot,
            &mut u_left,
            &mut u_right,
        );
        profiler.entry("set_inlets").and_modify(|x| *x += t0.elapsed().as_secs_f64()).or_insert(0.0);
        // 2. Compute macroscopic fields
        let t0 = Instant::now();
        lbmcore::lattice::compute_macroscopic(&mut rho, &g, &mut u, &c);
        profiler.entry("compute_macroscopic").and_modify(|x| *x += t0.elapsed().as_secs_f64()).or_insert(0.0);
        // TBD. Output field
        // 4. Compute equilibrium state
        let t0 = Instant::now();
        lbmcore::lattice::compute_equilibrium(&u, &c, &w, &rho, &mut g_eq);
        profiler.entry("compute_equilibrium").and_modify(|x| *x += t0.elapsed().as_secs_f64()).or_insert(0.0);
        // 5. Streaming
        let t0 = Instant::now();
        lbmcore::lattice::collision_and_streaming(
            &mut g, &g_eq, &mut g_up, om_p_lbm, om_m_lbm, &ns, lx, ly,
        );
        profiler.entry("collision_and_streaming").and_modify(|x| *x += t0.elapsed().as_secs_f64()).or_insert(0.0);
        // 6. Boundary conditions
        let t0 = Instant::now();
        lbmcore::boundary_conditions::zou_he_bottom_wall_velocity(&mut u, &u_bot, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_left_wall_velocity(&mut u, &u_left, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_right_wall_velocity(lx, &mut u, &u_right, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_top_wall_velocity(ly, &mut u, &u_top, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_bottom_left_corner_velocity(&mut u, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_top_left_corner_velocity(ly, &mut u, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_top_right_corner_velocity(lx, ly, &mut u, &mut rho, &mut g);
        lbmcore::boundary_conditions::zou_he_bottom_right_corner_velocity(lx, &mut u, &mut rho, &mut g);
        profiler.entry("boundary_conditions").and_modify(|x| *x += t0.elapsed().as_secs_f64()).or_insert(0.0);
        // TBD: Compute observables (drag, lift, etc)
    }

    let elapsed = start_time.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
    for (key, value) in &profiler {
        println!("{:>25}: {}", key, value);
    }
}

// https://github.com/ndbaker1/bloe
