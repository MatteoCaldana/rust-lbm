use macroquad::prelude::*;
use ndarray::*;

pub mod lbmcore;
pub mod render;

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

pub struct Lattice {
    u_lbm: f64,
    sigma: f64,
    rho_lbm: f64,
    om_p_lbm: f64,
    om_m_lbm: f64,

    lx: usize,
    ly: usize,
    pub it_max: usize,

    // # Density arrays
    g: Array3<f64>,
    g_eq: Array3<f64>,
    g_up: Array3<f64>,

    // # Boundary conditions
    u_left: Array2<f64>,
    u_right: Array2<f64>,
    u_top: Array2<f64>,
    u_bot: Array2<f64>,

    // # Physical fields
    pub rho: Array2<f64>,
    pub u: Array3<f64>,
}

pub fn init_lattice(lattice: &mut Lattice) {
    // Initialize and compute first equilibrium
    set_inlets(
        lattice.sigma,
        lattice.u_lbm,
        0,
        &mut lattice.u_top,
        &mut lattice.u_bot,
        &mut lattice.u_left,
        &mut lattice.u_right,
    );
    lattice.rho *= lattice.rho_lbm;
    lbmcore::lattice::compute_equilibrium(&lattice.u, &lattice.rho, &mut lattice.g_eq);
    lattice.g.assign(&lattice.g_eq);
}

pub fn step_lattice(lattice: &mut Lattice, it: usize) {
    set_inlets(
        lattice.sigma,
        lattice.u_lbm,
        it,
        &mut lattice.u_top,
        &mut lattice.u_bot,
        &mut lattice.u_left,
        &mut lattice.u_right,
    );
    // 2. Compute macroscopic fields
    lbmcore::lattice::compute_macroscopic(&mut lattice.rho, &lattice.g, &mut lattice.u);
    // 4. Compute equilibrium state
    lbmcore::lattice::compute_equilibrium(&lattice.u, &lattice.rho, &mut lattice.g_eq);
    // 5. Streaming
    lbmcore::lattice::collision_and_streaming(
        &mut lattice.g,
        &lattice.g_eq,
        &mut lattice.g_up,
        lattice.om_p_lbm,
        lattice.om_m_lbm,
        lattice.lx,
        lattice.ly,
    );
    // 6. Boundary conditions
    lbmcore::boundary_conditions::zou_he_bottom_wall_velocity(
        &mut lattice.u,
        &lattice.u_bot,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_left_wall_velocity(
        &mut lattice.u,
        &lattice.u_left,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_right_wall_velocity(
        lattice.lx,
        &mut lattice.u,
        &lattice.u_right,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_top_wall_velocity(
        lattice.ly,
        &mut lattice.u,
        &lattice.u_top,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_bottom_left_corner_velocity(
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_top_left_corner_velocity(
        lattice.ly,
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_top_right_corner_velocity(
        lattice.lx,
        lattice.ly,
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
    lbmcore::boundary_conditions::zou_he_bottom_right_corner_velocity(
        lattice.lx,
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
}

pub fn get_cavity_lattice() -> Lattice {
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

    let tau_p_lbm: f64 = tau_lbm;
    let lambda_trt: f64 = 1.0 / 4.0; // Best for stability
    let tau_m_lbm: f64 = lambda_trt / (tau_p_lbm - 0.5) + 0.5;
    let om_p_lbm: f64 = 1.0 / tau_p_lbm;
    let om_m_lbm: f64 = 1.0 / tau_m_lbm;
    //let om_lbm: f64 = 1.0 / tau_lbm;

    return Lattice {
        u_lbm: u_lbm,
        sigma: sigma,
        rho_lbm: rho_lbm,
        om_p_lbm: om_p_lbm,
        om_m_lbm: om_m_lbm,

        lx: lx,
        ly: ly,
        it_max: it_max,

        // # Density arrays
        g: Array3::<f64>::zeros((lbmcore::constants::Q, nx, ny)),
        g_eq: Array3::<f64>::zeros((lbmcore::constants::Q, nx, ny)),
        g_up: Array3::<f64>::zeros((lbmcore::constants::Q, nx, ny)),

        // # Boundary conditions
        u_left: Array2::<f64>::zeros((2, ny)),
        u_right: Array2::<f64>::zeros((2, ny)),
        u_top: Array2::<f64>::zeros((2, nx)),
        u_bot: Array2::<f64>::zeros((2, nx)),

        // # Physical fields
        rho: Array2::<f64>::ones((nx, ny)),
        u: Array3::<f64>::zeros((2, nx, ny)),
    };
}

pub fn run_main_ndarray() {
    use std::time::Instant;

    let start_time = Instant::now();

    let mut _lattice = get_cavity_lattice();
    init_lattice(&mut _lattice);

    for it in 0.._lattice.it_max {
        if it % 100 == 0 {
            println!("Iteration: {:>6}/{}", it, _lattice.it_max);
        }
        step_lattice(&mut _lattice, it);
    }

    let elapsed = start_time.elapsed();
    println!("Elapsed: {:.2?}", elapsed);
}
