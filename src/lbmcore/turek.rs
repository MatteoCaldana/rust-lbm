use super::boundary_conditions;
use super::constants;
use super::lattice;

use ndarray::*;

pub fn apply_bc(lattice: &mut lattice::Lattice) {
    boundary_conditions::bounce_back_obstacle(
        &lattice.bnd,
        &lattice.ibb,
        &lattice.g_up,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_bottom_wall_velocity(
        &mut lattice.u,
        &lattice.u_bot,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_left_wall_velocity(
        &mut lattice.u,
        &lattice.u_left,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_right_wall_velocity(
        lattice.lx,
        &mut lattice.u,
        &lattice.u_right,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_top_wall_velocity(
        lattice.ly,
        &mut lattice.u,
        &lattice.u_top,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_right_wall_pressure(
        lattice.lx,
        &mut lattice.u,
        &lattice.u_right,
        &mut lattice.rho,
        &lattice.rho_right,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_bottom_left_corner_velocity(
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_top_left_corner_velocity(
        lattice.ly,
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_top_right_corner_velocity(
        lattice.lx,
        lattice.ly,
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
    boundary_conditions::zou_he_bottom_right_corner_velocity(
        lattice.lx,
        &mut lattice.u,
        &mut lattice.rho,
        &mut lattice.g,
    );
}

pub fn get_lattice() -> lattice::Lattice {
    let re_lbm: f64 = 100.0;
    let npts: usize = 100;
    let u_lbm: f64 = 0.05;
    let rho_lbm: f64 = 1.0;
    let t_max: f64 = 0.02;
    let x_min: f64 = -0.2;
    let x_max: f64 = 2.0;
    let y_min: f64 = -0.2;
    let y_max: f64 = 0.21;
    let c_s: f64 = 1.0 / f64::sqrt(3.0);
    let ny: usize = npts;
    let nu_lbm: f64 = u_lbm * (npts as f64) / re_lbm;
    let tau_lbm: f64 = 0.5 + nu_lbm / (c_s * c_s);
    let dt: f64 = re_lbm * nu_lbm / ((npts * npts) as f64);
    let dy: f64 = (y_max - y_min) / ((ny - 1) as f64);
    let nx: usize = f64::round((ny as f64) * (x_max - x_min) / (y_max - y_min)) as usize;
    let dx: f64 = (x_max - x_min) / ((nx - 1) as f64);

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

    return lattice::Lattice {
        u_lbm: u_lbm,
        sigma: sigma,
        rho_lbm: rho_lbm,
        om_p_lbm: om_p_lbm,
        om_m_lbm: om_m_lbm,

        nx: nx,
        ny: ny,

        x_min: x_min,
        x_max: x_max,
        y_min: y_min,
        y_max: y_max,

        dx: dx,
        dy: dy,

        lx: lx,
        ly: ly,
        it_max: it_max,

        tag: Array2::<u8>::zeros((nx, ny)),
        obs: Vec::new(),
        bnd: Vec::new(),
        ibb: Vec::new(),
        // # Density arrays
        g: Array3::<f64>::zeros((constants::Q, nx, ny)),
        g_eq: Array3::<f64>::zeros((constants::Q, nx, ny)),
        g_up: Array3::<f64>::zeros((constants::Q, nx, ny)),

        // # Boundary conditions
        u_left: Array2::<f64>::zeros((2, ny)),
        u_right: Array2::<f64>::zeros((2, ny)),
        u_top: Array2::<f64>::zeros((2, nx)),
        u_bot: Array2::<f64>::zeros((2, nx)),
        rho_right: Array1::<f64>::zeros(ny),

        // # Physical fields
        rho: Array2::<f64>::ones((nx, ny)),
        u: Array3::<f64>::zeros((2, nx, ny)),
    };
}

pub fn set_inlets(lattice: &mut lattice::Lattice, it: usize) {
    let ret: f64 = 1.0 - f64::exp(-((it * it) as f64) / (2.0 * lattice.sigma * lattice.sigma));
    let h = lattice.y_max - lattice.y_min;
    for i in 0..lattice.u_top.shape()[1] {
        lattice.u_top[[0, i]] = 0.0;
        lattice.u_bot[[0, i]] = 0.0;
    }
    for j in 0..lattice.u_left.shape()[1] {
        lattice.u_right[[1, j]] = 0.0;
        lattice.rho_right[j] = lattice.rho_lbm;
        let y = lattice.y_min + (j as f64) * lattice.dy;
        lattice.u_left[[0, j]] =
            ret * lattice.u_lbm * 4.0 * (lattice.y_max - y) * (y - lattice.y_min) / h / h;
        lattice.u_left[[1, j]] = 0.0;
    }
}
