use ndarray::*;

use super::constants;

pub struct Lattice {
    pub u_lbm: f64,
    pub sigma: f64,
    pub tau_p_lbm: f64,
    pub tau_m_lbm: f64,
    pub rho_lbm: f64,
    pub om_p_lbm: f64,
    pub om_m_lbm: f64,

    pub nx: usize,
    pub ny: usize,

    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,

    pub dx: f64,
    pub dy: f64,
    pub dt: f64,

    pub lx: usize,
    pub ly: usize,
    pub it_max: usize,

    // shape tag
    pub tag: Array2<u8>,
    pub obs: Vec<[u32; 2]>, 
    pub bnd: Vec<[usize; 3]>, 
    pub ibb: Vec<f64>,

    // # Density arrays
    pub g: Array3<f64>,
    pub g_eq: Array3<f64>,
    pub g_up: Array3<f64>,

    // # Boundary conditions
    pub u_left: Array2<f64>,
    pub u_right: Array2<f64>,
    pub u_top: Array2<f64>,
    pub u_bot: Array2<f64>,
    pub rho_right: Array1<f64>,

    // # Physical fields
    pub rho: Array2<f64>,
    pub u: Array3<f64>,
}

pub fn compute_equilibrium(u: &Array3<f64>, rho: &Array2<f64>, g_eq: &mut Array3<f64>) {
    for q in 0..9 {
        for i in 0..u.shape()[1] {
            for j in 0..u.shape()[2] {
                let v = 1.5 * (u[[0, i, j]] * u[[0, i, j]] + u[[1, i, j]] * u[[1, i, j]]);
                let t = 3.0 * (u[[0, i, j]] * constants::CX[q] + u[[1, i, j]] * constants::CY[q]);
                g_eq[[q, i, j]] = (1.0 + t + 0.5 * t * t - v) * rho[[i, j]] * constants::W[q];
            }
        }
    }
}

pub fn collision_and_streaming(
    g: &mut Array3<f64>,
    g_eq: &Array3<f64>,
    g_up: &mut Array3<f64>,
    om_p: f64,
    om_m: f64,
    lx: usize,
    ly: usize,
) {
    // Take care of q=0 first
    for i in 0..g.shape()[1] {
        for j in 0..g.shape()[2] {
            g_up[[0, i, j]] = (1.0 - om_p) * g[[0, i, j]] + om_p * g_eq[[0, i, j]];
            g[[0, i, j]] = g_up[[0, i, j]];
        }
    }

    // Collide other indices
    for q in 1..9 {
        let qb = constants::NS[q];
        for i in 0..g.shape()[1] {
            for j in 0..g.shape()[2] {
                g_up[[q, i, j]] = (1.0 - 0.5 * (om_p + om_m)) * g[[q, i, j]]
                    - 0.5 * (om_p - om_m) * g[[qb, i, j]]
                    + 0.5 * (om_p + om_m) * g_eq[[q, i, j]]
                    + 0.5 * (om_p - om_m) * g_eq[[qb, i, j]];
            }
        }
    }

    // Stream
    for i in 0..lx {
        for j in 0..ly {
            g[[1, i + 1, j]] = g_up[[1, i, j]];
            g[[2, i, j]] = g_up[[2, i + 1, j]];
            g[[3, i, j + 1]] = g_up[[3, i, j]];
            g[[4, i, j]] = g_up[[4, i, j + 1]];
            g[[5, i + 1, j + 1]] = g_up[[5, i, j]];
            g[[6, i, j]] = g_up[[6, i + 1, j + 1]];
            g[[7, i, j + 1]] = g_up[[7, i + 1, j]];
            g[[8, i + 1, j]] = g_up[[8, i, j + 1]];
        }
    }

    for i in 0..lx {
        g[[1, i + 1, ly]] = g_up[[1, i, ly]];
        g[[2, i, ly]] = g_up[[2, i + 1, ly]];
    }

    for j in 0..ly {
        g[[3, lx, j + 1]] = g_up[[3, lx, j]];
        g[[4, lx, j]] = g_up[[4, lx, j + 1]];
    }
}

pub fn compute_macroscopic(rho: &mut Array2<f64>, g: &Array3<f64>, u: &mut Array3<f64>) {
    // TODO: rewrite with cache alignment
    for i in 0..rho.shape()[0] {
        for j in 0..rho.shape()[1] {
            rho[[i, j]] = 0.0;
            u[[0, i, j]] = 0.0;
            u[[1, i, j]] = 0.0;
            for q in 0..g.shape()[0] {
                rho[[i, j]] += g[[q, i, j]];
                u[[0, i, j]] += constants::CX[q] * g[[q, i, j]];
                u[[1, i, j]] += constants::CY[q] * g[[q, i, j]];
            }
            u[[0, i, j]] /= rho[[i, j]];
            u[[1, i, j]] /= rho[[i, j]];
        }
    }
}

pub fn compute_drag_lift(
    boundary: &Vec<[usize; 3]>,
    g_up: &Array3<f64>,
    g: &Array3<f64>,
    r_ref: f64,
    u_ref: f64,
    l_ref: f64,
) -> (f64, f64) {
    let mut fx = 0.0;
    let mut fy = 0.0;

    for k in 0..boundary.len() {
        let i = boundary[k][0];
        let j = boundary[k][1];
        let q = boundary[k][2];
        let qb = constants::NS[q];
        let g0 = g_up[[q, i, j]] + g[[qb, i, j]];
        fx += g0 * constants::CX[q];
        fy += g0 * constants::CY[q];
    }

    let rescale = -2.0 / (r_ref * l_ref * u_ref * u_ref);
    fx *= rescale;
    fy *= rescale;

    return (fx, fy);
}

pub fn init_lattice(lattice: &mut Lattice, set_inlets: fn(&mut Lattice, usize)) {
    // Initialize and compute first equilibrium
    set_inlets(&mut *lattice, 0);
    lattice.rho *= lattice.rho_lbm;
    compute_equilibrium(&lattice.u, &lattice.rho, &mut lattice.g_eq);
    lattice.g.assign(&lattice.g_eq);
}

pub fn step_lattice(
    lattice: &mut Lattice,
    it: usize,
    set_inlets: fn(&mut Lattice, usize),
    set_bc: fn(&mut Lattice),
    compute_observables: fn(&mut Lattice),
) {
    set_inlets(&mut *lattice, it);
    // 2. Compute macroscopic fields
    compute_macroscopic(&mut lattice.rho, &lattice.g, &mut lattice.u);
    // 4. Compute equilibrium state
    compute_equilibrium(&lattice.u, &lattice.rho, &mut lattice.g_eq);
    // 5. Streaming
    collision_and_streaming(
        &mut lattice.g,
        &lattice.g_eq,
        &mut lattice.g_up,
        lattice.om_p_lbm,
        lattice.om_m_lbm,
        lattice.lx,
        lattice.ly,
    );
    // 6. Boundary conditions
    set_bc(&mut *lattice);
    // 7. compute lift and drag
    compute_observables(&mut *lattice);
}
