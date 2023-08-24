use ndarray::*;

use super::constants;

pub struct Lattice {
    pub u_lbm: f64,
    pub sigma: f64,
    pub rho_lbm: f64,
    pub om_p_lbm: f64,
    pub om_m_lbm: f64,

    pub lx: usize,
    pub ly: usize,
    pub it_max: usize,

    // # Density arrays
    pub g: Array3<f64>,
    pub g_eq: Array3<f64>,
    pub g_up: Array3<f64>,

    // # Boundary conditions
    pub u_left: Array2<f64>,
    pub u_right: Array2<f64>,
    pub u_top: Array2<f64>,
    pub u_bot: Array2<f64>,

    // # Physical fields
    pub rho: Array2<f64>,
    pub u: Array3<f64>,
}

pub fn compute_equilibrium(
    u: &Array3<f64>,
    rho: &Array2<f64>,
    g_eq: &mut Array3<f64>,
) {
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

pub fn compute_macroscopic(
    rho: &mut Array2<f64>,
    g: &Array3<f64>,
    u: &mut Array3<f64>,
) {
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


    // rho.par_map_inplace(|x| *x = 0.0);
    // u.par_map_inplace(|x| *x = 0.0);

    // for i in 0..g.shape()[1] {
    //     rho.index_axis_mut(Axis(0), i)
    //         .assign::<ndarray::Dim<[usize; 1]>, OwnedRepr<f64>>(
    //             &g.slice(s![.., i, ..])
    //                 .axis_iter(Axis(1))
    //                 .into_iter()
    //                 .map(|row| row.sum())
    //                 .collect::<Vec<f64>>()
    //                 .into()
    //         );
    // }

    // for dim in 0..2{
    //     Zip::from(u.index_axis_mut(Axis(0), dim))
    //         .and(&mut *rho) // fresh reborrow to avoid move
    //         .par_for_each(|u_elem, rho_elem| {
    //             *u_elem /= *rho_elem;
    //         });
    // }
}
