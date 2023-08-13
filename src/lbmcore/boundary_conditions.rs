use ndarray::*;

pub fn zou_he_bottom_left_corner_velocity(
    u: &mut Array3<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    u[[0, 0, 0]] = u[[0, 1, 0]];
    u[[1, 0, 0]] = u[[1, 1, 0]];

    rho[[0, 0]] = rho[[1, 0]];

    g[[1, 0, 0]] = g[[2, 0, 0]] + (2.0 / 3.0) * rho[[0, 0]] * u[[0, 0, 0]];

    g[[3, 0, 0]] = g[[4, 0, 0]] + (2.0 / 3.0) * rho[[0, 0]] * u[[1, 0, 0]];

    g[[5, 0, 0]] = g[[6, 0, 0]]
        + (1.0 / 6.0) * rho[[0, 0]] * u[[0, 0, 0]]
        + (1.0 / 6.0) * rho[[0, 0]] * u[[1, 0, 0]];

    g[[7, 0, 0]] = 0.0;
    g[[8, 0, 0]] = 0.0;

    g[[0, 0, 0]] = rho[[0, 0]]
        - g[[1, 0, 0]]
        - g[[2, 0, 0]]
        - g[[3, 0, 0]]
        - g[[4, 0, 0]]
        - g[[5, 0, 0]]
        - g[[6, 0, 0]]
        - g[[7, 0, 0]]
        - g[[8, 0, 0]];
}

pub fn zou_he_top_left_corner_velocity(
    ly: usize,
    u: &mut Array3<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    u[[0, 0, ly]] = u[[0, 1, ly]];
    u[[1, 0, ly]] = u[[1, 1, ly]];

    rho[[0, ly]] = rho[[1, ly]];

    g[[1, 0, ly]] = g[[2, 0, ly]] + (2.0 / 3.0) * rho[[0, ly]] * u[[0, 0, ly]];

    g[[4, 0, ly]] = g[[3, 0, ly]] - (2.0 / 3.0) * rho[[0, ly]] * u[[1, 0, ly]];

    g[[8, 0, ly]] = g[[7, 0, ly]] + (1.0 / 6.0) * rho[[0, ly]] * u[[0, 0, ly]]
        - (1.0 / 6.0) * rho[[0, ly]] * u[[1, 0, ly]];

    g[[5, 0, ly]] = 0.0;
    g[[6, 0, ly]] = 0.0;

    g[[0, 0, ly]] = rho[[0, ly]]
        - g[[1, 0, ly]]
        - g[[2, 0, ly]]
        - g[[3, 0, ly]]
        - g[[4, 0, ly]]
        - g[[5, 0, ly]]
        - g[[6, 0, ly]]
        - g[[7, 0, ly]]
        - g[[8, 0, ly]];
}

pub fn zou_he_top_right_corner_velocity(
    lx: usize,
    ly: usize,
    u: &mut Array3<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    u[[0, lx, ly]] = u[[0, lx - 1, ly]];
    u[[1, lx, ly]] = u[[1, lx - 1, ly]];

    rho[[lx, ly]] = rho[[lx - 1, ly]];

    g[[2, lx, ly]] = g[[1, lx, ly]] - (2.0 / 3.0) * rho[[lx, ly]] * u[[0, lx, ly]];

    g[[4, lx, ly]] = g[[3, lx, ly]] - (2.0 / 3.0) * rho[[lx, ly]] * u[[1, lx, ly]];

    g[[6, lx, ly]] = g[[5, lx, ly]]
        - (1.0 / 6.0) * rho[[lx, ly]] * u[[0, lx, ly]]
        - (1.0 / 6.0) * rho[[lx, ly]] * u[[1, lx, ly]];

    g[[7, lx, ly]] = 0.0;
    g[[8, lx, ly]] = 0.0;

    g[[0, lx, ly]] = rho[[lx, ly]]
        - g[[1, lx, ly]]
        - g[[2, lx, ly]]
        - g[[3, lx, ly]]
        - g[[4, lx, ly]]
        - g[[5, lx, ly]]
        - g[[6, lx, ly]]
        - g[[7, lx, ly]]
        - g[[8, lx, ly]];
}

pub fn zou_he_bottom_right_corner_velocity(
    lx: usize,
    u: &mut Array3<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    u[[0, lx, 0]] = u[[0, lx - 1, 0]];
    u[[1, lx, 0]] = u[[1, lx - 1, 0]];

    rho[[lx, 0]] = rho[[lx - 1, 0]];

    g[[2, lx, 0]] = g[[1, lx, 0]] - (2.0 / 3.0) * rho[[lx, 0]] * u[[0, lx, 0]];

    g[[3, lx, 0]] = g[[4, lx, 0]] + (2.0 / 3.0) * rho[[lx, 0]] * u[[1, lx, 0]];

    g[[7, lx, 0]] = g[[8, lx, 0]] - (1.0 / 6.0) * rho[[lx, 0]] * u[[0, lx, 0]]
        + (1.0 / 6.0) * rho[[lx, 0]] * u[[1, lx, 0]];

    g[[5, lx, 0]] = 0.0;
    g[[6, lx, 0]] = 0.0;

    g[[0, lx, 0]] = rho[[lx, 0]]
        - g[[1, lx, 0]]
        - g[[2, lx, 0]]
        - g[[3, lx, 0]]
        - g[[4, lx, 0]]
        - g[[5, lx, 0]]
        - g[[6, lx, 0]]
        - g[[7, lx, 0]]
        - g[[8, lx, 0]];
}

pub fn zou_he_left_wall_velocity(
    u: &mut Array3<f64>,
    u_left: &Array2<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    for i in 0..u.shape()[2] {
        u[[0, 0, i]] = u_left[[0, i]];
        u[[1, 0, i]] = u_left[[1, i]];

        rho[[0, i]] = (g[[0, 0, i]]
            + g[[3, 0, i]]
            + g[[4, 0, i]]
            + 2.0 * g[[2, 0, i]]
            + 2.0 * g[[6, 0, i]]
            + 2.0 * g[[7, 0, i]])
            / (1.0 - u[[0, 0, i]]);

        g[[1, 0, i]] = g[[2, 0, i]] + (2.0 / 3.0) * rho[[0, i]] * u[[0, 0, i]];

        g[[5, 0, i]] = g[[6, 0, i]] - (1.0 / 2.0) * (g[[3, 0, i]] - g[[4, 0, i]])
            + (1.0 / 6.0) * rho[[0, i]] * u[[0, 0, i]]
            + (1.0 / 2.0) * rho[[0, i]] * u[[1, 0, i]];

        g[[8, 0, i]] = g[[7, 0, i]]
            + (1.0 / 2.0) * (g[[3, 0, i]] - g[[4, 0, i]])
            + (1.0 / 6.0) * rho[[0, i]] * u[[0, 0, i]]
            - (1.0 / 2.0) * rho[[0, i]] * u[[1, 0, i]];
    }
}

pub fn zou_he_right_wall_velocity(
    lx: usize,
    u: &mut Array3<f64>,
    u_right: &Array2<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    for i in 0..u.shape()[2] {
        u[[0, lx, i]] = u_right[[0, i]];
        u[[1, lx, i]] = u_right[[1, i]];

        rho[[lx, i]] = (g[[0, lx, i]]
            + g[[3, lx, i]]
            + g[[4, lx, i]]
            + 2.0 * g[[1, lx, i]]
            + 2.0 * g[[5, lx, i]]
            + 2.0 * g[[8, lx, i]])
            / (1.0 + u[[0, lx, i]]);

        g[[2, lx, i]] = g[[1, lx, i]] - (2.0 / 3.0) * rho[[lx, i]] * u[[0, lx, i]];

        g[[6, lx, i]] = g[[5, lx, i]] + (1.0 / 2.0) * (g[[3, lx, i]] - g[[4, lx, i]])
            - (1.0 / 6.0) * rho[[lx, i]] * u[[0, lx, i]]
            - (1.0 / 2.0) * rho[[lx, i]] * u[[1, lx, i]];

        g[[7, lx, i]] = g[[8, lx, i]]
            - (1.0 / 2.0) * (g[[3, lx, i]] - g[[4, lx, i]])
            - (1.0 / 6.0) * rho[[lx, i]] * u[[0, lx, i]]
            + (1.0 / 2.0) * rho[[lx, i]] * u[[1, lx, i]];
    }
}

pub fn zou_he_top_wall_velocity(
    ly: usize,
    u: &mut Array3<f64>,
    u_top: &Array2<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    for i in 0..u.shape()[1] {
        u[[0, i, ly]] = u_top[[0, i]];
        u[[1, i, ly]] = u_top[[1, i]];

        rho[[i, ly]] = (g[[0, i, ly]]
            + g[[1, i, ly]]
            + g[[2, i, ly]]
            + 2.0 * g[[3, i, ly]]
            + 2.0 * g[[5, i, ly]]
            + 2.0 * g[[7, i, ly]])
            / (1.0 + u[[1, i, ly]]);

        g[[4, i, ly]] = g[[3, i, ly]] - (2.0 / 3.0) * rho[[i, ly]] * u[[1, i, ly]];

        g[[8, i, ly]] = g[[7, i, ly]] - (1.0 / 2.0) * (g[[1, i, ly]] - g[[2, i, ly]])
            + (1.0 / 2.0) * rho[[i, ly]] * u[[0, i, ly]]
            - (1.0 / 6.0) * rho[[i, ly]] * u[[1, i, ly]];

        g[[6, i, ly]] = g[[5, i, ly]] + (1.0 / 2.0) * (g[[1, i, ly]] - g[[2, i, ly]])
            - (1.0 / 2.0) * rho[[i, ly]] * u[[0, i, ly]]
            - (1.0 / 6.0) * rho[[i, ly]] * u[[1, i, ly]];
    }
}

pub fn zou_he_bottom_wall_velocity(
    u: &mut Array3<f64>,
    u_bot: &Array2<f64>,
    rho: &mut Array2<f64>,
    g: &mut Array3<f64>,
) {
    for i in 0..u.shape()[1] {
        u[[0, i, 0]] = u_bot[[0, i]];
        u[[1, i, 0]] = u_bot[[1, i]];

        rho[[i, 0]] = (g[[0, i, 0]]
            + g[[1, i, 0]]
            + g[[2, i, 0]]
            + 2.0 * g[[4, i, 0]]
            + 2.0 * g[[6, i, 0]]
            + 2.0 * g[[8, i, 0]])
            / (1.0 - u[[1, i, 0]]);

        g[[3, i, 0]] = g[[4, i, 0]] + (2.0 / 3.0) * rho[[i, 0]] * u[[1, i, 0]];

        g[[5, i, 0]] = g[[6, i, 0]] - (1.0 / 2.0) * (g[[1, i, 0]] - g[[2, i, 0]])
            + (1.0 / 2.0) * rho[[i, 0]] * u[[0, i, 0]]
            + (1.0 / 6.0) * rho[[i, 0]] * u[[1, i, 0]];

        g[[7, i, 0]] = g[[8, i, 0]] + (1.0 / 2.0) * (g[[1, i, 0]] - g[[2, i, 0]])
            - (1.0 / 2.0) * rho[[i, 0]] * u[[0, i, 0]]
            + (1.0 / 6.0) * rho[[i, 0]] * u[[1, i, 0]];
    }
}

// pub fn zou_he_right_wall_pressure(lx, ly, u, rho_right, u_right, rho, g){
//   rho[lx,:] = rho_right[:]
//   u[1,lx,:] = u_right[1,:]

//   u[0,lx,:] = (g[0,lx,:] + g[3,lx,:] + g[4,lx,:] +
//                2.0*g[1,lx,:] + 2.0*g[5,lx,:] +
//                2.0*g[8,lx,:])/rho[lx,:] - 1.0

//   g[2,lx,:] = (g[1,lx,:] - (2.0/3.0)*rho[lx,:]*u[0,lx,:])

//   g[6,lx,:] = (g[5,lx,:] + (1.0/2.0)*(g[3,lx,:] - g[4,lx,:]) -
//                (1.0/6.0)*rho[lx,:]*u[0,lx,:] -
//                (1.0/2.0)*rho[lx,:]*u[1,lx,:] )

//   g[7,lx,:] = (g[8,lx,:] - (1.0/2.0)*(g[3,lx,:] - g[4,lx,:]) -
//                (1.0/6.0)*rho[lx,:]*u[0,lx,:] +
//                (1.0/2.0)*rho[lx,:]*u[1,lx,:] )
//   }
