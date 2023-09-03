use super::point_in_shape;
use crate::lbm_core;

use ndarray::*;

pub struct Shape {
    pub pts: Vec<[f32; 2]>,
}

pub fn build_circle(center: [f32; 2], radius: f32, npts: usize) -> Shape {
    let mut shape = Shape { pts: Vec::new() };
    for i in 0..npts {
        let t: f32 = 2.0 * std::f32::consts::PI * (i as f32) / (npts as f32);
        shape.pts.push([
            radius * f32::cos(t) + center[0],
            radius * f32::sin(t) + center[1],
        ]);
    }
    return shape;
}

pub fn build_square(center: [f32; 2], side: f32) -> Shape {
    let mut shape = Shape { pts: Vec::new() };

    let x1 = center[0] - side / 2.0;
    let x2 = center[0] + side / 2.0;
    let y1 = center[1] - side / 2.0;
    let y2 = center[1] + side / 2.0;

    shape.pts.push([x1, y1]);
    shape.pts.push([x2, y1]);
    shape.pts.push([x2, y2]);
    shape.pts.push([x1, y2]);
    return shape;
}

fn yt_naca4(x: f32, t: f32) -> f32 {
    assert!(x >= 0.0 && x <= 1.0);
    // if close to 1.0 we impose 0.0 for symmetry
    if f32::abs(x - 1.0) < 2.0 * f32::EPSILON {
        return 0.0;
    }
    return 5.0
        * t
        * (0.2969 * f32::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x
            - 0.1015 * x * x * x * x);
}

pub fn build_naca4_sym(xx: &str, npts: usize) -> Shape {
    assert!(xx.len() == 2);
    let tt: u32 = xx.parse().unwrap();
    let t: f32 = (tt as f32) / 100.;
    let mut shape = Shape { pts: Vec::new() };

    let n = npts / 2;
    for i in 0..(n + 1) {
        let x = (i as f32) / (n as f32);
        let y = -yt_naca4(x, t);
        shape.pts.push([x, y]);
    }

    for i in 0..n - 1 {
        shape
            .pts
            .push([shape.pts[n - i - 1][0], -shape.pts[n - i - 1][1]]);
    }

    return shape;
}

fn yc_naca4(x: f32, p: f32, m: f32) -> f32 {
    assert!(x >= 0.0 && x <= 1.0);
    if x <= p {
        return m * (2.0 * p * x - x * x) / p / p;
    } else {
        return m * ((1. - 2. * p) + 2. * p * x - x * x) / (1. - p) / (1. - p);
    }
}

fn dyc_naca4(x: f32, p: f32, m: f32) -> f32 {
    assert!(x >= 0.0 && x <= 1.0);
    if x <= p {
        return 2. * m * (p - x) / p / p;
    } else {
        return 2. * m * (p - x) / (1. - p) / (1. - p);
    }
}

pub fn build_naca4(xx: &str, npts: usize) -> Shape {
    assert!(xx.len() == 4);
    let tt: u32 = xx[2..4].parse().unwrap();
    let pp: u32 = xx[1..2].parse().unwrap();
    let mm: u32 = xx[0..1].parse().unwrap();
    let t: f32 = (tt as f32) / 100.;
    let p: f32 = (pp as f32) / 10.;
    let m: f32 = (mm as f32) / 100.;

    let mut shape = Shape { pts: Vec::new() };

    let n = npts / 2;
    for i in 0..(n + 1) {
        let x = (i as f32) / (n as f32);
        let yt = yt_naca4(x, t);
        let yc = yc_naca4(x, p, m);
        let theta = f32::atan(dyc_naca4(x, p, m));
        shape
            .pts
            .push([x + yt * f32::sin(theta), yc - yt * f32::cos(theta)]);
    }

    for i in 0..n - 1 {
        let x = ((n - i - 1) as f32) / (n as f32);
        let yt = yt_naca4(x, t);
        let yc = yc_naca4(x, p, m);
        let theta = f32::atan(dyc_naca4(x, p, m));
        shape
            .pts
            .push([x - yt * f32::sin(theta), yc + yt * f32::cos(theta)]);
    }

    return shape;
}

fn compute_shape_bounding_box(shape: &Shape) -> [[f32; 2]; 2] {
    let mut bbox: [[f32; 2]; 2] = [[f32::MAX, f32::MIN], [f32::MAX, f32::MIN]];
    for i in 0..shape.pts.len() {
        for d in 0..2 {
            if shape.pts[i][d] < bbox[0][d] {
                bbox[0][d] = shape.pts[i][d];
            }
            if shape.pts[i][d] > bbox[1][d] {
                bbox[1][d] = shape.pts[i][d];
            }
        }
    }
    return bbox;
}

pub fn intersect_lattice_and_shape(
    x_min: f32,
    y_min: f32,
    dx: f32,
    dy: f32,
    nx: usize,
    ny: usize,
    shape: &Shape,
) -> (
    Vec<[u32; 2]>,
    Vec<[usize; 3]>,
    Vec<f64>,
    ArrayBase<OwnedRepr<u8>, Dim<[usize; 2]>>,
) {
    let mut obs: Vec<[u32; 2]> = Vec::new();
    let mut bnd: Vec<[usize; 3]> = Vec::new();
    let mut ibb: Vec<f64> = Vec::new();
    let mut tag: ArrayBase<OwnedRepr<u8>, Dim<[usize; 2]>> = Array2::<u8>::zeros((nx, ny));

    let bbox = compute_shape_bounding_box(&shape);

    let ilen = obs.len();

    for i in 0..nx {
        for j in 0..ny {
            let pt: [f32; 2] = [x_min + (i as f32) * dx, y_min + (j as f32) * dy];
            if (pt[0] > bbox[0][0])
                && (pt[0] < bbox[1][0])
                && (pt[1] > bbox[0][1])
                && (pt[1] < bbox[1][1])
            {
                if point_in_shape::cn_pn_poly(&pt, &shape) != 0 {
                    tag[[i, j]] = 1;
                    obs.push([i as u32, j as u32]);
                }
            }
        }
    }

    let lx = (nx - 1) as i32;
    let ly = (ny - 1) as i32;
    for k in ilen..obs.len() {
        let i = obs[k][0] as i32;
        let j = obs[k][1] as i32;
        for q in 1..9 {
            let qb = lbm_core::constants::NS[q];
            let cx = lbm_core::constants::CX[q];
            let cy = lbm_core::constants::CY[q];
            let ii = i + (cx as i32);
            let jj = j + (cy as i32);

            if (ii > lx) || (jj > ly) || (ii < 0) || (jj < 0) {
                continue;
            }
            if tag[[ii as usize, jj as usize]] != 0 {
                bnd.push([ii as usize, jj as usize, qb]);

                let pt: [f32; 2] = [x_min + (i as f32) * dx, y_min + (j as f32) * dy];

                let mut min_d2 = f32::MAX;
                for k in 0..shape.pts.len() {
                    let dx = shape.pts[k][0] - pt[0];
                    let dy = shape.pts[k][1] - pt[1];
                    let d2 = dx * dx + dy * dy;
                    if d2 < min_d2 {
                        min_d2 = d2;
                    }
                }
                let cx = (lbm_core::constants::CX[qb] as f32) * dx;
                let cy = (lbm_core::constants::CY[qb] as f32) * dy;
                ibb.push((min_d2 / f32::sqrt(cx * cx + cy * cy)) as f64);
            }
        }
    }
    return (obs, bnd, ibb, tag);
}

pub fn rotate(degrees: f32, offset: &[f32; 2], shape: &mut Shape) {
    let phi = degrees * std::f32::consts::PI / 180.;
    let rot: [[f32; 2]; 2] = [
        [f32::cos(phi), f32::sin(phi)],
        [-f32::sin(phi), f32::cos(phi)],
    ];
    for i in 0..shape.pts.len() {
        let tmp = [shape.pts[i][0] - offset[0], shape.pts[i][1] - offset[1]];
        for j in 0..2 {
            shape.pts[i][j] = tmp[0] * rot[j][0] + tmp[1] * rot[j][1] + offset[j];
        }
    }
}

pub fn barycenter(shape: &Shape) -> [f32; 2] {
    let mut r = [0.0f32, 0.];
    for i in 0..shape.pts.len() {
        r[0] += shape.pts[i][0];
        r[1] += shape.pts[i][1];
    }
    r[0] /= shape.pts.len() as f32;
    r[1] /= shape.pts.len() as f32;
    return r;
}
