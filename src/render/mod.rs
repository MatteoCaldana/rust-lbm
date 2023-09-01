use macroquad::prelude::*;
use ndarray::*;
use std::fmt;

use crate::geometry;
use crate::lbmcore::{self, lattice};

#[allow(dead_code)]
pub enum Plot {
    Density,
    VelocityX,
    VelocityY,
    VelocityNorm,
    Vorticity,
}

impl fmt::Display for Plot {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Plot::Density => write!(f, "Density"),
            Plot::VelocityX => write!(f, "VelocityX"),
            Plot::VelocityY => write!(f, "VelocityY"),
            Plot::VelocityNorm => write!(f, "VelocityNorm"),
            Plot::Vorticity => write!(f, "Vorticity"),
        }
    }
}

pub struct RenderSettings {
    pub iter_per_frame: usize,
    pub plot_mode: Plot,
    pub show_info: bool,
    pub should_restart: bool,
    pub target: (f32, f32),
    pub zoom: f32,
}

pub fn handle_event(render_settings: &mut RenderSettings) {
    if is_key_pressed(KeyCode::R) {
        render_settings.should_restart = true;
    }
    if is_key_pressed(KeyCode::C) {
        render_settings.show_info = !render_settings.show_info;
    }
    // camera pan
    const PAN: f32 = 0.1;
    if is_key_down(KeyCode::Up) {
        render_settings.target.1 -= PAN / render_settings.zoom;
    }
    if is_key_down(KeyCode::Down) {
        render_settings.target.1 += PAN / render_settings.zoom;
    }
    if is_key_down(KeyCode::Left) {
        render_settings.target.0 += PAN / render_settings.zoom;
    }
    if is_key_down(KeyCode::Right) {
        render_settings.target.0 -= PAN / render_settings.zoom;
    }
    // camera zoom
    match mouse_wheel() {
        (_x, y) if y != 0.0 => {
            // Normalize mouse wheel values is browser (chromium: 53, firefox: 3)
            #[cfg(target_arch = "wasm32")]
            let y = if y < 0.0 {
                -1.0
            } else if y > 0.0 {
                1.0
            } else {
                0.0
            };
            render_settings.zoom *= 1.1f32.powf(y);
            
        }
        _ => (),
    }
    // plot
    if is_key_pressed(KeyCode::X) {
        render_settings.plot_mode = Plot::VelocityX;
    }
    if is_key_pressed(KeyCode::Y) {
        render_settings.plot_mode = Plot::VelocityY;
    }
    if is_key_pressed(KeyCode::N) {
        render_settings.plot_mode = Plot::VelocityNorm;
    }
    if is_key_pressed(KeyCode::W) {
        render_settings.plot_mode = Plot::Vorticity;
    }
    if is_key_pressed(KeyCode::D) {
        render_settings.plot_mode = Plot::Density;
    }
    //
    const FRAME_DELTA: usize = 5;
    if is_key_pressed(KeyCode::PageUp) {
        render_settings.iter_per_frame = render_settings.iter_per_frame + FRAME_DELTA;
    }
    if is_key_pressed(KeyCode::PageDown) {
        if render_settings.iter_per_frame > FRAME_DELTA {
            render_settings.iter_per_frame = render_settings.iter_per_frame - FRAME_DELTA;
        }
    }
    //
    #[cfg(not(target_arch = "wasm32"))]
    if is_key_pressed(KeyCode::Escape) {
        std::process::exit(0);
    }
}

fn draw_field<T: ndarray::RawData<Elem = f64> + ndarray::Data>(
    field: &ArrayBase<T, Dim<[usize; 2]>>,
) {
    let w = screen_width();
    let h = screen_height();
    let nx = field.shape()[0] as f32;
    let ny = field.shape()[1] as f32;
    let dx = w / nx;
    let dy = h / ny;
    let mut fmin = f64::MAX;
    let mut fmax = f64::MIN;
    for i in 0..field.shape()[0] {
        for j in 0..field.shape()[1] {
            let x: f64 = field[[i, j]];
            if x < fmin {
                fmin = x;
            }
            if x > fmax {
                fmax = x;
            }
        }
    }

    let gradient = colorous::MAGMA;

    for i in 0..field.shape()[0] {
        for j in 0..field.shape()[1] {
            let color = gradient.eval_continuous((field[[i, j]] - fmin) / (fmax - fmin));
            draw_rectangle(
                (i as f32) * dx,
                (j as f32) * dy,
                dx,
                dy,
                Color::from_rgba(color.r, color.g, color.b, 255u8),
            );
        }
    }
}

pub fn draw_velocity(lattice: &lattice::Lattice, settings: &RenderSettings) {
    if let Plot::VelocityX = settings.plot_mode {
        draw_field(&lattice.u.index_axis(ndarray::Axis(0), 0));
    } else if let Plot::VelocityY = settings.plot_mode {
        draw_field(&lattice.u.index_axis(ndarray::Axis(0), 1));
    } else if let Plot::VelocityNorm = settings.plot_mode {
        let mut u2: ndarray::Array2<f64> =
            lattice.u.index_axis(ndarray::Axis(0), 0).mapv(|a| a * a);
        // TODO: do this is a more rust-like manner
        for i in 0..u2.shape()[0] {
            for j in 0..u2.shape()[1] {
                u2[[i, j]] = f64::sqrt(u2[[i, j]] + lattice.u[[1, i, j]] * lattice.u[[1, i, j]]);
            }
        }
        draw_field(&u2);
    } else if let Plot::Vorticity = settings.plot_mode {
        let mut w = Array2::<f64>::zeros((lattice.u.shape()[1] - 2, lattice.u.shape()[2] - 2));
        // TODO: check implementation is correct
        for i in 0..w.shape()[0] {
            for j in 0..w.shape()[1] {
                w[[i, j]] = (lattice.u[[0, i, j + 1]] - lattice.u[[0, i + 2, j + 1]])
                    - (lattice.u[[1, i + 1, j]] - lattice.u[[1, i + 1, j + 2]]);
            }
        }
        draw_field(&w);
    } else if let Plot::Density = settings.plot_mode {
        draw_field(&lattice.rho);
    }
}

pub fn draw_info(it: usize, settings: &RenderSettings) {
    const FONT_SIZE: f32 = 25.;
    if settings.show_info {
        let infos = [
            "Toggle [C]ommands. [R]estart. Arrows to move around, scroll to zoom.",
            &format!("Iteration: {:>6}, FPS: {:>3}", it, get_fps()),
            &format!(
                "Iter/frame: {} ([PageUp]/[PageDown] to change)",
                settings.iter_per_frame
            ),
            &format!(
                "Showing: {} ([D]ensity, Velocity[X]/[Y]/[N]orm, [W]orticity)",
                settings.plot_mode
            ),
        ];
        for i in 0..infos.len() {
            draw_text(infos[i], 0., (i + 1) as f32 * FONT_SIZE, FONT_SIZE, RED);
        }
    }
}

pub fn draw_shape(shape: &geometry::shape::Shape, lattice: &lbmcore::lattice::Lattice, obs: &Vec<[u32; 2]>) {
    let x_min: f32 = lattice.x_min as f32;
    let x_max: f32 = lattice.x_max as f32;
    let y_min: f32 = lattice.y_min as f32;
    let y_max: f32 = lattice.y_max as f32;

    let w = screen_width();
    let h = screen_height();
    let x_tr = w / (x_max - x_min);
    let y_tr = h / (y_max - y_min);
    for i in 0..shape.pts.len() - 1 {
        let x1 = x_tr * (shape.pts[i][0] - x_min);
        let x2 = x_tr * (shape.pts[i + 1][0] - x_min);
        let y1 = y_tr * (shape.pts[i][1] - y_min);
        let y2 = y_tr * (shape.pts[i + 1][1] - y_min);
        draw_line(x1, y1, x2, y2, 1.0, WHITE);
    }
    let x1 = x_tr * (shape.pts[shape.pts.len() - 1][0] - x_min);
    let x2 = x_tr * (shape.pts[0][0] - x_min);
    let y1 = y_tr * (shape.pts[shape.pts.len() - 1][1] - y_min);
    let y2 = y_tr * (shape.pts[0][1] - y_min);
    draw_line(x1, y1, x2, y2, 1.0, WHITE);

    for k in 0..obs.len() {
        let x = (obs[k][0] as f32) * (lattice.dx as f32);
        let y = (obs[k][1] as f32) * (lattice.dy as f32);
        draw_circle(x_tr * x , y_tr * y, 1.0, RED);
    }
}

pub fn draw_lattice(lattice: &lbmcore::lattice::Lattice) {
    let w = screen_width();
    let h = screen_height();
    for i in 0..lattice.nx {
        let x = w * (i as f32) / ((lattice.nx - 1) as f32);
        draw_line(x, 0.0, x, h, 1.0, GRAY);
    }
    for j in 0..lattice.ny {
        let y = h * (j as f32) / ((lattice.ny - 1) as f32);
        draw_line(0.0 , y, w, y, 1.0, GRAY);
    }
}
