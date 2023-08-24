use macroquad::prelude::*;
use ndarray::*;
use std::fmt;

use crate::lbmcore::lattice;


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
    pub should_restart: bool
} 


pub fn handle_event(render_settings: &mut RenderSettings) {
    if is_key_pressed(KeyCode::R) {
        render_settings.should_restart = true;
    }
    if is_key_pressed(KeyCode::C) {
        render_settings.show_info = !render_settings.show_info;
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
    if is_key_pressed(KeyCode::Up) {
        render_settings.iter_per_frame = render_settings.iter_per_frame + FRAME_DELTA;
    }
    if is_key_pressed(KeyCode::Down) {
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

    let gradient = colorous::VIRIDIS;
    
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


pub fn draw_lattice(lattice: &lattice::Lattice, settings: &RenderSettings) {
    if let Plot::VelocityX = settings.plot_mode{
        draw_field(&lattice.u.index_axis(ndarray::Axis(0), 0));
    } else if let Plot::VelocityY = settings.plot_mode{ 
        draw_field(&lattice.u.index_axis(ndarray::Axis(0), 1));
    } else if let Plot::VelocityNorm = settings.plot_mode {
        let mut u2: ndarray::Array2<f64> =  lattice.u.index_axis(ndarray::Axis(0), 0).mapv(|a| a*a);
        // TODO: do this is a more rust-like manner
        for i in 0..u2.shape()[0]{
            for j in 0..u2.shape()[1]{
                u2[[i, j]] = f64::sqrt(u2[[i, j]] + lattice.u[[1, i, j]] * lattice.u[[1, i, j]]);
            }
        }
        draw_field(&u2);
    }   else if let Plot::Vorticity = settings.plot_mode {
        let mut w= Array2::<f64>::zeros((lattice.u.shape()[1] - 2, lattice.u.shape()[2] - 2));
        for i in 0..w.shape()[0] {
            for j in 0..w.shape()[1] {
                w[[i, j]] = (lattice.u[[0, i, j+1]] - lattice.u[[0, i + 2, j+1]]) - (lattice.u[[1, i+1, j]] - lattice.u[[1, i+1, j+2]]);
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
            "Toggle [C]ommands. [R]estart.",
            &format!("Iteration: {:>6}, FPS: {:>3}", it, get_fps()),
            &format!("Iter/frame {} ([Up]/[Down] to change)", settings.iter_per_frame),
            &format!("Showing {} ([D]ensity, Velocity[X]/[Y]/[N]orm, [W]orticity)", settings.plot_mode),
        ];
        for i in 0..infos.len() {
            draw_text(infos[i], 0., (i + 1) as f32 * FONT_SIZE, FONT_SIZE,RED);
        }
    }
}