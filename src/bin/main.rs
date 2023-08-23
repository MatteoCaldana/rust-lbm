use rust_lbm;

use macroquad::prelude::*;

#[macroquad::main("Macroquad test")]
async fn main() {
    let mut render_settings = rust_lbm::render::RenderSettings {
        iter_per_frame: 100,
        plot_mode: rust_lbm::render::Plot::VelocityNorm,
        show_info: true, 
        should_restart: true
    };

    let mut lattice = rust_lbm::get_cavity_lattice();
    let mut it: usize = 0;
    loop {
        rust_lbm::render::handle_event(&mut render_settings);

        if render_settings.should_restart {
            lattice = rust_lbm::get_cavity_lattice();
            rust_lbm::init_lattice(&mut lattice);
            it = 0;
            render_settings.should_restart = false;
        }
        for _ in 0..render_settings.iter_per_frame {
            rust_lbm::step_lattice(&mut lattice, it);
            it = it + 1;
        }
        rust_lbm::render::draw_lattice(&lattice, &render_settings);
        rust_lbm::render::draw_info(it, &render_settings);
        next_frame().await;
    }
    
}