use rust_lbm;

use macroquad::prelude::*;

#[macroquad::main("Macroquad test")]
async fn main() {
    let mut render_settings = rust_lbm::render::RenderSettings {
        iter_per_frame: 100,
        plot_mode: rust_lbm::render::Plot::VelocityNorm,
        show_info: true,
        should_restart: true,
        target: (screen_width() / 2.0, screen_height() / 2.0),
        zoom: 2.0 / screen_width(),
    };

    let circle = rust_lbm::geometry::shape::build_circle([0.0, 0.0], 0.1, 200);
    let naca15 = rust_lbm::geometry::shape::build_naca4_sym("15", 100);
    let naca2412 = rust_lbm::geometry::shape::build_naca4("2412", 1001);

    let mut lattice = rust_lbm::lbmcore::turek::get_lattice();

    let shape = &naca2412;
    let (obs, bnd, ibb) =
        rust_lbm::geometry::shape::intersect_lattice_and_shape(&mut lattice, shape);
    //
    let mut it: usize = 0;

    println!("{}, {}", lattice.nx, lattice.ny);
    println!("{}, {}", lattice.dx, lattice.dy);

    loop {
        rust_lbm::render::handle_event(&mut render_settings);

        // if render_settings.should_restart {
        //     lattice = rust_lbm::lbmcore::cavity::get_lattice();
        //     rust_lbm::lbmcore::cavity::init_lattice(&mut lattice);
        //     it = 0;
        //     render_settings.should_restart = false;
        // }
        // for _ in 0..render_settings.iter_per_frame {
        //     rust_lbm::lbmcore::cavity::step_lattice(&mut lattice, it);
        //     it = it + 1;
        // }
        // rust_lbm::render::draw_velocity(&lattice, &render_settings);

        set_camera(&Camera2D {
            target: vec2(render_settings.target.0, render_settings.target.1),
            zoom: vec2(
                render_settings.zoom,
                render_settings.zoom * screen_width() / screen_height(),
            ),
            ..Default::default()
        });

        rust_lbm::render::draw_lattice(&lattice);
        rust_lbm::render::draw_shape(shape, &lattice, &obs);

        set_default_camera();
        rust_lbm::render::draw_info(it, &render_settings);
        next_frame().await;
    }
}
