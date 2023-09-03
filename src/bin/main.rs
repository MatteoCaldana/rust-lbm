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
    //let naca15 = rust_lbm::geometry::shape::build_naca4_sym("15", 1000);
    //let naca6409 = rust_lbm::geometry::shape::build_naca4("6409", 1001);
    let shape = &circle;

    let mut lattice = rust_lbm::lbm_app::turek::get_lattice();
    rust_lbm::geometry::shape::intersect_lattice_and_shape(&mut lattice, shape);
    //
    let mut it: usize = 0;
    println!("u_lbm: {}", lattice.u_lbm);
    println!("nx/ny: {}, {}", lattice.nx, lattice.ny);
    println!("dx/dy: {}, {}, {}", lattice.dx, lattice.dy, lattice.dt);
    println!("tau:   {}, {}", lattice.tau_p_lbm, lattice.tau_m_lbm);
    println!("om:    {}, {}", lattice.om_p_lbm, lattice.om_m_lbm);

    loop {
        println!("time: {}", lattice.dt * (it as f64));
        rust_lbm::render::handle_event(&mut render_settings);
        // check if need to restart
        if render_settings.should_restart {
            lattice = rust_lbm::lbm_app::turek::get_lattice();
            rust_lbm::lbm_core::lattice::init_lattice(
                &mut lattice,
                rust_lbm::lbm_app::poiseuille::set_inlets,
            );
            rust_lbm::geometry::shape::intersect_lattice_and_shape(&mut lattice, shape);
            it = 0;
            render_settings.should_restart = false;
        }
        // simulate
        if let rust_lbm::render::Plot::Grid = render_settings.plot_mode {
        } else {
            for _ in 0..render_settings.iter_per_frame {
                rust_lbm::lbm_core::lattice::step_lattice(
                    &mut lattice,
                    it,
                    rust_lbm::lbm_app::poiseuille::set_inlets,
                    rust_lbm::lbm_app::turek::apply_bc,
                    |_| (),
                );
                it = it + 1;
            }
            //TODO: compute drag and lift
        }

        // draw physical
        set_camera(&Camera2D {
            target: vec2(render_settings.target.0, render_settings.target.1),
            zoom: vec2(
                render_settings.zoom,
                render_settings.zoom * screen_width() / screen_height(),
            ),
            ..Default::default()
        });
        rust_lbm::render::draw(&lattice, shape, &render_settings);

        // draw text
        set_default_camera();
        rust_lbm::render::draw_info(it, &render_settings);

        next_frame().await;
    }
}
