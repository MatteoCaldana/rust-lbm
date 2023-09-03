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

    let mut lattice = rust_lbm::lbm_app::web::get_lattice();
    //
    let mut it: usize = 0;
    println!("u_lbm: {}", lattice.u_lbm);
    println!(
        "nx/ny: {}, {}, {}",
        lattice.nx,
        lattice.ny,
        lattice.nx * lattice.ny
    );
    println!("dx/dy: {}, {}, {}", lattice.dx, lattice.dy, lattice.dt);
    println!("tau:   {}, {}", lattice.tau_p_lbm, lattice.tau_m_lbm);
    println!("om:    {}, {}", lattice.om_p_lbm, lattice.om_m_lbm);

    let (mut drag, mut lift) = (0., 0.);

    loop {
        println!("time: {}", lattice.dt * (it as f64));
        rust_lbm::render::handle_event(&mut render_settings);
        // check if need to restart
        if render_settings.should_restart {
            lattice = rust_lbm::lbm_app::web::get_lattice();
            rust_lbm::lbm_core::lattice::init_lattice(
                &mut lattice,
                rust_lbm::lbm_app::poiseuille::set_inlets,
            );
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
                    rust_lbm::lbm_app::web::apply_bc,
                    |_| (),
                );
                it = it + 1;
            }
            (drag, lift) = rust_lbm::lbm_core::lattice::compute_drag_lift(
                &lattice.bnd,
                &lattice.g_up,
                &lattice.g,
                lattice.rho_lbm,
                2.0 * lattice.u_lbm / 3.0,
                (lattice.ny as f64) * 0.1 / (lattice.y_max - lattice.y_min),
            );
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
        rust_lbm::render::draw(&lattice, &render_settings);
        if let rust_lbm::render::Plot::Grid = render_settings.plot_mode {
        } else {
            rust_lbm::render::draw_drag_lift(drag as f32, lift as f32, &lattice);
        }

        // draw text
        set_default_camera();
        rust_lbm::render::draw_info(it, &render_settings);

        next_frame().await;
    }
}
