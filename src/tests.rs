#[cfg(test)]
mod tests {
    use crate::lbm_core;
    use crate::lbm_app;

    #[test]
    fn run_cavity_test() {
        use std::time::Instant;
        let start_time = Instant::now();

        let mut lattice = lbm_app::cavity::get_lattice();
        lbm_core::lattice::init_lattice(&mut lattice, lbm_app::cavity::set_inlets);

        for it in 0..lattice.it_max {
            if it % 100 == 0 {
                println!("Iteration: {:>6}/{}", it, lattice.it_max);
            }
            lbm_core::lattice::step_lattice(
                &mut lattice,
                it,
                lbm_app::cavity::set_inlets,
                lbm_app::cavity::apply_bc,
                |_| (),
            );
        }

        let elapsed = start_time.elapsed();
        println!("Elapsed: {:.2?}", elapsed);

        assert!(f64::abs(lattice.u[[0, 50, 10]] / lattice.u_lbm + 0.05968571) < 1.0e-6);
        assert!(f64::abs(lattice.u[[0, 50, 50]] / lattice.u_lbm + 0.19792323) < 1.0e-6);
        assert!(f64::abs(lattice.u[[0, 50, 90]] / lattice.u_lbm - 0.45278683) < 1.0e-6);
        assert!(f64::abs(lattice.u[[1, 10, 50]] / lattice.u_lbm - 0.12493089) < 1.0e-6);
        assert!(f64::abs(lattice.u[[1, 50, 50]] / lattice.u_lbm - 0.05295104) < 1.0e-6);
        assert!(f64::abs(lattice.u[[1, 90, 50]] / lattice.u_lbm + 0.16641426) < 1.0e-6);
    }


    #[test]
    fn run_poiseuille_test() {
        use std::time::Instant;
        let start_time = Instant::now();

        let mut lattice = lbm_app::poiseuille::get_lattice();
        lbm_core::lattice::init_lattice(&mut lattice, lbm_app::poiseuille::set_inlets);

        for it in 0..lattice.it_max {
            if it % 100 == 0 {
                println!("Iteration: {:>6}/{}", it, lattice.it_max);
            }
            lbm_core::lattice::step_lattice(
                &mut lattice,
                it,
                lbm_app::poiseuille::set_inlets,
                lbm_app::poiseuille::apply_bc,
                |_| (),
            );
        }

        let elapsed = start_time.elapsed();
        println!("Elapsed: {:.2?}", elapsed);

        let mut err_ux_l1 = 0.0;
        let mut err_uy_l1 = 0.0;
        for j in 0..lattice.u.shape()[2] {
            err_ux_l1 += f64::abs(lattice.u[[0, lattice.lx, j]] - lattice.u[[0, 0, j]]);
            err_uy_l1 += f64::abs(lattice.u[[1, lattice.lx, j]]);
        }
        println!("err ux: {}", err_ux_l1);
        println!("err uy: {}", err_uy_l1);

        assert!((err_ux_l1 < 0.025));
        assert!((err_uy_l1 < 1e-6));


    }
}
