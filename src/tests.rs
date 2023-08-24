#[cfg(test)]
mod tests {
    use crate::lbmcore;

    #[test]
    fn run_cavity_test() {
        use std::time::Instant;
        let start_time = Instant::now();

        let mut lattice = lbmcore::cavity::get_lattice();
        lbmcore::cavity::init_lattice(&mut lattice);

        for it in 0..lattice.it_max {
            if it % 100 == 0 {
                println!("Iteration: {:>6}/{}", it, lattice.it_max);
            }
            lbmcore::cavity::step_lattice(&mut lattice, it);
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
}