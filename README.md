
# Lattice Boltzmann method in Rust with WASM translation
![RustPhysicsBenchmark badge](https://github.com/MatteoCaldana/rust-lbm/workflows/RustPhysicsBenchmark/badge.svg)
![DeployWASM badge](https://github.com/MatteoCaldana/rust-lbm/workflows/DeployWASM/badge.svg)

<h3 align="center">
	<a href="https://matteocaldana.github.io/rust-lbm/">
		Online Demo ☁
	</a>
</h2>


## Usage

### Main
To run the code locally, clone the project and run:
```bash
cargo run --release
```

This will execute the [visualization demo](./src/bin/main.rs) locally.

### Tests
To verify the correctness of the code we use the following unit tests:
1. Compare the velocity field of the lid-driven cavity benchmark at the center of the domain with reference data from ["U. Ghia, K. N. Ghia, C. T. Shin, High-Re solutions for incompressible flow using Navier-Stokes equations and multigrid method"](https://www.sciencedirect.com/science/article/pii/0021999182900584).
2. Run a channel flow simulaiton at low Reynolds (Poiseuille flow) and check the velocity profile at the outlet.
3. **[WIP]** Computation of drag and lift in the Turek cylinder benchmark CFD case is described in ["Schafer, M., Turek, S. Benchmark Computations of Laminar Flow Around a Cylinder"](https://link.springer.com/chapter/10.1007/978-3-322-89849-4_39).

```bash
cargo test --release
```

### WASM
To translate the code to WASM follow the [Macroquad instructions](https://github.com/not-fl3/macroquad#wasm)
```bash
rustup target add wasm32-unknown-unknown
cargo build --target wasm32-unknown-unknown --release
```
This will produce `.wasm` file in `target/debug/wasm32-unknown-unknown/main.wasm` or in `target/release/wasm32-unknown-unknown/main.wasm` if built with `--release`. Copy it into the `pages` folder.
```bash
cp target/wasm32-unknown-unknown/release/main.wasm pages
```
One of the ways to then server static `.wasm` and `.html` (blocked by CORS otherwise):
```bash
cargo install basic-http-server
cd pages
basic-http-server .
```

## Motivation

In this small project I had the opportuninty to learn about:
* Rust: a modern high-performace language
* WASM: a binary instruction format that, albeit some penalty, makes it possible to run performance critical code in a web browser.
* LBM: an application of Boltzmann particle (microscopic) principles in a lattice grid to simulate computation fluid dynamics (CFD) without directly solving the Navier-Stokes equations.
* Computational geometry algorithms to raster a polygon

### Possible Extensions
Add (experimental) WASM features:
* Use SIMD instructions
* Parition the lattice and use threads to parallelize (spawing a thread is very expensive in WASM, must be done carefully)

## References
* https://github.com/ndbaker1/bloe
* https://github.com/jviquerat/lbm
* https://en.wikipedia.org/wiki/NACA_airfoil