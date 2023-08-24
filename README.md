
# Lattice Boltzmann method in Rust with WASM translation
![Rust badge](https://github.com/MatteoCaldana/rust-lbm/workflows/RustCavityBenchmark/badge.svg)

<!-- <h3 align="center">
	<a href="https://matteocaldana.github.io/rust-lbm/">
		Online Demo ☁
	</a>
</h2> -->


## Usage

### Main
To run the code locally, clone the project and run:
```bash
cargo run --release
```

This will execute the [visualization demo](./src/bin/main.rs) locally.

### Tests
To verify the correctness of the code we compare the velocity field of the lid-driven cavity benchmark at the center of the domain with reference data from ["U. Ghia, K. N. Ghia, C. T. Shin, High-Re solutions for incompressible flow using Navier-Stokes equations and multigrid method"](https://www.sciencedirect.com/science/article/pii/0021999182900584).

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
cp target/release/wasm32-unknown-unknown/main.wasm pages
```
One of the ways to then server static `.wasm` and `.html` (blocked by CORS otherwise):
```bash
cargo install basic-http-server
cd pages
basic-http-server .
```

## Motivation

In this small project I had the opportuninty to learn about
* Rust: a modern high-performace language that is gain a lot of popularity
* WASM: a binary instruction format that, albeit some penalty, makes it possible to run performance critical code in a web browser.
* LBM: an application of Boltzmann particle (microscopic) principles in a lattice grid to simulate computation fluid dynamics (CFD) without directly solving the Navier-Stokes equations.

## References
* https://github.com/ndbaker1/bloe
* https://github.com/jviquerat/lbm