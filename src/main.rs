#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]

pub mod absolute_value_tract;
pub mod baker_tract;
pub mod cone;
pub mod cone_errors;
pub mod examples;
pub mod fan;
pub mod field_tract;
mod integer_arith;
pub mod krasner;
pub mod linear_combination;
pub mod motzkin;
pub mod phase_tract;
pub mod polytope;
pub mod regular_partial_field;
pub mod toric_ideal;
pub mod tract_grassmannian;
pub mod tract_projective_variety;
pub mod tropical_tract;

use crate::{examples::main_examples, toric_ideal::main_toric_ideal_example};

fn main() {
    main_examples();
    main_toric_ideal_example();
}
