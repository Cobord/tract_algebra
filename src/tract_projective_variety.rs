/// The tract analog of a projective variety cut out by homogeneous polynomials.
///
/// Over an honest field `F`, a set of homogeneous degree-`HOMOGENEITY` polynomials with
/// integer coefficients in `K_VARIABLES` variables defines an ideal whose vanishing locus
/// in `F^{K_VARIABLES}` is invariant under scaling by `F*`, giving a variety in `FP^{K-1}`.
///
/// Over a tract `T` there is no well-defined addition, only a null-set membership test.
/// Each polynomial `f = Σ_t a_t · x_{i_{1,t}} · … · x_{i_{d,t}}` (encoded as a
/// `LinearCombination<[usize; HOMOGENEITY], i64>`) is evaluated at a coordinate tuple
/// by computing the product of the relevant tract elements for each monomial, attaching
/// the sign via `T::neg_one()` and the magnitude as a multiplicity, then testing whether
/// the resulting element of `N[T]` lies in the null set.
///
/// Scaling all coordinates by a common `z ∈ T*` multiplies every monomial value by `z^HOMOGENEITY`,
/// and since the null set is invariant under the `T*`-action, the membership test is
/// projectively well-defined. The zero point (all coordinates absorbing) is excluded.
///
/// A monomial `[i_1, i_2, …, i_HOMOGENEITY]` denotes `x_{i_1} · x_{i_2} · … · x_{i_HOMOGENEITY}`;
/// repeated indices are allowed (they produce higher powers of a single variable).
/// All indices must be in `0..K_VARIABLES`.
use crate::{baker_tract::BakerBowlerTract, linear_combination::LinearCombination};

pub struct TractProjectiveVariety<
    const K_VARIABLES: usize,
    const HOMOGENEITY: usize,
    T: BakerBowlerTract,
> {
    coordinates: [T; K_VARIABLES],
    relations: Vec<LinearCombination<[usize; HOMOGENEITY], i64>>,
}

#[derive(Debug, thiserror::Error)]
pub enum TractProjectiveVarietyError {
    #[error("all coordinates are absorbing — the point does not lie in projective space")]
    AllAbsorbing,
    #[error("variable index {index} is out of range; expected < {k_variables}")]
    VariableOutOfRange { index: usize, k_variables: usize },
    #[error("relation {relation_index} is not satisfied at this point")]
    RelationNotSatisfied { relation_index: usize },
}

impl<const K_VARIABLES: usize, const HOMOGENEITY: usize, T: BakerBowlerTract>
    TractProjectiveVariety<K_VARIABLES, HOMOGENEITY, T>
{
    /// Construct a point on the variety, validating the relations unless `skip_valid` is set.
    ///
    /// Returns `Err` if all coordinates are absorbing, any monomial index is out of range,
    /// or (when `!skip_valid`) any relation fails the null-set test.
    pub fn new(
        coordinates: [T; K_VARIABLES],
        relations: Vec<LinearCombination<[usize; HOMOGENEITY], i64>>,
        skip_valid: bool,
    ) -> Result<Self, TractProjectiveVarietyError> {
        if coordinates.iter().all(T::is_absorbing_elt) {
            return Err(TractProjectiveVarietyError::AllAbsorbing);
        }
        // Variable-bounds check is cheap and always worth doing.
        for relation in &relations {
            for monomial in relation.map.keys() {
                for &var_idx in monomial {
                    if var_idx >= K_VARIABLES {
                        return Err(TractProjectiveVarietyError::VariableOutOfRange {
                            index: var_idx,
                            k_variables: K_VARIABLES,
                        });
                    }
                }
            }
        }
        let to_return = Self {
            coordinates,
            relations,
        };
        if skip_valid {
            return Ok(to_return);
        }
        for (relation_index, relation) in to_return.relations.iter().enumerate() {
            if !to_return.relation_satisfied(relation) {
                return Err(TractProjectiveVarietyError::RelationNotSatisfied { relation_index });
            }
        }
        Ok(to_return)
    }

    /// Returns the coordinate at position `idx` (cloned).
    #[must_use]
    pub fn get_coordinate(&self, idx: usize) -> T {
        assert!(idx < K_VARIABLES, "coordinate index out of range");
        self.coordinates[idx].clone()
    }

    /// Tests whether the stored point satisfies the given relation.
    ///
    /// A positive coefficient `a > 0` contributes `(monomial_value, a as usize)`.
    /// A negative coefficient `a < 0` contributes `(monomial_value * neg_one, |a| as usize)`.
    /// Zero coefficients are skipped. An empty sum (all coefficients zero) is trivially null.
    fn relation_satisfied(&self, relation: &LinearCombination<[usize; HOMOGENEITY], i64>) -> bool {
        let mut terms: Vec<(T, usize)> = Vec::new();
        for (monomial, &coeff) in &relation.map {
            if coeff == 0 {
                continue;
            }
            // Product of the tract coordinates named by the monomial.
            let mut val = T::one();
            for &var_idx in monomial {
                val *= self.coordinates[var_idx].clone();
            }
            if coeff < 0 {
                val *= T::neg_one();
            }
            terms.push((
                val,
                usize::try_from(coeff.unsigned_abs()).expect("The number is small"),
            ));
        }
        let as_null_set: Result<T::NullSet, _> = terms.try_into();
        if let Ok(ns) = as_null_set {
            T::in_null_set(ns)
        } else {
            true // empty sum is null
        }
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;

    use num::One;

    use crate::{
        absolute_value_tract::AbsoluteValueTract, baker_tract::BakerBowlerTract,
        field_tract::FieldTract, krasner::Krasner, linear_combination::LinearCombination,
    };

    use super::{TractProjectiveVariety, TractProjectiveVarietyError};

    // Coordinates: x = 0, y = 1, z = 2.
    // All curves live in projective 2-space: K_VARIABLES = 3.

    fn f(n: i128) -> FieldTract<i128> {
        FieldTract::from(n)
    }

    // -------------------------------------------------------------------------
    // Cuspidal cubic  x³ - y²z = 0  (HOMOGENEITY = 3)
    //
    // Relation: +1·[0,0,0]  −1·[1,1,2]
    // Parametric family of smooth points: [t² : t³ : 1] for any t.
    // -------------------------------------------------------------------------

    type Cubic = TractProjectiveVariety<3, 3, FieldTract<i128>>;

    fn cuspidal_cubic_relation() -> LinearCombination<[usize; 3], i64> {
        vec![([0, 0, 0], 1i64), ([1, 1, 2], -1i64)]
            .try_into()
            .expect("valid cuspidal cubic relation")
    }

    fn cubic_point(x: i128, y: i128, z: i128) -> Result<Cubic, TractProjectiveVarietyError> {
        Cubic::new([f(x), f(y), f(z)], vec![cuspidal_cubic_relation()], false)
    }

    #[test]
    fn cusp_at_origin_on_curve() {
        // [0:0:1]: x³ = 0, y²z = 0 — the cusp itself.
        assert!(cubic_point(0, 0, 1).is_ok());
    }

    #[test]
    fn point_at_infinity_on_curve() {
        // [0:1:0]: x³ = 0, y²z = 1·0 = 0.
        assert!(cubic_point(0, 1, 0).is_ok());
    }

    #[test]
    fn fixed_parametric_points_on_curve() {
        // [t²:t³:1] satisfies x³ = t⁶ = y²z for t = 1, 2, 3.
        assert!(cubic_point(1, 1, 1).is_ok());
        assert!(cubic_point(4, 8, 1).is_ok());
        assert!(cubic_point(9, 27, 1).is_ok());
    }

    #[test]
    fn off_curve_points_rejected() {
        // [1:0:0]: x³ = 1, y²z = 0.
        assert!(cubic_point(1, 0, 0).is_err());
        // [1:2:1]: x³ = 1, y²z = 4.
        assert!(cubic_point(1, 2, 1).is_err());
    }

    proptest! {
        // [t²:t³:1] lies on the cuspidal cubic for every integer t.
        // t⁶ fits comfortably in i128 for |t| ≤ 10⁶.
        #[test]
        fn parametric_family_on_curve(t in -1_000_000i64..=1_000_000i64) {
            let t = t as i128;
            prop_assert!(cubic_point(t * t, t * t * t, 1).is_ok());
        }

        // A generic off-curve point: x³ ≠ y²z whenever y²z − x³ ≠ 0.
        // We pick x arbitrary, y = x + 1 (distinct), z = 1 and verify algebraically.
        #[test]
        fn off_curve_when_not_equal(x in -500i64..=500i64) {
            let x = x as i128;
            let y = x + 1;
            let x3 = x * x * x;
            let y2z = y * y;
            prop_assume!(x3 != y2z);
            prop_assert!(cubic_point(x, y, 1).is_err());
        }
    }

    // -------------------------------------------------------------------------
    // Conic  z² - x² - y² = 0  (HOMOGENEITY = 2)
    //
    // Relation: +1·[2,2]  −1·[0,0]  −1·[1,1]
    // Pythagorean triples give integer points.
    // -------------------------------------------------------------------------

    type Conic = TractProjectiveVariety<3, 2, FieldTract<i128>>;

    fn conic_relation() -> LinearCombination<[usize; 2], i64> {
        vec![([2, 2], 1i64), ([0, 0], -1i64), ([1, 1], -1i64)]
            .try_into()
            .expect("valid conic relation")
    }

    fn conic_point(x: i128, y: i128, z: i128) -> Result<Conic, TractProjectiveVarietyError> {
        Conic::new([f(x), f(y), f(z)], vec![conic_relation()], false)
    }

    #[test]
    fn pythagorean_triples_on_conic() {
        assert!(conic_point(3, 4, 5).is_ok());
        assert!(conic_point(5, 12, 13).is_ok());
        assert!(conic_point(8, 15, 17).is_ok());
    }

    #[test]
    fn non_pythagorean_off_conic() {
        // 1² + 1² = 2 ≠ 1 = 1².
        assert!(conic_point(1, 1, 1).is_err());
        // 1² + 2² = 5 ≠ 9 = 3².
        assert!(conic_point(1, 2, 3).is_err());
    }

    proptest! {
        // Scaling a Pythagorean triple by any nonzero k stays on the conic.
        #[test]
        fn pythagorean_triple_scaled(k in -200i64..=200i64) {
            prop_assume!(k != 0);
            let k = k as i128;
            prop_assert!(conic_point(3 * k, 4 * k, 5 * k).is_ok());
        }
    }

    // -------------------------------------------------------------------------
    // Error cases
    // -------------------------------------------------------------------------

    #[test]
    fn all_absorbing_rejected() {
        let result = Cubic::new([f(0), f(0), f(0)], vec![cuspidal_cubic_relation()], false);
        assert!(result.is_err());
    }

    #[test]
    fn out_of_range_variable_rejected() {
        // Index 3 is out of range for K_VARIABLES = 3.
        let bad: LinearCombination<[usize; 3], i64> = vec![([0, 0, 3], 1i64)]
            .try_into()
            .expect("valid linear combination");
        let result = Cubic::new([f(1), f(1), f(1)], vec![bad], false);
        assert!(result.is_err());
    }

    #[test]
    fn skip_valid_bypasses_relation_check() {
        // [1:2:1] is not on the cuspidal cubic, but skip_valid lets it through.
        let result = Cubic::new([f(1), f(2), f(1)], vec![cuspidal_cubic_relation()], true);
        assert!(result.is_ok());
    }

    // =========================================================================
    // AbsoluteValueTract<i64>
    //
    // A ℂ-point [x:y:z] maps to [|x|:|y|:|z|] ∈ AbsoluteValueTract.
    // The null-set condition is the polygon inequality 2·max ≤ Σ rᵢ rather than
    // field zero, so the variety condition is coarser than the field condition.
    //
    // Since neg_one = one for AbsoluteValueTract, signs in the relation do not
    // flip the modulus; the condition for the cubic becomes |x|³ = |y|²|z| (the
    // two-term polygon inequality collapses to equality), and the conic condition
    // becomes the triangle inequality for (x², y², z²).
    // =========================================================================

    fn av(n: i64) -> AbsoluteValueTract<i64> {
        AbsoluteValueTract::modulus(n)
    }

    type CubicAV = TractProjectiveVariety<3, 3, AbsoluteValueTract<i64>>;
    type ConicAV = TractProjectiveVariety<3, 2, AbsoluteValueTract<i64>>;

    fn cubic_av(x: i64, y: i64, z: i64) -> Result<CubicAV, TractProjectiveVarietyError> {
        CubicAV::new(
            [av(x), av(y), av(z)],
            vec![cuspidal_cubic_relation()],
            false,
        )
    }

    fn conic_av(x: i64, y: i64, z: i64) -> Result<ConicAV, TractProjectiveVarietyError> {
        ConicAV::new([av(x), av(y), av(z)], vec![conic_relation()], false)
    }

    // The parametric family [t²:t³:1] satisfies |x|³ = t⁶ = |y|²|z|.
    #[test]
    fn av_cubic_parametric_points_on_curve() {
        assert!(cubic_av(1, 1, 1).is_ok());
        assert!(cubic_av(4, 8, 1).is_ok());
        assert!(cubic_av(9, 27, 1).is_ok());
    }

    // Negative inputs are treated as their absolute value, so [-4:8:1] = [4:8:1].
    #[test]
    fn av_cubic_negative_inputs_same_as_positive() {
        assert!(cubic_av(-4, -8, 1).is_ok());
    }

    // |x|³ = 1, |y|²|z| = 4: not equal → not on curve.
    #[test]
    fn av_cubic_off_curve() {
        assert!(cubic_av(1, 2, 1).is_err());
    }

    proptest! {
        // [t²:t³:1] lies on the AbsoluteValueTract cubic for every positive integer t.
        // Limit t to avoid i64 overflow in t⁶ during monomial multiplication.
        #[test]
        fn av_cubic_parametric_family(t in 1i64..=300) {
            prop_assert!(cubic_av(t * t, t * t * t, 1).is_ok());
        }
    }

    // A Pythagorean triple satisfies x²+y²=z², so (z², x², y²) has 2·z² = z²+x²+y²:
    // the polygon inequality holds with equality (degenerate collinear case).
    #[test]
    fn av_conic_pythagorean_triple_on_curve() {
        assert!(conic_av(3, 4, 5).is_ok());
        assert!(conic_av(5, 12, 13).is_ok());
    }

    // [1:1:1]: the three side lengths (1,1,1) form an equilateral triangle — null.
    // This point is NOT on the field conic (1+1 ≠ 1), illustrating the coarser condition.
    #[test]
    fn av_conic_equilateral_on_curve_but_not_field_conic() {
        assert!(conic_av(1, 1, 1).is_ok());
    }

    // z² = 100 > x²+y² = 2: violates the triangle inequality → not on curve.
    #[test]
    fn av_conic_triangle_violated_off_curve() {
        assert!(conic_av(1, 1, 10).is_err());
    }

    // =========================================================================
    // Krasner
    //
    // A ℂ-point [x:y:z] maps to its support via Krasner::final_morphism:
    // nonzero coordinate → Krasner::one(), zero coordinate → absorbing.
    // The null-set condition is: count of nonzero monomials ≠ 1 (i.e., zero or ≥ 2).
    //
    // For the cubic [0,0,0]:+1, [1,1,2]:-1:
    //   monomial x³ is nonzero iff x ≠ 0; monomial y²z is nonzero iff y ≠ 0 and z ≠ 0.
    // For the conic [2,2]:+1, [0,0]:-1, [1,1]:-1:
    //   each monomial is nonzero iff its variable is nonzero.
    // =========================================================================

    fn k(nonzero: bool) -> Krasner {
        if nonzero {
            Krasner::one()
        } else {
            Krasner::absorbing_element()
        }
    }

    type CubicK = TractProjectiveVariety<3, 3, Krasner>;
    type ConicK = TractProjectiveVariety<3, 2, Krasner>;

    fn cubic_k(x: bool, y: bool, z: bool) -> Result<CubicK, TractProjectiveVarietyError> {
        CubicK::new([k(x), k(y), k(z)], vec![cuspidal_cubic_relation()], false)
    }

    fn conic_k(x: bool, y: bool, z: bool) -> Result<ConicK, TractProjectiveVarietyError> {
        ConicK::new([k(x), k(y), k(z)], vec![conic_relation()], false)
    }

    // Both monomials nonzero (count = 2) → null.
    #[test]
    fn krasner_cubic_all_nonzero_on_curve() {
        assert!(cubic_k(true, true, true).is_ok());
    }

    // x=0, y=0: x³ = 0, y²z = 0 (count = 0) → null.  This is the cusp.
    #[test]
    fn krasner_cubic_cusp_support_on_curve() {
        assert!(cubic_k(false, false, true).is_ok());
    }

    // x=0, y ≠ 0 at infinity [0:1:0]: x³=0, y²z = y·y·0 = 0 (count = 0) → null.
    #[test]
    fn krasner_cubic_point_at_infinity_on_curve() {
        assert!(cubic_k(false, true, false).is_ok());
    }

    // x=0, y ≠ 0, z ≠ 0: x³=0, y²z nonzero (count = 1) → not null.
    #[test]
    fn krasner_cubic_only_y2z_nonzero_off_curve() {
        assert!(cubic_k(false, true, true).is_err());
    }

    // x ≠ 0, y=0, z ≠ 0: x³ nonzero, y²z=0 (count = 1) → not null.
    #[test]
    fn krasner_cubic_only_x3_nonzero_off_curve() {
        assert!(cubic_k(true, false, true).is_err());
    }

    // All three monomials nonzero (count = 3) → null.
    #[test]
    fn krasner_conic_all_nonzero_on_curve() {
        assert!(conic_k(true, true, true).is_ok());
    }

    // x ≠ 0, y=0, z ≠ 0: z² and x² nonzero, y²=0 (count = 2) → null.
    #[test]
    fn krasner_conic_two_nonzero_on_curve() {
        assert!(conic_k(true, false, true).is_ok());
    }

    // Only x² nonzero (count = 1) → not null.
    #[test]
    fn krasner_conic_only_one_monomial_nonzero_off_curve() {
        assert!(conic_k(true, false, false).is_err());
    }
}
