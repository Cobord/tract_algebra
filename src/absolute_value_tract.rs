use std::{
    collections::HashMap,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg},
};

use num::{One, Zero};

use crate::{baker_tract::BakerBowlerTract, linear_combination::LinearCombination};

/// Bounds for the carrier `R` of the absolute value tract.
///
/// `R` represents a positive real number (the modulus of a nonzero complex number).
/// - Multiplication is the group operation (product of moduli).
/// - Addition and `PartialOrd` are used only in `in_null_set`.
/// - `Zero` provides the additive identity for accumulating the total.
/// - `abs` returns the non-negative representative; automatically provided for
///   types that implement `Neg<Output = Self> + Ord`.
pub trait AbsoluteValueBounds:
    Mul<Self, Output = Self>
    + MulAssign<Self>
    + Div<Self, Output = Self>
    + DivAssign<Self>
    + Add<Self, Output = Self>
    + AddAssign<Self>
    + One
    + Zero
    + PartialOrd
    + std::hash::Hash
    + Eq
    + Clone
    + Sized
{
    #[must_use = "Replaced by the absolute value so it is the representative of the C^*/U(1) coset"]
    fn abs(self) -> Self;
}

impl<T> AbsoluteValueBounds for T
where
    T: Mul<T, Output = T>
        + MulAssign<T>
        + Div<T, Output = T>
        + DivAssign<T>
        + Add<T, Output = T>
        + AddAssign<T>
        + One
        + Zero
        + PartialOrd
        + std::hash::Hash
        + Eq
        + Clone
        + Neg<Output = T>
        + Ord,
{
    fn abs(self) -> T {
        if self < T::zero() { -self } else { self }
    }
}

/// Tract corresponding to `FieldTract(ℂ) / U(1)` — the quotient of the complex field
/// tract by the unit circle subgroup.
///
/// Elements are cosets of U(1) in ℂ*, each identified with a positive real (the modulus).
/// The multiplicative group is (ℝ>0, ×) under scaling.
/// The absorbing element represents 0 ∈ ℂ, which lies outside ℂ*.
///
/// **Null set**: a formal sum `{r₁, …, rₙ}` is null iff there exist directions
/// `u₁, …, uₙ ∈ U(1)` with `Σ uᵢrᵢ = 0` in ℂ, i.e. the `rᵢ` can be realised as side
/// lengths of a closed possibly degenerate polygon in ℝ².  This is equivalent to `2 · max(rᵢ) ≤ Σ rᵢ`
/// (weak inequality, including degenerate collinear cases).  The weak form is required
/// for the tract axioms: the strict version has no 2-element null sets, so `neg_one`
/// could not exist.
/// So calling it triangular or polygon tract can be misleading.
///
/// Compare: `PhaseTract` = `FieldTract(ℂ) / ℝ>0` (quotient by positive reals, keeping angle).
///
/// `in_null_set` panics when `PartialOrd` returns `None` (e.g. NaN).
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct AbsoluteValueTract<R: AbsoluteValueBounds> {
    /// `None` = absorbing element (zero-length, outside the multiplicative group).
    value: Option<R>,
}

impl<R: AbsoluteValueBounds> AbsoluteValueTract<R> {
    pub fn modulus(r: R) -> Self {
        Self {
            value: Some(r.abs()),
        }
    }

    #[must_use]
    pub fn absorbing() -> Self {
        Self { value: None }
    }
}

impl<R: AbsoluteValueBounds> MulAssign<Self> for AbsoluteValueTract<R> {
    fn mul_assign(&mut self, rhs: Self) {
        match (&mut self.value, rhs.value) {
            (Some(a), Some(b)) => *a *= b,
            _ => self.value = None,
        }
    }
}

impl<R: AbsoluteValueBounds> Mul<Self> for AbsoluteValueTract<R> {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<R: AbsoluteValueBounds> DivAssign<Self> for AbsoluteValueTract<R> {
    fn div_assign(&mut self, rhs: Self) {
        assert!(rhs.value.is_some(), "Division by absorbing element");
        if let (Some(a), Some(b)) = (&mut self.value, rhs.value) {
            *a /= b;
        }
        // absorbing / non-absorbing = absorbing: `value` stays `None`
    }
}

impl<R: AbsoluteValueBounds> Div<Self> for AbsoluteValueTract<R> {
    type Output = Self;

    fn div(mut self, rhs: Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<R: AbsoluteValueBounds> One for AbsoluteValueTract<R> {
    fn one() -> Self {
        Self {
            value: Some(R::one()),
        }
    }
}

impl<R: AbsoluteValueBounds> BakerBowlerTract for AbsoluteValueTract<R> {
    type NullSet = LinearCombination<Self, usize>;

    fn absorbing_element() -> Self {
        Self { value: None }
    }

    fn neg_one() -> Self {
        // {1, g} is null iff 2·max(1,g) ≤ 1+g.
        // For g = 1: 2 ≤ 2, true. For g < 1: max = 1, so 2 > 1+g, not null.
        // For g > 1: max = g, so 2g > 1+g, not null. So the only solution is g = 1.
        Self::one()
    }

    fn in_null_set(element: Self::NullSet) -> bool {
        let map: HashMap<Self, usize> = element.into();

        // Absorbing elements (zero-length) do not contribute.
        // In the degenerate polygon picture, you can put two
        // vertices on top of each other.
        let finite: Vec<(&R, usize)> = map
            .iter()
            .filter_map(|(t, &n)| t.value.as_ref().map(|r| (r, n)))
            .collect();

        if finite.is_empty() {
            return true;
        }

        // Accumulate total = Σ nᵢ·rᵢ and track the maximum side length.
        let mut total = R::zero();
        let mut max: &R = finite[0].0;
        for (r, n) in &finite {
            for _ in 0..*n {
                total += (*r).clone();
            }
            match (*r).partial_cmp(max) {
                Some(std::cmp::Ordering::Greater) => max = r,
                Some(_) => {}
                None => panic!("in_null_set: incomparable side lengths;"),
            }
        }

        // Possibly Degenerate Polygon inequality: 2·max ≤ total.
        let two_max = max.clone() + max.clone();
        match two_max.partial_cmp(&total) {
            Some(o) => !o.is_gt(),
            None => panic!("in_null_set: incomparable side lengths;"),
        }
    }
}

#[cfg(test)]
mod tests {
    use num::One;
    use ordered_float::OrderedFloat;
    use proptest::prelude::*;

    use crate::baker_tract::BakerBowlerTract;

    use super::AbsoluteValueTract;

    // i64 (positive range) for exact arithmetic in group-law tests.
    // OrderedFloat<f64> for null-set tests that only need ordering.
    type AbsValI = AbsoluteValueTract<i64>;
    type AbsValF = AbsoluteValueTract<OrderedFloat<f64>>;

    fn modulus_i(n: i64) -> AbsValI {
        AbsValI::modulus(n)
    }

    fn modulus_f(x: f64) -> AbsValF {
        AbsValF::modulus(OrderedFloat(x))
    }

    fn abs_val_i_finite() -> impl Strategy<Value = AbsValI> {
        (1i64..=100).prop_map(modulus_i)
    }

    fn abs_val_i_any() -> impl Strategy<Value = AbsValI> {
        prop_oneof![Just(AbsValI::absorbing()), abs_val_i_finite()]
    }

    // --- constant properties ---

    #[test]
    fn neg_one_equals_one() {
        assert_eq!(AbsValI::neg_one(), AbsValI::one());
        assert_eq!(AbsValF::neg_one(), AbsValF::one());
    }

    // --- group laws (exact i64 arithmetic) ---

    proptest! {
        #[test]
        fn mul_commutative(a in abs_val_i_any(), b in abs_val_i_any()) {
            prop_assert_eq!(a.clone() * b.clone(), b * a);
        }

        #[test]
        fn mul_associative(a in abs_val_i_any(), b in abs_val_i_any(), c in abs_val_i_any()) {
            prop_assert_eq!((a.clone() * b.clone()) * c.clone(), a * (b * c));
        }

        #[test]
        fn one_is_identity(a in abs_val_i_any()) {
            prop_assert_eq!(a.clone() * AbsValI::one(), a.clone());
            prop_assert_eq!(AbsValI::one() * a.clone(), a);
        }

        #[test]
        fn div_inverts_mul(a in abs_val_i_finite(), b in abs_val_i_finite()) {
            prop_assert_eq!((a.clone() * b.clone()) / b, a);
        }

        #[test]
        fn absorbing_annihilates(a in abs_val_i_any()) {
            prop_assert_eq!(a.clone() * AbsValI::absorbing(), AbsValI::absorbing());
            prop_assert_eq!(AbsValI::absorbing() * a, AbsValI::absorbing());
        }
    }

    // --- null set (OrderedFloat<f64> for ordering) ---

    #[test]
    fn equilateral_triangle_is_null() {
        let lc = vec![(modulus_f(1.0), 3usize)]
            .try_into()
            .expect("valid linear combination");
        assert!(AbsValF::in_null_set(lc));
    }

    #[test]
    fn degenerate_two_equal_sides_is_null() {
        // Two equal sides: 2·max = total → null (degenerate polygon).
        let lc = vec![(modulus_f(3.0), 2usize)]
            .try_into()
            .expect("valid linear combination");
        assert!(AbsValF::in_null_set(lc));
    }

    #[test]
    fn single_side_not_null() {
        let lc = vec![(modulus_f(5.0), 1usize)]
            .try_into()
            .expect("valid linear combination");
        assert!(!AbsValF::in_null_set(lc));
    }

    #[test]
    fn triangle_inequality_violation_not_null() {
        // 5 > 2 + 1: fails possibly degenerate polygon inequality.
        let lc = vec![
            (modulus_f(5.0), 1usize),
            (modulus_f(2.0), 1usize),
            (modulus_f(1.0), 1usize),
        ]
        .try_into()
        .expect("valid linear combination");
        assert!(!AbsValF::in_null_set(lc));
    }

    proptest! {
        // Any n ≥ 2 equal positive sides form a valid (possibly degenerate) polygon.
        #[test]
        fn equal_sides_always_null(x in 1e-3f64..1e3f64, n in 2usize..10) {
            let lc = vec![(modulus_f(x), n)].try_into().expect("valid linear combination");
            prop_assert!(AbsValF::in_null_set(lc));
        }

        // If the longest side exceeds the sum of all others, not null.
        #[test]
        fn strict_violation_not_null(
            rest in prop::collection::vec(1e-3f64..10.0f64, 2..5),
            excess in 1e-3f64..1.0f64,
        ) {
            let rest_sum: f64 = rest.iter().sum();
            let big = rest_sum + excess;
            let mut v: Vec<(AbsValF, usize)> = rest.iter().map(|&r| (modulus_f(r), 1)).collect();
            v.push((modulus_f(big), 1));
            let lc = v.try_into().expect("valid linear combination");
            prop_assert!(!AbsValF::in_null_set(lc));
        }

        // Absorbing elements do not affect the verdict.
        #[test]
        fn absorbing_ignored(x in 1e-3f64..1e3f64, n_abs in 1usize..5) {
            let mut v = vec![(modulus_f(x), 1usize)];
            v.extend((0..n_abs).map(|_| (AbsValF::absorbing(), 1)));
            let lc = v.try_into().expect("valid linear combination");
            prop_assert!(!AbsValF::in_null_set(lc)); // single finite side: never null
        }
    }
}
