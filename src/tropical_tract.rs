use std::cmp::Ordering;
use std::{
    collections::HashMap,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign},
};

use num::{One, Zero};

use crate::{baker_tract::BakerBowlerTract, linear_combination::LinearCombination};

/// Minimal bounds for the carrier set `R` of a tropical semiring.
/// Supports tropical multiplication (`R`-addition) and comparison for the null set,
/// but does not require invertibility (no `Sub`/`SubAssign`).
pub trait TropicalBounds:
    Add<Self, Output = Self>
    + AddAssign<Self>
    + PartialOrd
    + Zero
    + std::hash::Hash
    + Eq
    + Clone
    + Sized
{
}

impl<T> TropicalBounds for T where
    T: Add<T, Output = T> + AddAssign<T> + PartialOrd + Zero + std::hash::Hash + Eq + Clone
{
}

/// Stronger bounds required when `R`-subtraction is available, i.e. when the
/// tropical multiplication is invertible and `TropicalTract` can form a tract
/// (group structure on non-absorbing elements).
pub trait TropicalGroupBounds: TropicalBounds + Sub<Self, Output = Self> + SubAssign<Self> {}

impl<T> TropicalGroupBounds for T where T: TropicalBounds + Sub<T, Output = T> + SubAssign<T> {}

/// Tropical tract over carrier `R` in **log-coordinates** (exponents).
///
/// An element `a ∈ R` represents the physical value `e^{a/ℏ}` in the
/// Maslov dequantization picture.  Because the map `a ↦ e^{a/ℏ}` covers
/// all of ℝ>0, the exponent space is unrestricted: `R` should be a group
/// under addition (e.g. ℤ, ℚ, ℝ).  In particular, **ℝ≥0 under addition
/// is not a valid `R` for `TropicalGroupBounds`**: subtraction can produce
/// negative exponents, exiting the set.  If the physical values live in
/// ℝ>0 and you want to work in linear (non-log) coordinates, that is a
/// different multiplicative tropical structure, not this type.
///
/// `MAX_VS_MIN = true`  → (max, +): absorbing element is -∞.
/// `MAX_VS_MIN = false` → (min, +): absorbing element is +∞.
///
/// Tract multiplication is addition in `R`.
/// The multiplicative identity is `R::zero()` (the exponent of 1 = e^0).
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct TropicalTract<R: TropicalBounds, const MAX_VS_MIN: bool> {
    /// `None` = absorbing element (±∞ depending on convention).
    value: Option<R>,
}

impl<R: TropicalBounds, const MAX_VS_MIN: bool> TropicalTract<R, MAX_VS_MIN> {
    pub fn finite(r: R) -> Self {
        Self { value: Some(r) }
    }

    #[must_use = "the absorbing element is the formal ±∞ adjoined to R: it annihilates tropical multiplication and is the identity for tropical addition"]
    pub fn absorbing() -> Self {
        Self { value: None }
    }
}

// --- Add/AddAssign/Zero: tropical semiring addition (max or min), only need TropicalBounds ---

impl<R: TropicalBounds, const MAX_VS_MIN: bool> AddAssign<Self> for TropicalTract<R, MAX_VS_MIN> {
    fn add_assign(&mut self, rhs: Self) {
        self.value = match (self.value.take(), rhs.value) {
            // Absorbing element is the additive identity (-∞ for max, +∞ for min).
            (None, other) | (other, None) => other,
            (Some(a), Some(b)) => Some(match a.partial_cmp(&b) {
                Some(Ordering::Greater) => {
                    if MAX_VS_MIN {
                        a
                    } else {
                        b
                    }
                }
                Some(Ordering::Less) => {
                    if MAX_VS_MIN {
                        b
                    } else {
                        a
                    }
                }
                // Equal or incomparable: keep a.
                _ => a,
            }),
        };
    }
}

impl<R: TropicalBounds, const MAX_VS_MIN: bool> Add<Self> for TropicalTract<R, MAX_VS_MIN> {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<R: TropicalBounds, const MAX_VS_MIN: bool> Zero for TropicalTract<R, MAX_VS_MIN> {
    fn zero() -> Self {
        Self::absorbing()
    }

    fn is_zero(&self) -> bool {
        self.value.is_none()
    }
}

// --- Mul/MulAssign/One: only need TropicalBounds ---

#[allow(clippy::suspicious_op_assign_impl)]
impl<R: TropicalBounds, const MAX_VS_MIN: bool> MulAssign<Self> for TropicalTract<R, MAX_VS_MIN> {
    fn mul_assign(&mut self, rhs: Self) {
        match (&mut self.value, rhs.value) {
            (Some(a), Some(b)) => *a += b,
            _ => self.value = None,
        }
    }
}

impl<R: TropicalBounds, const MAX_VS_MIN: bool> Mul<Self> for TropicalTract<R, MAX_VS_MIN> {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<R: TropicalBounds, const MAX_VS_MIN: bool> One for TropicalTract<R, MAX_VS_MIN> {
    fn one() -> Self {
        Self {
            value: Some(R::zero()),
        }
    }
}

// --- Div/DivAssign/BakerBowlerTract: need TropicalGroupBounds ---

impl<R: TropicalGroupBounds, const MAX_VS_MIN: bool> DivAssign<Self>
    for TropicalTract<R, MAX_VS_MIN>
{
    fn div_assign(&mut self, rhs: Self) {
        assert!(rhs.value.is_some(), "Division by absorbing element");
        if let (Some(a), Some(b)) = (&mut self.value, rhs.value) {
            *a -= b;
        }
        // absorbing / finite = absorbing: self.value stays None
    }
}

impl<R: TropicalGroupBounds, const MAX_VS_MIN: bool> Div<Self> for TropicalTract<R, MAX_VS_MIN> {
    type Output = Self;

    fn div(mut self, rhs: Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<R: TropicalGroupBounds, const MAX_VS_MIN: bool> BakerBowlerTract
    for TropicalTract<R, MAX_VS_MIN>
{
    type NullSet = LinearCombination<Self, usize>;

    fn absorbing_element() -> Self {
        Self::absorbing()
    }

    fn neg_one() -> Self {
        // The null set requires the extremal value to appear with multiplicity >= 2,
        // so {one, x} in N_T forces x = one, hence neg_one = one.
        Self::one()
    }

    fn in_null_set(element: Self::NullSet) -> bool {
        let map: HashMap<Self, usize> = element.into();

        // Absorbing elements do not contribute; collect finite entries.
        let finite: Vec<(&R, usize)> = map
            .iter()
            .filter_map(|(t, &n)| t.value.as_ref().map(|r| (r, n)))
            .collect();

        if finite.is_empty() {
            return true;
        }

        // Find the extremal (max or min) R-value.
        let extremal = finite.iter().fold(None::<&R>, |best, (r, _)| match best {
            None => Some(r),
            Some(b) => {
                let r_is_better = if MAX_VS_MIN {
                    (*r).partial_cmp(b).is_some_and(std::cmp::Ordering::is_gt)
                } else {
                    (*r).partial_cmp(b).is_some_and(std::cmp::Ordering::is_lt)
                };
                Some(if r_is_better { r } else { b })
            }
        });

        let Some(extremal) = extremal else {
            return true;
        };

        // Null iff the extremal value has total multiplicity >= 2.
        let count: usize = finite
            .iter()
            .filter(|(r, _)| *r == extremal)
            .map(|(_, n)| n)
            .sum();

        count >= 2
    }
}

#[cfg(test)]
mod tests {
    use num::{One, Zero};
    use ordered_float::OrderedFloat;
    use proptest::prelude::*;

    use crate::baker_tract::BakerBowlerTract;

    use super::TropicalTract;

    // f64 with OrderedFloat for ordering-based tests (null set, tropical addition).
    type MaxPlusF = TropicalTract<OrderedFloat<f64>, true>;
    type MinPlusF = TropicalTract<OrderedFloat<f64>, false>;

    // i64 for arithmetic-identity tests (mul associativity, div/mul roundtrip):
    // f64 addition is not exact, so those laws only hold for exact types.
    type MaxPlusI = TropicalTract<i64, true>;

    fn of(x: f64) -> OrderedFloat<f64> {
        OrderedFloat(x)
    }

    fn max_f_finite() -> impl Strategy<Value = MaxPlusF> {
        (-1e3f64..1e3f64).prop_map(|x| MaxPlusF::finite(of(x)))
    }

    fn max_f_any() -> impl Strategy<Value = MaxPlusF> {
        prop_oneof![Just(MaxPlusF::zero()), max_f_finite()]
    }

    fn min_f_finite() -> impl Strategy<Value = MinPlusF> {
        (-1e3f64..1e3f64).prop_map(|x| MinPlusF::finite(of(x)))
    }

    fn min_f_any() -> impl Strategy<Value = MinPlusF> {
        prop_oneof![Just(MinPlusF::zero()), min_f_finite()]
    }

    fn max_i_finite() -> impl Strategy<Value = MaxPlusI> {
        // Divide range by 3 so that a + b + c never overflows in associativity tests.
        (i64::MIN / 3..i64::MAX / 3).prop_map(MaxPlusI::finite)
    }

    fn max_i_any() -> impl Strategy<Value = MaxPlusI> {
        prop_oneof![Just(MaxPlusI::zero()), max_i_finite()]
    }

    // --- constant properties (no random input needed) ---

    #[test]
    fn neg_one_equals_one() {
        assert_eq!(MaxPlusF::neg_one(), MaxPlusF::one());
        assert_eq!(MinPlusF::neg_one(), MinPlusF::one());
    }

    // --- semiring laws ---
    // Ordering-only laws (max/min is exact): use OrderedFloat<f64>.
    // Arithmetic-identity laws (addition in R): use i64 to avoid f64 rounding.

    proptest! {
        #[test]
        fn zero_additive_identity(a in max_f_any()) {
            prop_assert_eq!(a.clone() + MaxPlusF::zero(), a.clone());
            prop_assert_eq!(MaxPlusF::zero() + a.clone(), a);
        }

        #[test]
        fn one_multiplicative_identity(a in max_i_any()) {
            prop_assert_eq!(a.clone() * MaxPlusI::one(), a.clone());
            prop_assert_eq!(MaxPlusI::one() * a.clone(), a);
        }

        #[test]
        fn absorbing_annihilates_mul(a in max_i_any()) {
            prop_assert_eq!(a.clone() * MaxPlusI::zero(), MaxPlusI::zero());
            prop_assert_eq!(MaxPlusI::zero() * a, MaxPlusI::zero());
        }

        #[test]
        fn add_commutative(a in max_f_any(), b in max_f_any()) {
            prop_assert_eq!(a.clone() + b.clone(), b + a);
        }

        #[test]
        fn mul_commutative(a in max_i_any(), b in max_i_any()) {
            prop_assert_eq!(a.clone() * b.clone(), b * a);
        }

        #[test]
        fn add_associative(a in max_f_any(), b in max_f_any(), c in max_f_any()) {
            prop_assert_eq!((a.clone() + b.clone()) + c.clone(), a + (b + c));
        }

        // Exact arithmetic: use i64 to avoid f64 rounding making (a+b)+c != a+(b+c).
        #[test]
        fn mul_associative(a in max_i_any(), b in max_i_any(), c in max_i_any()) {
            prop_assert_eq!((a.clone() * b.clone()) * c.clone(), a * (b * c));
        }

        // Tropical addition is idempotent: a + a == a
        #[test]
        fn add_idempotent(a in max_f_any()) {
            prop_assert_eq!(a.clone() + a.clone(), a);
        }

        // Distributivity: a * (b + c) == (a * b) + (a * c)
        #[test]
        fn distributive(a in max_i_any(), b in max_i_any(), c in max_i_any()) {
            prop_assert_eq!(
                a.clone() * (b.clone() + c.clone()),
                (a.clone() * b) + (a * c),
            );
        }

        // Same laws for MinPlus
        #[test]
        fn min_add_commutative(a in min_f_any(), b in min_f_any()) {
            prop_assert_eq!(a.clone() + b.clone(), b + a);
        }

        #[test]
        fn min_distributive(a in min_f_any(), b in min_f_any(), c in min_f_any()) {
            prop_assert_eq!(
                a.clone() * (b.clone() + c.clone()),
                (a.clone() * b) + (a * c),
            );
        }
    }

    // --- tract: division ---

    proptest! {
        // Exact arithmetic: use i64 so (a * b) / b == a holds without rounding.
        #[test]
        fn div_inverts_mul(a in max_i_finite(), b in max_i_finite()) {
            prop_assert_eq!((a.clone() * b.clone()) / b, a);
        }

        #[test]
        fn absorbing_div_finite_stays_absorbing(b in max_i_finite()) {
            prop_assert_eq!(MaxPlusI::zero() / b, MaxPlusI::zero());
        }
    }

    // --- tract: null set ---

    proptest! {
        // Any finite element with multiplicity >= 2 is null.
        #[test]
        fn null_set_repeated_finite(x in -1e3f64..1e3f64, n in 2usize..10) {
            let lc = vec![(MaxPlusF::finite(of(x)), n)].try_into().expect("valid linear combination");
            prop_assert!(MaxPlusF::in_null_set(lc));
        }

        // In max-plus, if hi > lo, the pair {hi:1, lo:1} is not null (max achieved once).
        #[test]
        fn null_set_strict_max_not_null(lo in -1e3f64..1e3f64, offset in 1e-3f64..1e3f64) {
            let hi = lo + offset;
            let lc = vec![
                (MaxPlusF::finite(of(hi)), 1usize),
                (MaxPlusF::finite(of(lo)), 1usize),
            ].try_into().expect("valid linear combination");
            prop_assert!(!MaxPlusF::in_null_set(lc));
        }

        // In min-plus, if lo < hi, the pair {lo:1, hi:1} is not null (min achieved once).
        #[test]
        fn null_set_strict_min_not_null(lo in -1e3f64..1e3f64, offset in 1e-3f64..1e3f64) {
            let hi = lo + offset;
            let lc = vec![
                (MinPlusF::finite(of(lo)), 1usize),
                (MinPlusF::finite(of(hi)), 1usize),
            ].try_into().expect("valid linear combination");
            prop_assert!(!MinPlusF::in_null_set(lc));
        }

        // Max achieved with total multiplicity >= 2 is null, even with lower-valued terms.
        #[test]
        fn null_set_tied_max_is_null(lo in -1e3f64..1e3f64, offset in 1e-3f64..1e3f64, n in 2usize..10) {
            let hi = lo + offset;
            let lc = vec![
                (MaxPlusF::finite(of(hi)), n),
                (MaxPlusF::finite(of(lo)), 1usize),
            ].try_into().expect("valid linear combination");
            prop_assert!(MaxPlusF::in_null_set(lc));
        }

        // Absorbing elements never affect the null set verdict.
        #[test]
        fn null_set_absorbing_ignored(x in -1e3f64..1e3f64, n_absorbing in 1usize..10) {
            // Single finite element with multiplicity 1: not null regardless of absorbing count.
            let mut v = vec![(MaxPlusF::finite(of(x)), 1usize)];
            v.extend((0..n_absorbing).map(|_| (MaxPlusF::zero(), 1usize)));
            let lc = v.try_into().expect("valid linear combination");
            prop_assert!(!MaxPlusF::in_null_set(lc));
        }
    }
}
