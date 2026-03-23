use std::{
    collections::HashMap,
    ops::{Div, DivAssign, Mul, MulAssign, Neg},
};

use num::One;

use crate::{baker_tract::BakerBowlerTract, linear_combination::LinearCombination};

/// Bounds required for a type `S` to serve as the phase group element.
///
/// `S` represents an element of U(1): the multiplicative group of unit complex numbers.
/// The group operation is multiplication; `Div` is the group inverse on the right.
///
/// `Neg` is the **antipodal map**: `-s` is the element at angle `s.angle() + π`.
/// Implementors must ensure `(-S::one()).angle() ≡ π` (mod 2π).
/// This is used to define `neg_one` for the tract.
///
/// `angle()` returns the argument of `s` in `(-π, π]`, embedding `S` into U(1) ⊂ ℝ².
/// It is used for the null set check: 0 lies in the convex hull of the points iff
/// the maximum angular gap between consecutive points (sorted by angle) is at most π.
pub trait PhaseBounds:
    Mul<Self, Output = Self>
    + MulAssign<Self>
    + Div<Self, Output = Self>
    + DivAssign<Self>
    + Neg<Output = Self>
    + One
    + PartialEq
    + Eq
    + std::hash::Hash
    + Clone
    + Sized
{
    /// The argument of this element in `(-π, π]`.
    fn angle(&self) -> f64;
}

/// Phase hyperfield tract over phase group `S`.
///
/// The multiplicative group is U(1) (represented by `S`); the absorbing element is
/// the formal zero outside U(1), represented as `None`.
///
/// Null set: a formal sum is null iff 0 lies in the (closed) convex hull of the
/// corresponding points on the unit circle.  Equivalently, the maximum angular gap
/// between consecutive points sorted by angle is at most π.
#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct PhaseTract<S: PhaseBounds> {
    /// `None` = absorbing element (0 in ℂ, which lies outside U(1)).
    value: Option<S>,
}

impl<S: PhaseBounds> PhaseTract<S> {
    pub fn from_phase(s: S) -> Self {
        Self { value: Some(s) }
    }

    #[must_use]
    pub fn absorbing() -> Self {
        Self { value: None }
    }
}

impl<S: PhaseBounds> MulAssign<Self> for PhaseTract<S> {
    fn mul_assign(&mut self, rhs: Self) {
        match (&mut self.value, rhs.value) {
            (Some(a), Some(b)) => *a *= b,
            _ => self.value = None,
        }
    }
}

impl<S: PhaseBounds> Mul<Self> for PhaseTract<S> {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<S: PhaseBounds> DivAssign<Self> for PhaseTract<S> {
    fn div_assign(&mut self, rhs: Self) {
        assert!(rhs.value.is_some(), "Division by absorbing element");
        if let (Some(a), Some(b)) = (&mut self.value, rhs.value) {
            *a /= b;
        }
        // absorbing / non-absorbing = absorbing: `self.value` stays `None`
    }
}

impl<S: PhaseBounds> Div<Self> for PhaseTract<S> {
    type Output = Self;

    fn div(mut self, rhs: Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<S: PhaseBounds> One for PhaseTract<S> {
    fn one() -> Self {
        Self {
            value: Some(S::one()),
        }
    }
}

impl<S: PhaseBounds> BakerBowlerTract for PhaseTract<S> {
    type NullSet = LinearCombination<Self, usize>;

    fn absorbing_element() -> Self {
        Self { value: None }
    }

    fn neg_one() -> Self {
        // The antipodal element of 1 ∈ U(1), at angle π.
        // `Neg` on `S` is the antipodal map, so `-S::one()` has angle π.
        Self {
            value: Some(-S::one()),
        }
    }

    fn in_null_set(element: Self::NullSet) -> bool {
        let map: HashMap<Self, usize> = element.into();

        // Absorbing elements lie outside U(1) and do not correspond to points on the circle.
        let mut angles: Vec<f64> = map
            .keys()
            .filter_map(|t| t.value.as_ref().map(|s| to_nonneg_angle(s.angle())))
            .collect();

        if angles.is_empty() {
            return true;
        }

        if angles.len() == 1 {
            // A single unit vector is never 0, so 0 ∉ conv({p}).
            return false;
        }

        angles.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        // Maximum angular gap between consecutive sorted points, including the wrap-around gap.
        let wrap_gap = angles[0] + 2.0 * std::f64::consts::PI - angles.last().expect("non-empty");
        let max_gap = angles
            .windows(2)
            .map(|w| w[1] - w[0])
            .chain(std::iter::once(wrap_gap))
            .fold(0.0f64, f64::max);

        // 0 ∈ conv iff no open half-plane contains all points,
        // i.e., the maximum angular gap is ≤ π.
        // The tolerance of 1e-9 absorbs floating-point rounding in `angle()` implementations
        // that compute via degree-to-radian conversion (error ≲ a few ULPs of π ≈ 1e-15).
        //
        // **Error analysis**: the tolerance 1e-9 is chosen to be much larger than any
        // floating-point rounding error in `angle()` (≲ a few ULPs ≈ 1e-15), so a
        // configuration that is truly null (max gap ≤ π) will never be computed as
        // > π + 1e-9 and will never be incorrectly rejected. False negatives are not
        // possible at this tolerance.
        //
        // False positives — reporting null when the max gap is in (π, π + 1e-9] — are
        // possible. Geometrically: two consecutive sorted points span an arc just over π,
        // all other points (however many, at whatever angles) lie in the complementary arc,
        // and we incorrectly report the configuration as null. False positives are the
        // tolerable error direction: a caller that needs exact membership can apply further
        // checks. False negatives — incorrectly rejecting a truly null configuration — would
        // be unrecoverable, and they do not occur here.
        max_gap <= std::f64::consts::PI + 1e-9
    }
}

/// Normalizes an angle from `(-π, π]` to `[0, 2π)` for gap-sorting.
fn to_nonneg_angle(a: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    ((a % two_pi) + two_pi) % two_pi
}

#[cfg(test)]
mod tests {
    use std::{
        f64::consts::PI,
        ops::{Div, DivAssign, Mul, MulAssign, Neg},
    };

    use num::One;
    use proptest::prelude::*;

    use crate::baker_tract::BakerBowlerTract;

    use super::{PhaseBounds, PhaseTract, to_nonneg_angle};

    // -------------------------------------------------------------------------
    // `UnitAngleDeg`: exact integer-degree arithmetic (mod 360).
    // Used for group-law propertests where f64 rounding would cause spurious failures.
    // -------------------------------------------------------------------------

    #[derive(PartialEq, Eq, Hash, Clone, Debug)]
    struct UnitAngleDeg(u16); // angle in whole degrees, 0..360

    impl UnitAngleDeg {
        fn new(deg: i64) -> Self {
            Self(((deg % 360 + 360) % 360) as u16)
        }
    }

    impl Mul for UnitAngleDeg {
        type Output = Self;
        fn mul(self, rhs: Self) -> Self {
            Self::new(self.0 as i64 + rhs.0 as i64)
        }
    }
    impl MulAssign for UnitAngleDeg {
        fn mul_assign(&mut self, rhs: Self) {
            *self = self.clone() * rhs;
        }
    }
    impl Div for UnitAngleDeg {
        type Output = Self;
        fn div(self, rhs: Self) -> Self {
            Self::new(self.0 as i64 - rhs.0 as i64)
        }
    }
    impl DivAssign for UnitAngleDeg {
        fn div_assign(&mut self, rhs: Self) {
            *self = self.clone() / rhs;
        }
    }
    impl Neg for UnitAngleDeg {
        type Output = Self;
        fn neg(self) -> Self {
            Self::new(self.0 as i64 + 180)
        }
    }
    impl One for UnitAngleDeg {
        fn one() -> Self {
            Self(0)
        }
    }
    impl PhaseBounds for UnitAngleDeg {
        fn angle(&self) -> f64 {
            (self.0 as f64).to_radians()
        }
    }

    type PhaseInt = PhaseTract<UnitAngleDeg>;

    fn deg(d: u16) -> PhaseInt {
        PhaseInt::from_phase(UnitAngleDeg(d % 360))
    }

    fn phase_int_strat() -> impl Strategy<Value = PhaseInt> {
        prop_oneof![Just(PhaseInt::absorbing()), (0u16..360).prop_map(deg),]
    }

    fn finite_int_strat() -> impl Strategy<Value = PhaseInt> {
        (0u16..360).prop_map(deg)
    }

    // -------------------------------------------------------------------------
    // Group laws (exact arithmetic via UnitAngleDeg)
    // -------------------------------------------------------------------------

    proptest! {
        #[test]
        fn mul_commutative(a in phase_int_strat(), b in phase_int_strat()) {
            prop_assert_eq!(a.clone() * b.clone(), b * a);
        }

        #[test]
        fn mul_associative(a in phase_int_strat(), b in phase_int_strat(), c in phase_int_strat()) {
            prop_assert_eq!((a.clone() * b.clone()) * c.clone(), a * (b * c));
        }

        #[test]
        fn one_is_identity(a in phase_int_strat()) {
            prop_assert_eq!(a.clone() * PhaseInt::one(), a.clone());
            prop_assert_eq!(PhaseInt::one() * a.clone(), a);
        }

        #[test]
        fn div_inverts_mul(a in finite_int_strat(), b in finite_int_strat()) {
            prop_assert_eq!((a.clone() * b.clone()) / b, a);
        }

        #[test]
        fn absorbing_annihilates(a in phase_int_strat()) {
            prop_assert_eq!(a.clone() * PhaseInt::absorbing(), PhaseInt::absorbing());
            prop_assert_eq!(PhaseInt::absorbing() * a, PhaseInt::absorbing());
        }
    }

    // -------------------------------------------------------------------------
    // neg_one
    // -------------------------------------------------------------------------

    #[test]
    fn neg_one_has_angle_pi() {
        let n = PhaseInt::neg_one();
        let angle = n
            .value
            .as_ref()
            .expect("neg_one is always non-absorbing")
            .angle();
        assert!((angle - PI).abs() < 1e-10);
    }

    // -------------------------------------------------------------------------
    // Null set (fixed cases)
    // -------------------------------------------------------------------------

    #[test]
    fn single_element_not_null() {
        let lc = vec![(deg(45), 1usize)]
            .try_into()
            .expect("valid linear combination");
        assert!(!PhaseInt::in_null_set(lc));
    }

    proptest! {
        // Any pair of antipodal points {a°, (a+180)°} has max gap = π → null.
        #[test]
        fn antipodal_pair_is_null(a in 0u16..180) {
            let lc = vec![(deg(a), 1usize), (deg(a + 180), 1usize)]
                .try_into()
                .expect("valid linear combination");
            prop_assert!(PhaseInt::in_null_set(lc));
        }
    }

    #[test]
    fn nearby_pair_not_null() {
        // 10° and 20°: max gap = 340° > π → not null.
        let lc = vec![(deg(10), 1usize), (deg(20), 1usize)]
            .try_into()
            .expect("valid linear combination");
        assert!(!PhaseInt::in_null_set(lc));
    }

    #[test]
    fn three_balanced_is_null() {
        // 0°, 120°, 240°: max gap = 120° < π → null.
        let lc = vec![(deg(0), 1usize), (deg(120), 1usize), (deg(240), 1usize)]
            .try_into()
            .expect("valid linear combination");
        assert!(PhaseInt::in_null_set(lc));
    }

    proptest! {
        // Two *distinct* points are null iff their angular separation is exactly π.
        // Equal points merge into one entry in LinearCombination; that case is excluded.
        #[test]
        fn two_distinct_points_null_iff_gap_at_least_pi(
            a in 0u16..360,
            offset in 1u16..360,
        ) {
            let b = (a + offset) % 360;
            let lc = vec![(deg(a), 1usize), (deg(b), 1usize)].try_into().expect("valid linear combination");
            let gap = {
                let mut angles = [
                    to_nonneg_angle(a as f64 * PI / 180.0),
                    to_nonneg_angle(b as f64 * PI / 180.0),
                ];
                angles.sort_by(|x, y| x.partial_cmp(y).expect("test angles are finite"));
                let d = angles[1] - angles[0];
                d.max(2.0 * PI - d)
            };
            prop_assert_eq!(PhaseInt::in_null_set(lc), gap <= PI + 1e-9);
        }

        // Absorbing elements do not affect the null set verdict.
        #[test]
        fn absorbing_ignored_in_null_set(a in 0u16..360, n in 1usize..5) {
            let mut v = vec![(deg(a), 1usize)];
            v.extend((0..n).map(|_| (PhaseInt::absorbing(), 1)));
            let lc = v.try_into().expect("valid linear combination");
            prop_assert!(!PhaseInt::in_null_set(lc));
        }
    }
}
