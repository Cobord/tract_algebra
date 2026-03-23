use std::ops::{Div, DivAssign, Mul, MulAssign};

use num::One;

use crate::{baker_tract::BakerBowlerTract, linear_combination::LinearCombination};

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct Krasner(bool);

impl MulAssign<Self> for Krasner {
    fn mul_assign(&mut self, rhs: Self) {
        if self.0 {
            *self = rhs;
        }
    }
}

impl Mul<Self> for Krasner {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl DivAssign<Self> for Krasner {
    fn div_assign(&mut self, rhs: Self) {
        assert!(rhs.0, "Division by zero");
    }
}

impl Div<Self> for Krasner {
    type Output = Self;

    fn div(mut self, rhs: Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl One for Krasner {
    fn one() -> Self {
        Self(true)
    }
}

impl BakerBowlerTract for Krasner {
    type NullSet = LinearCombination<Krasner, usize>;

    fn absorbing_element() -> Self {
        Self(false)
    }

    fn neg_one() -> Self {
        Self(true)
    }

    fn in_null_set(element: Self::NullSet) -> bool {
        let mut count_one = 0;
        for (key, value) in element.map {
            if key.0 {
                count_one += value;
            }
        }
        count_one == 0 || count_one > 1
    }
}

impl Krasner {
    #[allow(clippy::needless_pass_by_value)]
    pub fn final_morphism<T: BakerBowlerTract>(t: T) -> Self {
        if t.is_absorbing_elt() {
            Self::absorbing_element()
        } else {
            Self::one()
        }
    }

    pub fn final_morphism_ref<T: BakerBowlerTract>(t: &T) -> Self {
        if t.is_absorbing_elt() {
            Self::absorbing_element()
        } else {
            Self::one()
        }
    }
}
