use std::ops::{Div, DivAssign, Mul, MulAssign};

use num::One;

use crate::{baker_tract::BakerBowlerTract, linear_combination::LinearCombination};

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub struct RegularPartialField(Option<bool>);

impl MulAssign<Self> for RegularPartialField {
    fn mul_assign(&mut self, rhs: Self) {
        match self.0 {
            Some(true) => {
                *self = rhs;
            }
            Some(false) => match rhs.0 {
                Some(true) => *self = Self(Some(false)),
                Some(false) => *self = Self(Some(true)),
                None => *self = Self(None),
            },
            None => {}
        }
    }
}

impl Mul<Self> for RegularPartialField {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}

impl DivAssign<Self> for RegularPartialField {
    fn div_assign(&mut self, rhs: Self) {
        assert!(rhs.0.is_some(), "Division by zero");
        *self *= rhs;
    }
}

impl Div<Self> for RegularPartialField {
    type Output = Self;

    fn div(mut self, rhs: Self) -> Self::Output {
        self /= rhs;
        self
    }
}

impl One for RegularPartialField {
    fn one() -> Self {
        Self(Some(true))
    }
}

impl BakerBowlerTract for RegularPartialField {
    type NullSet = LinearCombination<RegularPartialField, usize>;

    fn absorbing_element() -> Self {
        Self(None)
    }

    fn neg_one() -> Self {
        Self(Some(false))
    }

    fn in_null_set(element: Self::NullSet) -> bool {
        let mut count_one = 0;
        let mut count_neg_one = 0;
        for (key, value) in element.map {
            match key.0 {
                Some(true) => count_one += value,
                Some(false) => count_neg_one += value,
                None => {}
            }
        }
        count_neg_one == count_one
    }
}

impl RegularPartialField {
    #[must_use = "You have applied the initial morphism from initial object to some target T. How is T utilized?"]
    pub fn initial_morphism<T: BakerBowlerTract>(self) -> T {
        match self.0 {
            Some(true) => T::one(),
            Some(false) => T::neg_one(),
            None => T::absorbing_element(),
        }
    }
}
