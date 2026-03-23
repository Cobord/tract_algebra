use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{One, Zero};

use crate::{
    baker_tract::{BakerBowlerTract, MulGroupBounds},
    linear_combination::LinearCombination,
};

pub trait FieldBounds:
    MulGroupBounds
    + Add<Self, Output = Self>
    + AddAssign<Self>
    + Sub<Self, Output = Self>
    + SubAssign<Self>
    + Neg<Output = Self>
    + Zero
{
}

impl<T> FieldBounds for T where
    T: MulGroupBounds
        + Add<Self, Output = Self>
        + AddAssign<Self>
        + Sub<Self, Output = Self>
        + SubAssign<Self>
        + Neg<Output = Self>
        + Zero
{
}

#[derive(Hash, Clone, PartialEq, Eq, Debug)]
pub struct FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    field_elt: F,
}

impl<F> From<F> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn from(field_elt: F) -> Self {
        Self { field_elt }
    }
}

impl<F> Add<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            field_elt: self.field_elt + rhs.field_elt,
        }
    }
}

impl<F> AddAssign<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn add_assign(&mut self, rhs: Self) {
        self.field_elt += rhs.field_elt;
    }
}

impl<F> Sub<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            field_elt: self.field_elt - rhs.field_elt,
        }
    }
}

impl<F> SubAssign<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.field_elt -= rhs.field_elt;
    }
}

impl<F> Mul<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            field_elt: self.field_elt * rhs.field_elt,
        }
    }
}

impl<F> Mul<usize> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type Output = Self;

    fn mul(mut self, rhs: usize) -> Self::Output {
        if rhs == 0 {
            Self::zero()
        } else if rhs % 2 == 1 {
            let mut answer = self.clone();
            self = self * ((rhs - 1) / 2);
            answer += self.clone();
            answer += self;
            answer
        } else {
            self = self * (rhs / 2);
            self.clone() + self
        }
    }
}

impl<F> MulAssign<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn mul_assign(&mut self, rhs: Self) {
        self.field_elt *= rhs.field_elt;
    }
}

impl<F> Div<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        assert!(!rhs.field_elt.is_zero(), "Division by zero");
        Self {
            field_elt: self.field_elt / rhs.field_elt,
        }
    }
}

impl<F> DivAssign<Self> for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn div_assign(&mut self, rhs: Self) {
        assert!(!rhs.field_elt.is_zero(), "Division by zero");
        self.field_elt /= rhs.field_elt;
    }
}

impl<F> Zero for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn is_zero(&self) -> bool {
        self.field_elt.is_zero()
    }

    fn zero() -> Self {
        Self {
            field_elt: F::zero(),
        }
    }
}

impl<F> One for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    fn one() -> Self {
        Self {
            field_elt: F::one(),
        }
    }
}

impl<F> Neg for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            field_elt: -self.field_elt,
        }
    }
}

impl<F> BakerBowlerTract for FieldTract<F>
where
    F: FieldBounds + PartialEq + Eq + Clone + std::hash::Hash + Sized,
{
    type NullSet = LinearCombination<Self, usize>;

    fn absorbing_element() -> Self {
        Self {
            field_elt: F::zero(),
        }
    }

    fn neg_one() -> Self {
        -Self::one()
    }

    fn in_null_set(element: Self::NullSet) -> bool {
        let sum = element
            .map
            .into_iter()
            .fold(F::zero(), |acc, next| acc + (next.0 * next.1).field_elt);
        sum.is_zero()
    }
}
