use std::{
    collections::HashMap,
    ops::{Add, AddAssign, Mul, MulAssign},
};

use num::Zero;

#[derive(Clone, Debug)]
pub struct LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
{
    pub(crate) map: HashMap<Set, Coeffs>,
}

impl<Set, Coeffs> AddAssign<Self> for LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
{
    fn add_assign(&mut self, rhs: Self) {
        for (key, v) in rhs.map {
            self.map
                .entry(key)
                .and_modify(|old_value| *old_value += v.clone())
                .or_insert(v);
        }
    }
}

impl<Set, Coeffs> Add<Self> for LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
{
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self += rhs;
        self
    }
}

impl<Set, Coeffs> Zero for LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
{
    fn zero() -> Self {
        Self {
            map: HashMap::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.map.is_empty()
    }
}

impl<Set, Coeffs> Mul<Set> for LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
    Set: MulAssign<Set> + Clone,
{
    type Output = Self;

    fn mul(self, rhs: Set) -> Self::Output {
        let mut new_self = Self::zero();
        #[allow(clippy::suspicious_arithmetic_impl)]
        for (mut key, value) in self.map {
            key *= rhs.clone();
            new_self
                .map
                .entry(key)
                .and_modify(|old_value| *old_value += value.clone())
                .or_insert(value);
        }
        new_self
    }
}

impl<Set, Coeffs> MulAssign<Set> for LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
    Set: MulAssign<Set> + Clone,
{
    fn mul_assign(&mut self, rhs: Set) {
        let mut dummy = Self::zero();
        core::mem::swap(self, &mut dummy);
        *self = dummy * rhs;
    }
}

#[allow(clippy::implicit_hasher)]
impl<Set, Coeffs> From<LinearCombination<Set, Coeffs>> for HashMap<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
{
    fn from(value: LinearCombination<Set, Coeffs>) -> Self {
        value.map
    }
}

impl<Set, Coeffs> TryFrom<Vec<(Set, Coeffs)>> for LinearCombination<Set, Coeffs>
where
    Set: std::hash::Hash + Eq,
    Coeffs: Add<Coeffs, Output = Coeffs> + AddAssign<Coeffs> + Clone,
{
    type Error = ();

    fn try_from(value: Vec<(Set, Coeffs)>) -> Result<Self, Self::Error> {
        let mut to_return = HashMap::new();
        for (key, v) in value {
            to_return
                .entry(key)
                .and_modify(|old_value| *old_value += v.clone())
                .or_insert(v);
        }
        Ok(Self { map: to_return })
    }
}
