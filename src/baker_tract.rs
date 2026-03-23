use std::{
    collections::HashMap,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign},
};

use num::{One, Zero};

pub trait MulMonoidBounds: Mul<Self, Output = Self> + MulAssign<Self> + Sized {}

impl<T> MulMonoidBounds for T where T: Mul<T, Output = T> + MulAssign<T> {}

pub trait MulGroupBounds:
    MulMonoidBounds + Div<Self, Output = Self> + DivAssign<Self> + One
{
}

impl<T> MulGroupBounds for T where T: MulMonoidBounds + Div<T, Output = T> + DivAssign<T> + One {}

#[derive(PartialEq, Eq, Clone, Copy, Debug, thiserror::Error)]
pub enum BakerBowlerTractError {
    #[error("absorbing element is not part of the multiplicative group")]
    AbsorbingElementNotPartOfGroup,
    #[error("expected b == -a to form a null set element")]
    AssumingANegB,
}

pub trait BakerBowlerTract:
    MulMonoidBounds + MulGroupBounds + PartialEq + Eq + Clone + std::hash::Hash
{
    type NullSet: Zero
        + Add<Self::NullSet, Output = Self::NullSet>
        + AddAssign<Self::NullSet>
        + TryFrom<Vec<(Self, usize)>>
        + Mul<Self, Output = Self::NullSet>
        + MulAssign<Self>
        + Into<HashMap<Self, usize>>
        + Clone;

    fn absorbing_element() -> Self;

    fn is_absorbing_elt(&self) -> bool {
        self == &Self::absorbing_element()
    }

    fn safe_div(self, other: Self) -> Result<Self, BakerBowlerTractError> {
        if self.is_absorbing_elt() || other.is_absorbing_elt() {
            Err(BakerBowlerTractError::AbsorbingElementNotPartOfGroup)
        } else {
            Ok(self / other)
        }
    }

    fn produce_null_set_element(self, b: Self) -> Result<Self::NullSet, BakerBowlerTractError> {
        if self.is_absorbing_elt() || b.is_absorbing_elt() {
            return Err(BakerBowlerTractError::AbsorbingElementNotPartOfGroup);
        }
        let a = self;
        let neg_a = Self::neg_one() * a.clone();
        if b == neg_a {
            vec![(a, 1), (b, 1)]
                .try_into()
                .map_err(|_| BakerBowlerTractError::AbsorbingElementNotPartOfGroup)
        } else {
            Err(BakerBowlerTractError::AssumingANegB)
        }
    }

    fn neg_one() -> Self;

    fn in_null_set(element: Self::NullSet) -> bool;
}
