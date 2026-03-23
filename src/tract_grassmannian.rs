/// Consider any set of homogeneous degree d polynomials
/// with integer coefficients in some number K of variables.
/// In honest fields, we would construct the ideal they define
/// and vanishing locus of this set would be that everything in the
/// ideal vanishes on this locus of F^K and that locus would be invariant under
/// scaling and so it would define some variety in FP^{K-1}
/// However with tracts we are looking at each polynomial instead of just the ideal.
/// `f = \sum_t a_t x_{i_1,t} ... x_{i_degree,t}`
/// The tract does not really have enough information to say yes a sum is zero
/// but it does have enough information to say if something in N[F] is in a
/// null set which is to say if we could sum properly, could this sum be 0.
/// Asking for this condition instead gives a subset of F^K again.
/// Scaling by a common thing z in F^* picks up z^d for the value and
/// because the null set \subset N[F] is invariant under F^* multiplication action
/// we can again consider the locus in `(F^K - 0) / F^*` which is
/// the analog of `FP^{K-1}`. Here 0 means the all absorbing element of F because you don't
/// have additive identity without the addition anymore.
use std::collections::HashMap;

use itertools::Itertools;

use crate::{baker_tract::BakerBowlerTract, krasner::Krasner};

pub struct TractGrassmannian<
    const SUBSPACE_ISH: usize,
    const AMBIENT_ISH: usize,
    T: BakerBowlerTract,
> {
    plucker_coordinates: HashMap<[usize; SUBSPACE_ISH], T>,
}

impl<const SUBSPACE_ISH: usize, const AMBIENT_ISH: usize, T: BakerBowlerTract>
    TractGrassmannian<SUBSPACE_ISH, AMBIENT_ISH, T>
{
    #[must_use = "This is the number of Plucker coordinates, though they might not all be stored if some were the absorbing element"]
    #[allow(clippy::many_single_char_names)]
    pub const fn how_many_pluckers() -> usize {
        assert!(SUBSPACE_ISH <= AMBIENT_ISH);
        let mut n = AMBIENT_ISH;
        let mut k = SUBSPACE_ISH;
        if k > n - k {
            k = n - k;
        }
        let mut r = 1usize;
        let mut d = 1usize;
        while d <= k {
            r = match r.checked_mul(n) {
                Some(v) => v,
                None => panic!("Dimensions are too large"),
            };
            r /= d;
            n -= 1;
            d += 1;
        }
        r
    }

    #[allow(clippy::result_unit_err)]
    pub fn new(
        plucker_coordinates: HashMap<[usize; SUBSPACE_ISH], T>,
        skip_valid: bool,
    ) -> Result<Self, ()> {
        let plucker_coordinates: HashMap<_, _> = plucker_coordinates
            .into_iter()
            .filter(|(_, v)| !v.is_absorbing_elt())
            .collect();
        if plucker_coordinates.is_empty() {
            return Err(());
        }
        let to_return = Self {
            plucker_coordinates,
        };
        if skip_valid {
            return Ok(to_return);
        }
        if to_return.plucker_relations_satisfied() {
            Ok(to_return)
        } else {
            Err(())
        }
    }

    #[must_use = "What are you doing with this Plucker coordinate as an element of the tract"]
    pub fn get_plucker(&self, mut which_plucker: [usize; SUBSPACE_ISH]) -> T {
        if which_plucker.is_sorted() {
            self.plucker_coordinates
                .get(&which_plucker)
                .cloned()
                .unwrap_or(T::absorbing_element())
        } else {
            let parity: bool = {
                let mut even = true;
                for idx in 0..SUBSPACE_ISH {
                    for jdx in idx + 1..SUBSPACE_ISH {
                        if which_plucker[idx] > which_plucker[jdx] {
                            even = !even;
                        }
                    }
                }
                even
            };
            which_plucker.sort_unstable();
            let without_parity = self
                .plucker_coordinates
                .get(&which_plucker)
                .cloned()
                .unwrap_or(T::absorbing_element());
            if parity {
                without_parity
            } else {
                without_parity * T::neg_one()
            }
        }
    }

    fn plucker_relations_satisfied(&self) -> bool {
        for mut i1_through_ikm1 in (0..AMBIENT_ISH).combinations(SUBSPACE_ISH - 1) {
            i1_through_ikm1.sort_unstable();
            for mut j1_through_jkp1 in (0..AMBIENT_ISH).combinations(SUBSPACE_ISH + 1) {
                j1_through_jkp1.sort_unstable();
                if !self.plucker_relation_satisfied(&i1_through_ikm1, &j1_through_jkp1) {
                    return false;
                }
            }
        }
        true
    }

    fn plucker_relation_satisfied(
        &self,
        i1_through_ik: &[usize],
        j1_through_jkp1: &[usize],
    ) -> bool {
        let mut all_terms = Vec::new();
        let mut terms_chosen = Vec::new();
        for l in 0..=SUBSPACE_ISH {
            let i_set: [usize; SUBSPACE_ISH] = core::array::from_fn(|idx| {
                if idx < SUBSPACE_ISH - 1 {
                    i1_through_ik[idx]
                } else {
                    j1_through_jkp1[l]
                }
            });
            let j_set: [usize; SUBSPACE_ISH] = core::array::from_fn(|idx| {
                if idx < l {
                    j1_through_jkp1[idx]
                } else {
                    j1_through_jkp1[idx + 1]
                }
            });
            terms_chosen.push((l, i_set, j_set));
            let i_set_plucker = self.get_plucker(i_set);
            let j_set_plucker = self.get_plucker(j_set);
            let mut contribution_now = if l % 2 == 0 { T::neg_one() } else { T::one() };
            contribution_now *= i_set_plucker;
            contribution_now *= j_set_plucker;
            all_terms.push((contribution_now, 1usize));
        }
        let as_null_set: Result<T::NullSet, _> = all_terms.try_into();
        if let Ok(as_null_set) = as_null_set {
            T::in_null_set(as_null_set)
        } else {
            true
        }
    }

    /// Returns the underlying matroid as a `TractGrassmannian<Krasner>`.
    ///
    /// Applying `Krasner::final_morphism` to every Plücker coordinate forgets all
    /// tract-specific data and retains only the support: which index sets have a
    /// non-absorbing value.  The resulting Krasner Grassmannian encodes exactly the
    /// bases of the matroid — the collection of `SUBSPACE_ISH`-element subsets of
    /// `{0, …, AMBIENT_ISH − 1}` whose Plücker coordinate was nonzero.
    /// No separate matroid type is needed: `TractGrassmannian<Krasner>` is the matroid.
    #[must_use]
    pub fn underlying_matroid(&self) -> TractGrassmannian<SUBSPACE_ISH, AMBIENT_ISH, Krasner> {
        let plucker_coordinates = self
            .plucker_coordinates
            .iter()
            .filter_map(|(&index_set, value)| {
                let k = Krasner::final_morphism_ref(value);
                (!k.is_absorbing_elt()).then_some((index_set, k))
            })
            .collect();
        TractGrassmannian {
            plucker_coordinates,
        }
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashMap;

    use crate::krasner::Krasner;

    #[test]
    fn gr24() {
        use super::TractGrassmannian;
        use crate::field_tract::FieldTract;
        /*
        [1 0 4 5]
        [3 8 7 2]
         */
        let mut plucker_coordinates = HashMap::<[usize; 2], FieldTract<i128>>::new();
        plucker_coordinates.insert([0, 1], (8).into());
        plucker_coordinates.insert([0, 2], (-5).into());
        plucker_coordinates.insert([0, 3], (-13).into());
        plucker_coordinates.insert([1, 2], (-32).into());
        plucker_coordinates.insert([1, 3], (-40).into());
        plucker_coordinates.insert([2, 3], (8 - 35).into());

        let t = TractGrassmannian::<2, 4, FieldTract<i128>>::new(plucker_coordinates, false)
            .expect("We chose good coordinates");

        assert_eq!(t.get_plucker([0, 1]), (8).into());
        assert_eq!(t.get_plucker([1, 0]), (-8).into());
    }

    #[test]
    fn gr24_krasner() {
        use super::TractGrassmannian;
        use crate::field_tract::FieldTract;
        /*
        [1 0 4 5]
        [3 8 7 2]
         */
        let mut plucker_coordinates = HashMap::<[usize; 2], FieldTract<i128>>::new();
        plucker_coordinates.insert([0, 1], (8).into());
        plucker_coordinates.insert([0, 2], (-5).into());
        plucker_coordinates.insert([0, 3], (-13).into());
        plucker_coordinates.insert([1, 2], (-32).into());
        plucker_coordinates.insert([1, 3], (-40).into());
        plucker_coordinates.insert([2, 3], (8 - 35).into());
        let plucker_coordinates: HashMap<[usize; 2], Krasner> = plucker_coordinates
            .into_iter()
            .map(|(key, value)| (key, Krasner::final_morphism(value)))
            .collect();

        let t = TractGrassmannian::<2, 4, Krasner>::new(plucker_coordinates, false)
            .expect("We chose good coordinates");

        assert_eq!(
            t.get_plucker([0, 1]),
            Krasner::final_morphism::<FieldTract<i128>>((8).into())
        );
        assert_eq!(
            t.get_plucker([1, 0]),
            Krasner::final_morphism::<FieldTract<i128>>((-8).into())
        );
    }
}
