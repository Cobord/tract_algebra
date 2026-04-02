use crate::{
    baker_tract::BakerBowlerTract,
    linear_combination::LinearCombination,
    tract_projective_variety::{TractProjectiveVariety, TractProjectiveVarietyError},
};

pub struct BinomialEdgeIdeal<T, const N: usize, const TWO_N: usize>
where
    T: BakerBowlerTract,
{
    as_variety: TractProjectiveVariety<TWO_N, 2, T>,
    #[allow(dead_code)]
    edges: Vec<(usize, usize)>,
}

impl<T, const N: usize, const TWO_N: usize> BinomialEdgeIdeal<T, N, TWO_N>
where
    T: BakerBowlerTract,
{
    pub fn new(
        coords: [T; TWO_N],
        mut edges: Vec<(usize, usize)>,
        skip_valid: bool,
    ) -> Result<Self, TractProjectiveVarietyError> {
        let mut relations: Vec<LinearCombination<_, i64>> = Vec::with_capacity(edges.len());
        edges.retain_mut(|item| {
            if item.0 == item.1 {
                false
            } else {
                if item.0 > item.1 {
                    std::mem::swap(&mut item.0, &mut item.1);
                }
                true
            }
        });
        edges.sort_unstable();
        edges.dedup();
        for (i, j) in &edges {
            let relation_ij: Vec<([usize; 2], i64)> = vec![([*i, j + N], 1), ([*j, *i + N], -1)]
                .into_iter()
                .collect();
            relations.push(relation_ij.try_into().expect("This should never fail since the relation is always a linear combination of two monomials"));
        }
        let as_variety = TractProjectiveVariety::new(coords, relations, skip_valid)?;
        Ok(Self { as_variety, edges })
    }

    pub fn as_variety(&self) -> &TractProjectiveVariety<TWO_N, 2, T> {
        &self.as_variety
    }
}

#[cfg(test)]
mod tests {

    use num::Zero;
    use proptest::{prop_assert, proptest};

    use crate::linear_combination::LinearCombination;

    #[test]
    fn small_binomomial_edge_ideal() {
        use super::BinomialEdgeIdeal;
        use crate::field_tract::FieldTract;
        use crate::tract_projective_variety::TractProjectiveVarietyError;
        use num::One;
        let coords = core::array::from_fn(|_| FieldTract::<i64>::one());
        let edges = vec![(0, 1), (1, 2), (2, 3)];
        let ideal = BinomialEdgeIdeal::<_, 4, 8>::new(coords, edges, false)
            .expect("Failed to create binomial edge ideal");
        assert_eq!(ideal.edges.len(), 3);

        let alt_coordinates = core::array::from_fn(|idx| FieldTract::<i64>::from(idx as i64));
        let bad_coord =
            BinomialEdgeIdeal::<_, 4, 8>::new(alt_coordinates, vec![(0, 1), (1, 2), (2, 3)], false)
                .map(|_| ())
                .expect_err("This is not a point on the variety");
        assert_eq!(
            matches!(
                bad_coord,
                TractProjectiveVarietyError::RelationNotSatisfied { relation_index: _ }
            ),
            true
        );

        let alt_coordinates = core::array::from_fn(|idx| FieldTract::<i64>::from(idx as i64));
        let bad_coord =
            BinomialEdgeIdeal::<_, 4, 8>::new(alt_coordinates, vec![(0, 1), (1, 2), (2, 3)], true)
                .expect(
                    "This is not a point on the variety but we skipped the check for relations",
                );
        let invalid_relation = bad_coord.as_variety().invalid_relation().expect_err(
            "This is not a point on the variety, so there should be an invalid relation",
        );
        let mut expected_failed_relation: LinearCombination<[usize; 2], i64> =
            vec![([1, 4], 1), ([0, 5], -1)]
                .try_into()
                .expect("valid linear combination");
        expected_failed_relation.rescale(&-1);
        let mut diff = invalid_relation.clone() - expected_failed_relation;
        diff.simplify();
        assert!(diff.is_zero(), "{:?}", diff);
    }

    proptest! {

        #[test]
        fn small_binomomial_edge_ideal_prop(
            a in -1000i64..1000,
            b in -1000i64..1000,
            c in -1000i64..1000,
            d in -1000i64..1000,
            edges in proptest::collection::vec((0usize..4, 0usize..4), 1..10)
        ) {
            use super::BinomialEdgeIdeal;
            use crate::field_tract::FieldTract;

            let foo = |x : usize| -> i64 {
                match x {
                    0 => a,
                    1 => b,
                    2 => c,
                    3 => d,
                    _ => unreachable!(),
                }
            };

            let alt_coordinates = core::array::from_fn(|idx| FieldTract::<i64>::from(foo(idx % 4)));
            let num_edges = edges.len();
            let good_coord = BinomialEdgeIdeal::<_, 4, 8>::new(alt_coordinates, edges, false)
                .expect("This is a point on the variety");
            prop_assert!(good_coord.edges.len() <= num_edges);
        }
    }
}
