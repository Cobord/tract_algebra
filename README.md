# tract_algebra

A Rust library for **Baker-Bowler tracts** — an algebraic structure that generalizes fields, hyperfields, and partial fields — together with analogs of projective varieties defined over them.

---

## Baker-Bowler Tracts

### The Concept

A **Baker-Bowler tract** (or simply a *tract*) is a set `T` equipped with:

- A **multiplicative group** structure on the non-absorbing elements `T \ {0}`, where `0` denotes the *absorbing element* (the analog of the field zero). The absorbing element satisfies `0 · t = 0` for all `t`, and it is the unique element excluded from the multiplicative group. The set `T` as a whole has the structure of a multiplicative semigroup.
- A designated element **`neg_one`** — the unique `g ∈ T \ {0}` such that the formal sum `{1, g}` belongs to the *null set*.
- A **null set** `N[T]`, a collection of formal finite multisets of non-absorbing elements, which plays the role of "sums that equal zero". In a field, a finite sum `a₁ + … + aₙ = 0` is the null set condition; in a tract the null set generalizes this without requiring a well-defined addition. The null set must be invariant under the `T*`-action: if `{a₁, …, aₙ} ∈ N[T]` then `{λa₁, …, λaₙ} ∈ N[T]` for any `λ ∈ T \ {0}`.

### The `BakerBowlerTract` Trait

```rust
pub trait BakerBowlerTract:
    MulMonoidBounds + MulGroupBounds + PartialEq + Eq + Clone + std::hash::Hash
{
    type NullSet: Zero
        + Add<Self::NullSet, Output = Self::NullSet>
        + TryFrom<Vec<(Self, usize)>>
        + Mul<Self, Output = Self::NullSet>
        + Into<HashMap<Self, usize>>
        + Clone;

    fn absorbing_element() -> Self;
    fn neg_one() -> Self;
    fn in_null_set(element: Self::NullSet) -> bool;
}
```

- **`absorbing_element()`** — the unique element that annihilates multiplication, giving the `0`.
- **`neg_one()`** — the unique `g` such that `{1, g}` is null. In a field this is `-1`; in tracts without a sign distinction it may equal `1`.
- **`in_null_set(element)`** — the membership predicate. The argument is a `Self::NullSet` which holds a formal multiset (each tract element appears with a multiplicity) representing a formal sum.

The associated `NullSet` type is concretely `LinearCombination<T, usize>` in all implementations here: a `HashMap` from tract elements to their multiplicities, supporting addition (merging multisets by summing multiplicities) and right-multiplication by a tract element (shifting all keys). However in other contexts it may be better to use a different sort of structure to store such linear combinations.

---

## Specific Tracts

### `FieldTract<F>`

The **field tract** wraps an ordinary field `F` (or any type satisfying `FieldBounds`: a commutative ring with division).

- **Elements**: elements of `F`.
- **Absorbing element**: `F::zero()`.
- **`neg_one`**: the additive inverse of `F::one()`, i.e. `-1 ∈ F`.
- **Null set**: a formal sum `Σ nᵢ aᵢ` (with multiplicities `nᵢ ∈ ℕ`) is null iff `Σ nᵢ · aᵢ = 0` in `F` (where `n · a` is repeated field addition).

`FieldTract<i128>` is used in the tests for the projective variety, giving classical algebraic geometry over a field of characteristic `0` but testing with integer points to avoid floating point errors.

### `Krasner`

The **Krasner tract** (also known as the *Krasner hyperfield*) is the smallest non-trivial tract:

```rust
pub struct Krasner(bool);
```

- **Elements**: `{false, true}`, representing `{0, 1}` — the absorbing element and a single non-absorbing element.
- **Absorbing element**: `Krasner(false)`.
- **Multiplicative group**: the trivial group `{true}`, since `1 · 1 = 1`.
- **`neg_one`**: `Krasner(true) = 1`, because in the Krasner hyperfield `-1 = 1` (there is no sign distinction).
- **Null set**: a formal multiset of elements is null iff the total multiplicity of the non-absorbing element `1` is **not equal to 1** — i.e., is zero or at least two. Geometrically: a sum is "zero" unless exactly one nonzero term survives.

`Krasner` is the **terminal object** in the category of tracts: there is a unique tract morphism from any tract `T` to `Krasner`, given by `final_morphism(t) = 0` if `t` is absorbing, and `= 1` otherwise. This "forgets" all multiplicative structure and retains only the support (zero/nonzero distinction).

### `RegularPartialField`

The **regular partial field** (also called the *initial tract*) is the partial field `{1, -1, 0}`:

```rust
pub struct RegularPartialField(Option<bool>);
```

- **Elements**: `Some(true) = 1`, `Some(false) = -1`, `None = 0` (absorbing).
- **Multiplication**: `1·x = x`, `(-1)·(-1) = 1`, `0·x = 0`.
- **`neg_one`**: `Some(false) = -1`.
- **Null set**: a formal multiset is null iff the total multiplicity of `1` equals the total multiplicity of `-1`. This is the sign-cancellation condition: there are equal numbers of positive and negative terms.

`RegularPartialField` is the **initial object** in the category of tracts: there is a unique tract morphism `initial_morphism: RegularPartialField → T` for any tract `T`, sending `1 ↦ T::one()`, `-1 ↦ T::neg_one()`, `0 ↦ T::absorbing_element()`. Every tract receives a canonical map from `RegularPartialField`.

### `TropicalTract<R, MAX_VS_MIN>`

The **tropical tract** implements the tropical semiring in log-coordinates, as a Baker-Bowler tract when `R` supports subtraction.

```rust
pub struct TropicalTract<R: TropicalBounds, const MAX_VS_MIN: bool> {
    value: Option<R>,
}
```

- **Absorbing element**: `None` (representing `±∞` depending on `MAX_VS_MIN`).
- **Multiplication** (tract multiplication = tropical-semiring multiplication): addition of exponents, `a ⊙ b = a + b`. The identity is `R::zero()`. The absorbing element `±∞` annihilates.
- **Addition** (tropical-semiring addition, used only for the `NullSet` type): `a ⊕ b = max(a, b)` (or `min`, depending on `MAX_VS_MIN`).
- **`neg_one`**: `TropicalTract::one()`There is no sign distinction.
- **Null set**: a formal multiset is null iff the extremal value (maximum for `MAX_VS_MIN = true`, minimum for `MAX_VS_MIN = false`) appears with total multiplicity **at least 2**. This is the tropical analog of zero: in a tropical polynomial, a value is a "root" when the maximum is achieved by at least two terms — the tropical cancellation condition.

### `AbsoluteValueTract<R>`

The **absolute value tract** represents the quotient `ℂ* / U(1)`, i.e. the positive reals `ℝ>0` under multiplication, equipped with the (possibly degenerate) polygon null set.

```rust
pub struct AbsoluteValueTract<R: AbsoluteValueBounds> {
    value: Option<R>,
}
```

Elements are cosets of U(1) in `ℂ*`, identified with their modulus (a positive real). The carrier `R` stores the modulus.

- **Absorbing element**: `None` (representing `0 ∈ ℂ`).
- **Multiplication**: product of moduli.
- **`neg_one`**: `AbsoluteValueTract::one()` (modulus 1), because `|-1| = |1| = 1`. The pair `{1, 1}` satisfies `2·max(1,1) = 2 = 1 + 1`, which satisfies the null set condition (degenerate polygon with two equal sides).
- **Null set** (the *possibly degenerate polygon inequality*): a multiset of moduli `{r₁, …, rₙ}` (with multiplicities) is null iff there exist directions on the unit circle such that the corresponding complex numbers sum to zero.


### `PhaseTract<S>`

The **phase tract** (or *phase hyperfield*) is the quotient `ℂ* / ℝ>0 = U(1)`, keeping only the argument (angle) and discarding the modulus.

```rust
pub struct PhaseTract<S: PhaseBounds> {
    value: Option<S>,
}
```

- **Absorbing element**: `None` (representing `0 ∈ ℂ`).
- **Multiplication**: multiplication of phases in U(1) (addition of angles).
- **`neg_one`**: the element at angle `π` (the antipodal point of `1 ∈ U(1)`).
- **Null set**: a formal multiset of phases `{e^{iθ₁}, …, e^{iθₙ}}` is null iff `0` lies in the (closed) convex hull of the corresponding points on the unit circle. This is the condition that the phases can be realized as directions of complex numbers summing to zero.

---

## Projective Varieties over Tracts

### The Problem with Ideals

Over a field `k`, an algebraic variety is classically defined as the vanishing locus of an ideal in the polynomial ring `k[x₀, …, xₙ]`. Over a tract `T`, there is no polynomial ring and no well-defined addition, so there is no notion of an ideal. The sum of two polynomials is not a polynomial — only the null set membership of a formal multiset of monomials makes sense.

### `TractProjectiveVariety<K_VARIABLES, HOMOGENEITY, T>`

```rust
pub struct TractProjectiveVariety<
    const K_VARIABLES: usize,
    const HOMOGENEITY: usize,
    T: BakerBowlerTract,
> {
    coordinates: [T; K_VARIABLES],
    relations: Vec<LinearCombination<[usize; HOMOGENEITY], i64>>,
}
```

A point in **projective `(K_VARIABLES - 1)`-space** over the tract `T` is given by:
- A coordinate tuple `[t₀, …, t_{K-1}] ∈ T^K` with at least one non-absorbing entry (since the all-absorbing point does not exist in projective space).
- A list of **relations**: homogeneous polynomial equations of degree `HOMOGENEITY`, each encoded as a `LinearCombination<[usize; HOMOGENEITY], i64>`.

### Encoding Polynomials as Relations

A monomial of degree `d` in variables `x₀, …, x_{K-1}` is encoded as an array `[i₁, i₂, …, i_d]` of variable indices (repeated indices give higher powers). An integer coefficient `c ∈ ℤ` records the signed multiplicity of this monomial. A **relation** is a formal `ℤ`-linear combination of monomials: a `LinearCombination<[usize; HOMOGENEITY], i64>`.

For example, the cuspidal cubic `x³ − y²z = 0` in `ℙ²` with variables `(x, y, z) = (0, 1, 2)` is encoded as:
```
[0,0,0] ↦ +1    (the monomial x³, coefficient +1)
[1,1,2] ↦ −1    (the monomial y²z, coefficient −1)
```

### Evaluating a Relation at a Point

Given coordinates `[t₀, …, t_{K-1}]` and a relation, each monomial `[i₁, …, i_d]` evaluates to the tract element `t_{i₁} · t_{i₂} · … · t_{i_d}` (product in the multiplicative group of `T`). The sign of the coefficient determines whether to contribute the raw value or `val · T::neg_one()`. The absolute value `|c|` is the multiplicity in the formal multiset.

The relation is **satisfied** at the point iff the resulting element of `LinearCombination<T, usize>` lies in the null set of `T`, i.e. `T::in_null_set(formal_sum)` returns `true`. A point lies on the variety iff all relations are satisfied.

### Projective Well-Definedness

Scaling all coordinates by `λ ∈ T \ {0}` multiplies every degree-`d` monomial value by `λ^d`. The null set is invariant under the `T*`-action (this is a tract axiom), so null set membership is preserved. The variety condition is therefore **projectively well-defined**: it does not depend on the choice of representative within a projective equivalence class.

### The Same Equations, Different Tracts

The same set of encoded relations defines geometrically different loci depending on the tract, because each tract has a different null set predicate. The tract variety condition is in general a *combinatorial thickening* of the field variety condition:

**Cuspidal cubic `x³ − y²z = 0` in `ℙ²`:**

| Tract | Condition on `[x:y:z]` |
|---|---|
| `FieldTract<ℤ>` | `x³ = y²z` exactly |
| `AbsoluteValueTract` | `2·max(|x|³, |y|²|z|) ≤ |x|³ + |y|²|z|`, i.e. `|x|³ = |y|²|z|` (two-term polygon collapses to equality) |
| `Krasner` | the count of nonzero monomials among `{x³, y²z}` is not equal to 1 |

**Conic `z² − x² − y² = 0` in `ℙ²`:**

| Tract | Condition on `[x:y:z]` |
|---|---|
| `FieldTract<ℤ>` | `z² = x² + y²` (Pythagorean triples) |
| `AbsoluteValueTract` | `2·max(z², x², y²) ≤ z² + x² + y²` (polygon inequality — includes `[1:1:1]`, which is not a field point) |
| `Krasner` | count of nonzero monomials among `{z², x², y²}` is not equal to 1 |

---

## Summary of Tract Morphisms

```
RegularPartialField  (initial object)
         │
         │  initial_morphism  (canonical map to every tract)
         ↓
      FieldTract<F>
      TropicalTract<R>
      AbsoluteValueTract<R>
      PhaseTract<S>
         │
         │  final_morphism  (canonical map from every tract)
         ↓
       Krasner  (terminal object)
```

Every tract receives a canonical map from `RegularPartialField` (sending `1 ↦ 1`, `-1 ↦ neg_one`, `0 ↦ 0`) and maps canonically to `Krasner` (sending non-absorbing elements to `1` and the absorbing element to `0`).
