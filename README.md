# dissociation

Acid–base equilibrium without solving a quadratic by hand.

Asking a solution for its pH is asking an equilibrium where it settled.
A strong acid answers straight away, because it gives up every proton it carries and `[H+]` is the concentration itself.
A weak acid keeps most of them, so only a fraction dissociates,
and that fraction moves with the concentration it was dissolved at.
Reading the pH off the bottle is impossible; the equilibrium has to be solved.

Its origin is a shortcut found in a university chemistry course around 2007,
where the assigned method wanted the quadratic worked out for every exercise,
usually after dropping a term to keep the arithmetic short.
The shortcut turned out to be published — in 1993, and again in 2023 — which
[the shortcut was published first](#the-shortcut-was-published-first) cites in full.
This module solves that quadratic once, symbolically, and ships the root as a function.
Every exercise after it is a call: exact, with nothing approximated away,
and the same treatment carries on into polyprotic acids, buffers and titration curves.

The public API lives in [`dissociation/`](dissociation/), split by equilibrium,
buffers, titration, species and examples. PyYAML loads the packaged constants catalog.

## The doors

```
./build.sh    ./test.sh     # unix
build.cmd     test.cmd      # windows
```

`build.sh` fetches its own tooling into `.venv`, then runs lint, types, tests and the wheel.
It refuses to build what does not pass, and a refusal names the line that works.
The lint catches bugs and never taste: the aligned tables are deliberate, so no formatter reflows them.
`test.sh` runs that complete build and then the worked demo.

The `.cmd` files are shims. They name, call and propagate the exit code, and nothing else;
`run-bash.cmd` finds the interpreter once, for both of them.
They are read from here and proven only over there, so until someone runs them on Windows,
call them unverified.

## The coefficients never change

A weak monoprotic acid always yields the same quadratic, so it is worth solving once and for all.
Dissolving `HA` at analytical concentration `C` leaves three unknowns — `[HA]`, `[A−]` and `[H+]` —
and three statements that pin them down.

The [acid dissociation constant](https://en.wikipedia.org/wiki/Acid_dissociation_constant) defines the first:

```
Ka = [H+]·[A−] / [HA]
```

Mass conserves the acid, whichever of the two forms it ended up in:

```
C = [HA] + [A−]
```

Charge balances the solution, and while water stays negligible each proton answers to one anion:

```
[H+] = [A−]
```

Substituting the last two into the first leaves one unknown, `x = [H+]`,
in `Ka = x² / (C − x)`, which rearranges to

```
x² + Ka·x − Ka·C = 0
```

and the coefficients are `a = 1`, `b = Ka`, `c = −Ka·C` for every acid and every concentration.
Nothing there depends on which exercise you were given.
The quadratic formula answers it in general, and the negative root goes,
because a concentration below zero means nothing:

```
x = (−Ka + √(Ka² + 4·Ka·C)) / 2
  = (Ka/2) · (√(1 + 4C/Ka) − 1)
```

That leaves a direct operation, exact, with no `C − x ≈ C` approximation and no crank to turn.

Floating point disagrees with algebra on how to evaluate it, so `monoprotic_h` ships the rationalized form.
Multiplying numerator and denominator by `Ka + √(Ka² + 4·Ka·C)` turns the difference into a sum:

```
[H+] = 2·Ka·C / (Ka + √(Ka² + 4·Ka·C))
```

Both expressions are the same number. The first one subtracts two nearly equal quantities when `C ≪ Ka`,
so floating-point arithmetic can lose significant digits in that limit
— see [loss of significance](https://en.wikipedia.org/wiki/Loss_of_significance).

### Precision and simplicity benchmark

Let `q = C/Ka`. The exact root, in units of `Ka`, is
`x/Ka = 2q / (1 + √(1 + 4q))`. It gives a stable reference for comparing the direct formula,
the rationalized formula, and the classroom approximation `x ≈ √(Ka·C)`.

The approximation assumes `x` is small beside `C`, which holds as `q` grows.
Its relative overestimate is measured against the exact root; the floating-point columns compare
the two exact forms with an 80-digit `Decimal` reference using binary64 values for `Ka = 1` and `C = q`.

| `C/Ka` | Approximation overestimate | Direct-form relative error | Rationalized-form relative error |
| ---: | ---: | ---: | ---: |
| `1e-16` | Outside its valid regime | `+11.0223%` | `+1.00e-14%` |
| `1e-12` | Outside its valid regime | `−0.002212%` | `−4.84e-15%` |
| `1e2` | `+5.1249%` | `+9.25e-15%` | `−9.42e-15%` |
| `1e4` | `+0.50125%` | `−1.97e-15%` | `−1.97e-15%` |
| `1e6` | `+0.05001%` | `−6.43e-16%` | `−6.43e-16%` |

For a simple calculation benchmark, count arithmetic operators and square roots in the normalized formulas,
excluding constants: the direct exact form uses four arithmetic operations and one square root;
the rationalized exact form uses five and one; `√q` uses one square root and no arithmetic operations.
That last count buys simplicity by assuming the approximation. The direct factored form is the shortest exact expression;
the rationalized form adds one operation and avoids cancellation when `C ≪ Ka`.
The approximation's error depends on `C/Ka`; it is within 5% here once `C/Ka` is about `105` or larger.
For the common weak-acid regime `Ka ≪ C`, the direct exact formula is already well-conditioned,
and the approximation improves as `C/Ka` increases.

## Many constants become a sum

A polyprotic acid generalizes the closed form into a sum, because `n` constants raise the polynomial to degree `n+1`.
Solving that symbolically buys nothing past the quadratic, so the species carry the structure instead.

Cumulative constants collapse the stepwise ones into a single product per species,
`βⱼ = Ka₁·Ka₂···Kaⱼ` with `β₀ = 1`.
Each `βⱼ` weighs the form that has lost `j` protons, so at a given `h = [H+]` the fractions normalize to one:

```
αⱼ = βⱼ·h^(n−j) / Σₖ βₖ·h^(n−k)
```

A species that lost `j` protons handed `j` of them to the solution, so the acid contributes `Σⱼ j·αⱼ·C`,
water contributes `Kw/h`, and the proton balance closes the system:

```
h = Σⱼ j·αⱼ·C + Kw/h
```

One constant recovers the quadratic, plus the `Kw/h` term the closed form leaves out,
which is why `exact=True` changes the answer only where water matters.

`polyprotic_h` solves that balance by [bisection](https://en.wikipedia.org/wiki/Bisection_method)
on a logarithmic scale, which is monotone and always converges.
The [self-ionization of water](https://en.wikipedia.org/wiki/Self-ionization_of_water) enters the balance here,
so extreme dilutions approach pH 7 instead of drifting past it.
`fractions` exposes the αⱼ on their own, for a distribution diagram.

## A base closes against pKw

A weak base runs the same algebra on `[OH−]` and finishes with `pH = pKw − pOH`.
`base_ph` keeps `kw` as a parameter rather than hardcoding the 14,
so a temperature other than 25 °C costs one argument.

## The table answers by name

Twenty-three common species carry their constants, so a calculation needs no lookup elsewhere.
`ACIDS` holds fifteen acids and `BASES` holds eight bases, all at 25 °C,
rounded to two significant figures from the standard tables
(see the [CRC Handbook of Chemistry and Physics](https://en.wikipedia.org/wiki/CRC_Handbook_of_Chemistry_and_Physics)).
The values, formulas and aliases live in
[`common-species.yaml`](dissociation/constants/common-species.yaml); Python only validates them
and turns them into immutable `Species` values. Aliases are grouped under `en`, `es` and `pt` in
the YAML, then flattened into `Species.alias` for lookup compatibility.

Every species is named in English, and Spanish and Portuguese ride along as aliases.
`find_species` matches on formula, name or alias, and folds accents and case before comparing,
so `"acetic acid"`, `"CH3COOH"`, `"ácido acético"`, `"ACIDO ACETICO"`, `"vinegar"` and `"vinagre"`
all reach the same entry. The folding is why most Spanish and Portuguese names need no second alias:
`"ácido carbónico"` and `"ácido carbônico"` both normalize to `acido carbonico`.
Only the forms that stay distinct after folding are listed separately — the Portuguese
`fluorídrico`, `cianídrico` and `sulfídrico` against the Spanish `fluorhídrico`, `cianhídrico`
and `sulfhídrico`, and `amônia` against `amoníaco`.

| English | Spanish | Portuguese | Formula |
| --- | --- | --- | --- |
| hydrofluoric acid | ácido fluorhídrico | ácido fluorídrico | `HF` |
| hydrocyanic acid | ácido cianhídrico | ácido cianídrico | `HCN` |
| hydrosulfuric acid | ácido sulfhídrico | ácido sulfídrico | `H2S` |
| carbonic acid | ácido carbónico | ácido carbônico | `H2CO3` |
| ammonia | amoníaco | amônia | `NH3` |
| acetic acid | ácido acético, vinagre | ácido acético, vinagre | `CH3COOH` |
| ascorbic acid | ácido ascórbico, vitamina c | ácido ascórbico, vitamina c | `C6H8O6` |

## A buffer is the same balance

A conjugate pair adds one spectator cation, and nothing else about the balance changes.
`buffer_ph` feeds `base_c` in as that cation and solves exactly,
so the answer stays right where
[Henderson–Hasselbalch](https://en.wikipedia.org/wiki/Henderson%E2%80%93Hasselbalch_equation) drifts.
Both ship, because the approximation is the one people quote:
a 1:1 acetate pair at 0.1 M gives 4.745 either way, and the same pair at 1×10⁻⁵ M gives 5.268 exactly
against 4.745 approximated.

## One balance draws the whole curve

A titration needs no piecewise treatment, because dilution and the strong titrant fit inside the same equation.
`titration_ph` dilutes analyte and titrant by the combined volume and enters the strong base as its cation,
so the buffer region, the equivalence point and the excess-titrant tail all come from one solve.
`titration_curve` samples it, `equivalence_volume` locates the jumps,
and volumes carry whatever unit you pass in.

Acetic acid, 0.1 M in 50 mL, against 0.1 M NaOH:
pH 2.875 at the start, 4.745 at half-equivalence, 8.722 at equivalence, 11.959 at 60 mL.
Phosphoric acid shows its first two `pKa` plateaus at 25 and 75 mL and flattens on the third,
which is what a `Ka₃` of 4.2×10⁻¹³ does in water.

## The shortcut was published first

The idea here was derived around 2007, in a university course, and de Levie had published it in 1993.
Nobody asked at the time recognized it — one professor urged publishing it, which is the honest
reaction of someone who had not seen it either.
Fourteen years of priority and a full literature sat behind a vocabulary nobody in that room had:
*closed-form*, *without segmentation*, *proton condition*.
Finding it took until 2026, and finding it settles the question: this is a rediscovery, not a discovery.

The work stands where it was published, and the citations go here so nobody repeats the search:

- One equation for the whole titration curve is de Levie,
  *Explicit expressions of the general form of the titration curve in terms of concentration:
  Writing a single closed-form expression for the titration curve for a variety of titrations
  without using approximations or segmentation*,
  [J. Chem. Educ. 70, 209 (1993)](https://doi.org/10.1021/ed070p209).
  That subtitle names what `titration_ph` does, thirty-three years earlier.
- The monoprotic closed form, the buffer and the titration all come out of one expression in
  Caicedo, Morales and Arango, *Closed-form expressions for monoprotic weak acid aqueous solutions*,
  [ChemTexts 9, 3 (2023)](https://doi.org/10.1007/s40828-023-00181-w).
  They derive it with a computer algebra system and then measure what the classroom approximations cost,
  which is the comparison `henderson_hasselbalch` keeps alive here.
- The `αⱼ` and their cumulative `βⱼ` are the standard
  [distribution diagram](https://en.wikipedia.org/wiki/Acid_dissociation_constant#Polyprotic_acids) formalism,
  borrowed from the [stability constants](https://en.wikipedia.org/wiki/Stability_constants_of_complexes) of complexes,
  where a cumulative constant is always the product of the stepwise ones.
- Closing the system with mass balance, charge balance and `Kw` is the
  [systematic treatment](https://en.wikipedia.org/wiki/PH#pH_calculations) that analytical chemistry
  has taught for decades — Butler's *Ionic Equilibrium* and Harris's *Quantitative Chemical Analysis*
  both carry it. The approximation is the shortcut; solving is the default it departed from.
- Rationalizing the root is the textbook fix for the unstable quadratic formula,
  not a chemistry result at all — see [loss of significance](https://en.wikipedia.org/wiki/Loss_of_significance).

What is left unclaimed by those papers is the software, and that is all this module claims.
de Levie published an equation and Caicedo published a derivation; neither published an implementation.
A small package, one declarative constants catalog, and the exact path taken by default where a
course would have reached for the approximation.

That a correct method reached print in 1993 and still failed to reach a classroom in 2007
is a fact about how chemistry is taught, never a fact about the method.

## Use

The table answers by name, and `ph` reads the entry to decide whether it is an acid or a base:

```python
from dissociation import ph, acid_ph, base_ph, find_species

ph(
    species="acetic acid",
    c=0.1,
)  # 2.875 — 0.1 M
ph(
    species="CH3COOH",
    c=0.1,
)  # 2.875 — the formula reaches the same entry
ph(
    species="ammonia",
    c=0.1,
)  # 11.125 — a base, closed against pKw
ph(
    species="phosphoric acid",
    c=0.1,
)  # 1.622 — three constants, solved as one balance

find_species("acetic acid")
# <acetic acid (CH3COOH) Ka=[1.8e-05]>
```

Spanish and Portuguese resolve to the same entries, accented or not:

```python
ph(
    species="ácido acético",
    c=0.1,
)  # 2.875 — es and pt spell this one alike
ph(
    species="ACIDO ACETICO",
    c=0.1,
)  # 2.875 — accents and case are folded
ph(
    species="vinagre",
    c=0.1,
)  # 2.875 — es and pt nickname; "vinegar" in English

ph(
    species="ácido fluorhídrico",
    c=0.1,
)  # 2.102 — es
ph(
    species="ácido fluorídrico",
    c=0.1,
)  # 2.102 — pt, one letter apart
ph(
    species="hydrofluoric acid",
    c=0.1,
)  # 2.102 — en

ph(
    species="amoníaco",
    c=0.1,
)  # 11.125 — es
ph(
    species="amônia",
    c=0.1,
)  # 11.125 — pt
ph(
    species="ammonia",
    c=0.1,
)  # 11.125 — en

find_species("ácido cianídrico")
# <hydrocyanic acid (HCN) Ka=[6.2e-10]>
```

A constant you bring yourself skips the table entirely:

```python
from dissociation import pk, k_from_pk

acid_ph(
    ka=1.8e-5,
    c=0.1,
)  # 2.875 — the closed form
base_ph(
    kb=1.8e-5,
    c=0.1,
)  # 11.125 — the same number read as a Kb
base_ph(
    kb=1.8e-5,
    c=0.1,
    kw=k_from_pk(13.6),
)  # 10.725 — water at another temperature

pk(1.8e-5)                  # 4.745
k_from_pk(4.76)             # 1.7378e-05 — the round trip
```

`exact=True` trades the closed form for the full balance, which matters once water competes:

```python
acid_ph(
    ka=4.9e-10,
    c=1e-6,
)  # 7.660 — the closed form, past neutral
acid_ph(
    ka=4.9e-10,
    c=1e-6,
    exact=True,
)  # 6.990 — water included, where it belongs
```

`fractions` answers what the solution is made of at a given `[H+]`, which draws a distribution diagram:

```python
from dissociation import fractions

fractions(find_species("H3PO4").k, 1e-7)
# [0.0, 0.6173, 0.3827, 0.0] — speciation at pH 7

fractions(find_species("H2CO3").k, 10 ** -7.4)
# [0.0846, 0.9143, 0.0011] — carbonate at blood pH, almost all bicarbonate
```

A buffer solves the same balance with the spectator cation added, so it stays right where the approximation drifts:

```python
from dissociation import buffer_ph, henderson_hasselbalch

buffer_ph(
    ka=1.8e-5,
    acid_c=0.1,
    base_c=0.1,
)  # 4.745 — exact
henderson_hasselbalch(
    ka=1.8e-5,
    acid_c=0.1,
    base_c=0.1,
)  # 4.745 — approximated
buffer_ph(
    ka=1.8e-5,
    acid_c=0.1,
    base_c=0.05,
)  # 4.444, against 4.444 approximated
buffer_ph(
    ka=1.8e-5,
    acid_c=1e-5,
    base_c=1e-5,
)  # 5.268, approximation says 4.745
```

A titration comes out of one solve per point, whichever region the point falls in:

```python
from dissociation import titration_ph, titration_curve, equivalence_volume

equivalence_volume(
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
)  # 50.0 mL
equivalence_volume(
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    protons=3,
)  # 150.0 mL — all three of H3PO4

titration_ph(
    kas=[1.8e-5],
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    titrant_volume=0.0,
)  # 2.875 — the free acid
titration_ph(
    kas=[1.8e-5],
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    titrant_volume=25.0,
)  # 4.745 — half-equivalence, pH = pKa
titration_ph(
    kas=[1.8e-5],
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    titrant_volume=50.0,
)  # 8.722 — equivalence, basic as it should be
titration_ph(
    kas=[1.8e-5],
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    titrant_volume=60.0,
)  # 11.959 — excess titrant

titration_curve(
    kas=[1.8e-5],
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    points=5,
)
# [(0.0, 2.875), (25.0, 4.745), (50.0, 8.722), (75.0, 12.301), (100.0, 12.523)]
```

A polyprotic curve needs nothing new, only more constants:

```python
kas = find_species("H3PO4").k
titration_ph(
    kas=kas,
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    titrant_volume=25.0,
)  # 2.267 — the first plateau
titration_curve(
    kas=kas,
    analyte_c=0.1,
    analyte_volume=50.0,
    titrant_c=0.1,
    up_to=200.0,
    points=201,
)
```

The module runs its own demo, through the door, the installed command or the file:

```
./test.sh
dissociation
python dissociation.py
python -m dissociation
```

## Scope

Concentrations stand in for activities, so
[activity coefficients](https://en.wikipedia.org/wiki/Activity_coefficient)
and [ionic strength](https://en.wikipedia.org/wiki/Ionic_strength) are absent,
which is the usual textbook footing and holds well enough below roughly 0.1 M.
A titration assumes a strong titrant; a weak-against-weak one is not covered.
What is here reaches one acid or base in water, its conjugate pair and its titration curve,
and within that boundary it solves rather than approximates.

## Install

The package ships, imports and runs under one name:

```
pip install dissociation
```

`pyproject.toml` packages the modules and YAML catalog together and declares PyYAML as the loader.
`import dissociation` exposes the same public API from a checkout or an installation.
Installing also puts a `dissociation` command on the PATH, which runs `demo()` — the same worked
example, reachable once the file is no longer in front of you.

## Tests

`test_dissociation.py` holds 66 tests on `unittest`; the suite needs no extra framework to be
believed. It enters through the public functions only
and spells every expectation by hand, so it never agrees with the code it watches.
Eight mutations were applied and reverted; the one that survived the first pass
(`14.0` hardcoded in place of `pKw`) bought the test that now kills it.

```
python -m unittest discover -p 'test_*.py'
```

The floor in `requires-python` is measured, not guessed.
The suite passes on 3.9 through 3.14, and on 3.8 as well, which is older than the floor admits:

```
for v in 3.9 3.10 3.11 3.12 3.13 3.14; do
	uv run --isolated --no-project --python $v python -m unittest discover -p 'test_*.py'
done
```

## State

Identifiers, prose and tabulated names are English throughout, and Spanish and Portuguese
survive only in the aliases.
The file was renamed from `disosiacion.py`, whose spelling missed a `c` and disagreed with its own docstring;
that file is now gone.

## License

[MIT](LICENSE). Permissive on purpose: research code earns its keep by being reused,
and a copyleft term would bar it from the BSD-licensed scientific stack it belongs beside.
