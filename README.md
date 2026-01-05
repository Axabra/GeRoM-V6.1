# GeRoM - Geometric Reasoning Machine

**V6.1 released – January 2026**

GeRoM is a deterministic symbolic audit engine for Navier-Stokes fields.

**Not a "thinking AI"**. No hallucinations. No magic discoveries.  
Just rigorous, traceable verification of mathematical coherence under explicit assumptions.

Open-source code: Python + SymPy + NetworkX. Runs locally (even on phone). Zero cloud.

## What GeRoM V6.1 does (and does very well)

- Checks incompressibility (∇·u = 0 – hard gate)
- Infers pressure in 2D if missing (conservative gradient check + symbolic integration)
- Validates momentum residual = 0
- Structural Poisson check
- Vorticity computation
- Pattern detection (rigid rotation, potential flow...)

All symbolic, exact, auditable.

## Key principle

Strict separation of facts vs assumptions.  
GeRoM never assumes silently. Every hypothesis (steady, energy_inequality, domain/BC) must be declared – or flagged as missing.  
No "phantom proofs". No intellectual noise.

## Pressure as consequence (in 2D)

GeRoM computes required ∇p, checks if conservative, reconstructs p up to constant.  
Pure classical differential calculus – executed cleanly.

## Output

Auditable reports: Why? Because of what? Under which assumptions? What blocks?  
Findings with cause + suggested fix.

## What GeRoM is NOT

- Not a thinking AI (no abstract reasoning)
- Not an NS solver (doesn't find solutions)
- Not automatic proof (doesn't replace human)

Deliberately limited. And that's its strength.

## Why it's strategic

Cleans noise (sign errors, implicit assumptions, inconsistent pressure).  
Turns NS into an auditable, versionable object.  
Prepares the ground for real reasoning (V7: dependency graphs).

## V7 coming soon – exclusive early testing on Patreon!

V7 adds DAG dependency graphs: causal paths, minimal proofs, blockages explained.  
Support on Patreon for early access + custom runs: https://patreon.com/AxelAbdelRahamaneMeghezzi

## Long-term vision: Geometric AI (~6 months)

Not an oracle. Role: heuristic guidance (propose ansätze/transformations) – ALWAYS validated by deterministic audit (V6.1/V7).  
AI proposes → Symbolic engine verifies → Only verified survives.

## Install & Run

```bash
pip install sympy networkx
python gerom_v6_1.py
