# A04 — Results & Discussion / hypothesis framework (pre-specified)

Fill the [...] with actual numbers once `A04_downstream_DE.ipynb` has run. Pre-specifying
the predictions BEFORE seeing results keeps the interpretation honest (less HARKing).

## The hypothesis (inherited from A06)
A06 found a **non-canonical UPR**: ER chaperones (HSPA5, CALR, WFS1) induced *without*
canonical sensor activation (XBP1/ATF6/ATF4 flat), with a candidate mechanism of
**butyrate → HDAC inhibition → transcriptional chaperone induction**, and a parallel
**NF-κB / chemokine attenuation**. A06 could not isolate butyrate (Bt makes none; Cd does).
**A04 dials butyrate directly (0/5/10 mM), so it is the test.**

## Pre-specified predictions, per contrast

| Contrast (from notebook) | If butyrate→HDAC mechanism is real | If NOT (e.g. AR-driven, commensal-specific) |
|---|---|---|
| `butyrate_in_Mock` (Mock 10 vs 0) | HSPA5/CALR/WFS1 rise **dose-dependently**, XBP1/ATF6/ATF4 stay flat | chaperones flat or non-monotonic |
| `butyrate_in_Cd` (Cd 10 vs 0) | same chaperone induction, plus CXCL1/2/3/8 **fall** as butyrate rises | chemokines unchanged by butyrate |
| `infection_at_0mM` (Cd vs Mock, 0 mM) | strong inflammatory/chemokine induction (baseline Cd response) | weak/absent response (check infection worked) |
| `infection_at_5/10mM` | infection-induced chemokines **progressively attenuated** vs 0 mM | infection response unchanged across doses |

**Interaction (the headline):** the *difference* between `infection_at_0mM` and
`infection_at_10mM` effect sizes. Shrinking Cd response at higher butyrate = the A06
attenuation reproduced and now causally tied to butyrate.

## Decision tree for the discussion
1. Chaperones up with butyrate **and** canonical sensors flat → confirms non-canonical UPR;
   supports HDAC-inhibition transcriptional route. State it.
2. Chaperones up **and** CXCL down with dose **and** Cd response attenuated at 10 mM →
   butyrate is sufficient to reproduce the A06 protective signature → strong cross-model claim.
3. Chaperones up but canonical sensors **also** up → it's (at least partly) canonical UPR here;
   revise the A06 "non-canonical" framing for the organoid context.
4. No butyrate effect on chaperones/chemokines → butyrate is **not** the A06 driver; points back
   to AR-disinhibition / commensal-specific mechanism (MIST1/AR remain candidates).
5. MIST1 (BHLHA15) and AR (+TMPRSS2) behaviour across doses → note as supporting evidence for
   whichever regulator tracks the chaperones.

## Pathway / GSEA (downstream of DE, for the discussion)
- Run GSEA (hallmark + KEGG + reactome) or PROGENy on each contrast's ranked LFC.
- Expect (if hypothesis holds): UPR / "unfolded protein response" hallmark up with butyrate;
  TNFα-via-NFκB and inflammatory-response hallmarks **down** in Cd as butyrate rises.
- Same databases/approach as A06 → keep the two comparable.

## Results-table skeleton (fill from results/tables/DE_*.csv)
```
Contrast              | #DEG padj<0.05 | HSPA5 LFC | CALR | WFS1 | XBP1 | ATF4 | AR | BHLHA15 | CXCL8
butyrate_in_Mock      |  [...]         | [...]     | ...
butyrate_in_Cd        |  [...]         |
infection_at_0mM      |  [...]         |
infection_at_5mM      |  [...]         |
infection_at_10mM     |  [...]         |
```

## Cross-model framing (for the writeup)
Single-pathogen host response modulated by (a) a commensal partner [A06: Bt],
(b) a phage on the bacterial side [A03], (c) a metabolite [A04: butyrate]. A04's role is to
test whether the **metabolite alone** reproduces the commensal's protective UPR/chemokine
signature from A06.

## Caveats to state honestly
- Single donor (patient 35) → no biological-donor replication; effects are within-donor.
- First-look DE only; regulator/network inference (A06-style) is a later step.
- polyA host-only → no direct read on Cd transcriptional state; bacterial behaviour inferred.
- Model/database limits (OmniPath coverage) as in A06 — candidates, not proof.
