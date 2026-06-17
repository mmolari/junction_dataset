# Extending the MGE annotation toolset

This document tracks the addition of three new detectors to the MGE sub-pipeline
(`rules/mges.smk`). They are implemented **one at a time, in order**: abricate →
Integron_Finder → anti-defense.

## The per-tool template

Every MGE tool already in the pipeline (geNomad, DefenseFinder, ISEScan) follows the same shape:

1. a per-accession **run rule** → `data/<tool>/{acc}/...`
2. an optional **location helper** (DefenseFinder only, to attach coordinates from the protein
   FASTA)
3. a **preformat script** under `scripts/annotations/` that reshapes the tool's native output to the
   unified MGE schema `id,iso,beg,end,type` → `results/mges/<tool>.csv`
4. registration in the `ZERO_BASED` dict + `MGE_TOOLS` list in `rules/mges.smk`.

`MGE_TOOLS = list(ZERO_BASED)`, so adding a key to `ZERO_BASED` auto-wires the tool into
`mge_assign_positions`, `mge_to_junction_gff`, and `mges_all`. The shared downstream
(`scripts/assign_junctions.py`, `scripts/mges_to_gff.py`) is tool-agnostic — adding a tool means
adding its env + run rule + preformat script and registering its name; **no downstream changes**.

`ZERO_BASED[tool]` is `True` when the tool reports 0-based element coordinates, `False` for 1-based
(it keys the `--zero_based` flag passed to `assign_junctions.py`).

## Tool 1 — abricate (resistance genes)

Not MGEs strictly, but valuable MGE cargo that tends to live in MGE-rich regions.

- **Env**: `config/conda_envs/abricate.yaml` (`abricate` from bioconda). abricate ships its
  databases in the package — no download rule.
- **Database**: `card`.
- **Run rule** `abricate_run` (per `{acc}`): `abricate --db card` over `data/fasta/{acc}.fa`,
  table to `data/abricate/{acc}/{acc}.tsv`.
- **Preformat** `scripts/annotations/abricate_df_preformat.py`: `iso = SEQUENCE`, `beg = START`,
  `end = END` (abricate is **1-based**), `type = "AMR"` (constant, mirroring `prophage`/`IS`/
  `defense_system`), `id = iso|GENE|index` (gene identity is preserved only here, since the unified
  schema has no gene column).
- **Register**: `ZERO_BASED["abricate"] = False`.

## Tool 2 — Integron_Finder (integrons) — DONE

- **Env**: `config/conda_envs/integron_finder.yaml` (`integron_finder` from bioconda). `--func-annot`
  uses the bundled `NCBIfam-AMRFinder.hmm` — no download rule.
- **Run rule** `integron_finder_run` (per `{acc}`): `integron_finder --local-max --func-annot --circ`
  over the isolate FASTA. Integron_Finder nests its output under `Results_Integron_Finder_<name>/`,
  so the rule flattens the `.integrons` (via the `Results_Integron_Finder_*` glob) to a stable
  `data/integron_finder/{acc}/{acc}.integrons`. SLURM resources bumped (cpus 4, 6h) — `--local-max`
  cmsearch is the slow step.
- **Preformat** `scripts/annotations/integron_finder_df_preformat.py`: skip the leading `#` comment
  line(s), drop rows with no `ID_integron` (no-integron placeholders), group per
  `(ID_replicon, ID_integron)`, `beg = min(pos_beg)`, `end = max(pos_end)` (**1-based**),
  `iso = ID_replicon`, `id = iso|n|type` (`n` parsed from `integron_01` -> `1`). **Decision:**
  `type` = the per-integron class column (`complete`/`CALIN`/`In0`), following the prior-analysis
  precedent (`data/junctions/annotations/loc/integronfinder.csv`) — these three classes are
  biologically distinct, so they are kept rather than flattened to a constant.
- **Register**: `ZERO_BASED["integron_finder"] = False`.

## Tool 3 — anti-defense systems

DefenseFinder with the `--antidefensefinder-only` flag — the **same binary**, so it reuses
`config/conda_envs/defensefinder.yaml` (no new env) and the existing
`defensefinder_models_download`.

- **Run rule** `antidefense_find`: mirrors `defensefinder_find` + `--antidefensefinder-only`,
  writing to a **separate** dir `data/antidefense_finder/{acc}/...` (the output filenames are
  identical to DefenseFinder's, so they must not share a directory).
- **Reuse** `defensefinder_gene_location.py` unchanged.
- **Reuse** `defensefinder_df_preformat.py`, parameterized with a new `--type` arg (default
  `defense_system`); anti-defense passes `--type antidefense_system` → `results/mges/antidefense.csv`.
- **Register**: `ZERO_BASED["antidefense"] = True` (0-based, mirroring DefenseFinder).
