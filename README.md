# Resistance evolution in a single species reorganizes community responses across a resource-antibiotic landscape

Data and analysis code for the study by Inga-Katariina Aapalampi, Anna M. Bischofberger,
Dovydas Kičiatovas, Mikko Kivikoski, Ville Mustonen, Teppo Hiltunen and Johannes Cairns.

Three- and 16-species bacterial communities containing either ancestral or
streptomycin-resistance-evolved *Aeromonas caviae* HAMBI 1972 were grown across 16
glucose-derived carbon concentrations (0.07–2.00 g C L⁻¹) crossed with four streptomycin
concentrations (0, 2, 20 and 200 µg ml⁻¹), with eight replicate communities per community type.
Endpoint composition was measured by 16S rRNA gene amplicon sequencing. Monoculture growth of
the three core species was measured across the same 64 environments.

The notebooks reproduce every statistic, table and figure of the article and its Supplementary
Information from the files in `data/`.

## Repository layout

```
data/                       input data (described below)
notebooks/                  analyses, run in numerical order
  00_setup_and_functions    paths, constants, data loading, statistical tests, segmentation,
                            compositional analyses and plotting helpers; run by notebooks 01-05
  01-05                     analyses (see Analyses below)
results/figures/            Figures 2–5 and S1–S12 (PDF and PNG)
results/tables/             Supplementary Tables S1–S14 (CSV)
results/statistics/         further statistics reported in the text (CSV)
```

## Data

All files are tab-separated. Streptomycin is in µg ml⁻¹ and carbon in g of glucose-derived
carbon per litre.

**`data/community_composition.tsv`** — endpoint composition of the 2,048 experimental
communities (4 community types × 8 replicates × 64 environments), one row per sample and species.

| Column | Description |
|---|---|
| `sample` | sample name, `C<1–4>_Gl<carbon level>_St<streptomycin>_R<replicate>`, with C1–C4 = `ANC_THREE`, `EVO_THREE`, `ANC_SIXTEEN`, `EVO_SIXTEEN` |
| `community` | `ANC_THREE`, `EVO_THREE`, `ANC_SIXTEEN`, `EVO_SIXTEEN`: three- or 16-species community containing ancestral (ANC) or resistance-evolved (EVO) *A. caviae* |
| `replicate` | replicate community (1–8); a replicate is one constructed landscape |
| `carbon_level`, `carbon_g_L` | carbon level (1–16) and concentration |
| `streptomycin` | streptomycin concentration |
| `species` | HAMBI strain code (see `species.tsv`) |
| `read_count` | reads assigned to the species by Rbec |
| `copy_number_corrected_count` | read count divided by the strain's 16S rRNA gene copy number, rounded |
| `relative_abundance` | copy-number-corrected relative abundance among the 16 species |

**`data/species.tsv`** — species names, membership of the three-species community and
streptomycin MIC class (µg ml⁻¹) of each species in the ancestral and resistance-evolved
communities (only *A. caviae* differs).

**`data/monoculture_growth.tsv`** — growth traits of each strain in monoculture across the 64
environments, six replicates per strain and environment, from the third 48-h growth cycle.

| Column | Description |
|---|---|
| `strain` | `HAMBI_1287` (*Citrobacter koseri*), `HAMBI_1972_ANC` / `HAMBI_1972_EVO` (ancestral / resistance-evolved *A. caviae*), `HAMBI_1977` (*Pseudomonas chlororaphis*) |
| `carbon_level`, `carbon_g_L`, `streptomycin` | environment |
| `replicate`, `plate`, `well` | replicate (1–6) and its position |
| `lag_time_h` | lag time (h) |
| `max_growth_rate` | maximum growth rate (OD₆₀₀ min⁻¹) |
| `yield` | yield (OD₆₀₀) |
| `auc` | area under the growth curve (OD₆₀₀ h) |

Raw sequencing reads are in the European Nucleotide Archive under accession **PRJEB129114**
(sample aliases `GLSTR_<sample>`, with `R<n>` written as `rep<n>`). Reads were processed to
species counts with the hambiAmplicon workflow (gitlab.utu.fi/slhogl/hambiAmplicon).

## Analyses

| Notebook | Content | Outputs |
|---|---|---|
| `01_community_states.ipynb` | segmentation of each landscape into community states; boundary configuration, region dependence and orientation; number and size of states; dominant species at inhibitory doses; boundary change versus focal-species sensitivity | Figures 2, S1–S4; Table S1 |
| `02_alpha_diversity.ipynb` | Shannon Index responses | Figures 3, S5, S6; Tables S2, S3 |
| `03_low_dose_response.ipynb` | antibiotic responses in log-ratio space; focal-species response at 2 µg ml⁻¹ along carbon; focal and non-focal responses; species-level contrasts | Figures 4, S7; Tables S4–S7, S13, S14 |
| `04_resistance_filtering.ipynb` | enrichment of high-MIC taxa among non-focal species; MIC class and species response | Figures S8, S9; Tables S8, S9 |
| `05_growth_traits.ipynb` | monoculture growth traits and competitive position in three-species communities | Figures 5, S10–S12; Tables S10–S12 |

`results/statistics/` holds the statistics that are reported in the text but not in a
Supplementary Table: dominant species at 20 and 200 µg ml⁻¹ (`dominance_at_inhibitory_doses.csv`),
boundary change versus focal-species sensitivity (`boundary_change_vs_focal_sensitivity.csv`),
species detection frequencies (`species_detection_16_species.csv`), the whole-community version
of the species contrasts (`lowdose_species_contrasts_whole_community.csv`) and the pseudocount
sensitivity of log-ratio contrasts (`clr_pseudocount_sensitivity.csv`).

### Statistical approach

The inferential unit is the constructed landscape: one replicate community measured across all
64 environments. Histories are compared by exact randomization over all 12,870 allocations of
the 16 landscapes of a richness level to the two histories; one-sample tests use all 256 sign
flips. Confidence intervals are percentile bootstrap intervals from resampling landscapes within
history. Random seeds are fixed, so every result is reproduced exactly.

Community states are found by entropy-penalized spatial segmentation, developed by Dovydas
Kičiatovas: adjacent environments are merged greedily while the increase in
within-state compositional heterogeneity (a generalized Jensen–Shannon loss) is smaller than the
penalty for the boundary removed, subject to a maximum Bray–Curtis diameter of a state. The
method is described in `notebooks/00_setup_and_functions.ipynb`. Generative AI tools were used
to assist in developing the algorithm and writing the code.

## Running the analysis

Python 3.13 with the packages in `requirements.txt`:

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cd notebooks
for nb in 0*.ipynb; do jupyter nbconvert --to notebook --execute --inplace "$nb"; done
```

The notebooks run in about two minutes in total and overwrite the files in `results/`. Each of
notebooks 01-05 first runs `00_setup_and_functions.ipynb`, so any of them can also be run on its own.

## Use of AI

Generative AI tools were used to assist with writing the code. All code was reviewed and approved
by the authors.

## Licence

Code: MIT licence (`LICENSE`). Data in `data/` and results in `results/`: Creative Commons
Attribution 4.0 International (`LICENSE-DATA`).
