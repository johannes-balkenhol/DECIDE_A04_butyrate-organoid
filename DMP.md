---
title: "DMP Template — VRC Würzburg (DFG-aligned)"
summary: "A fill-in Data Management Plan template with VRC infrastructure pre-filled. Generate the formal version in DataPLAN/RDMO; use this as the working reference."
tags: [DMP, template, RDM, VRC]
last_updated: 2026-06-07
permalink: dmp_template.html
sidebar: <FILL IN sidebar>
---

# DMP Template — VRC Würzburg

**How to use.** Generate the funder-formatted DMP in **DataPLAN** (`plan.nfdi4plants.org`, export to Word + JSON) or the **JMU RDMO** instance, and save it to `docs/DMP/` as a stable, versioned file. This page is the working reference: it shows the standard sections with the VRC answer already filled in, so you only edit the project-specific parts marked `<...>`. Review at the annual cycle; bump the version when scope changes.

## 1. Project information
- Project / acronym: `<...>` · Consortium / funding: `<e.g. DFG SFB 1583 (DECIDE)>`
- PI / Sprecher:in: `<...>` · Data steward: `<...>` · Start–end: `<...>`

## 2. Data summary — what data, what formats, what volume
- Data types: `<e.g. confocal microscopy, scRNA-seq, clinical questionnaires>`
- Raw formats: `<.czi/.lif, .fastq.gz, .csv>` → **primary (open) formats**: `<OME-TIFF, BAM/VCF, CSV/Parquet>`
- Estimated volume: raw `<GB>` / processed `<GB>`
- Re-use of existing data: `<yes/no — sources>`

## 3. Documentation & metadata
- Metadata standard (by domain): `<MIxS / OME-XML+REMBI / FHIR+SNOMED-LOINC / RRID+ChEBI>`
- Template source: CoreUnitRDM MetaDataTemplates + project `data/metadata/data_dictionary.md`
- Controlled vocabularies via EBI OLS; missing values coded with standard tokens (not blanks)
- Provenance: README + Labfolder ELN + VRC GitLab commit history

## 4. Storage & backup (active phase)
- Active storage: **cRDM Nextcloud group folder** (`coreunitrdm.biozentrum.uni-wuerzburg.de`)
- Backup: nightly to RZ; **3-2-1 rule** (3 copies, 2 media, 1 off-site) for raw
- Compute working data: HPC scratch (not backed up) — results synced back to Nextcloud
- Raw data: read-only + checksummed; never edited

## 5. Legal & ethics
- Personal/patient data? `<yes/no>`. If yes:
  - Ethics approval (Ethikkommission JMU/UKW): `<ref / pending>`
  - Informed consent covering planned use + reuse/publication: `<yes/no>`
  - DPIA for high-risk processing + entry in the VVT: `<status>`
  - Pseudonymisation key held only at UKW, separate from data; access logged
  - **Patient-identifying data stays at UKW — not in the cRDM VRC**
- IP / third-party restrictions: `<...>`

## 6. Preservation & sharing
- Long-term archive: **RZ Archive Server**, 10-year retention
- Repository at publication: `<WueData / Zenodo / ENA / EGA / IDR / PRIDE / NFDI4Health Study Hub>`
- Access level: `<open / controlled-access under DAA>`; embargo: `<none / until publication>`
- Licence: data `<CC BY 4.0 / CC0>`, code `<MIT / Apache-2.0 / GPL-3.0>`
- DOI: assigned at Gate G2; referenced in the manuscript

## 7. Responsibilities & resources
- Data steward / who maintains the DMP: `<...>`
- Costs (storage tier, archive, repository fees): `<...>`
- cRDM service package & start date: `<...>`

*Living document — review annually (Service Lifecycle SOP, Phase 3) and on any scope change. Related: SOP — RDM Lifecycle · Metadata Template Suggestions.*
