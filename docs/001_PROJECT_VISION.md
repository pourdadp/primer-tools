# File: D:\projects\primer-tools-main\docs\001_PROJECT_VISION.md

# Project Vision

Version: 0.1

Status: Draft

Last Updated: 2026-07-03

---

# QuickNGS Vision

## Mission

QuickNGS is an intelligent open-source platform for Sanger sequence analysis.

Our goal is to transform raw sequencing files into reliable biological conclusions with the minimum user interaction.

---

# Philosophy

Biology drives software.

Software should adapt to biological workflows, not force users to adapt to software limitations.

---

# Design Principles

## 1. Auto First

Auto Mode is the default and recommended workflow.

The software should automatically determine the best analysis strategy.

Users should obtain meaningful results with a single click whenever possible.

---

## 2. Expert Always

Every automatic decision can be reviewed and overridden.

Advanced users always have complete control over analysis parameters.

---

## 3. One Click Analysis

A new user should be able to analyze sequencing files within 30 seconds.

Complexity should remain hidden unless requested.

---

## 4. Explain Every Decision

Every automatic decision must include:

- Reason
- Evidence
- Confidence

No black-box analysis.

---

## 5. Scientific Transparency

All biological conclusions must be reproducible.

Every report should clearly describe how results were generated.

---

## 6. Modular Architecture

Every major function must be implemented as an independent module.

Modules should communicate through clearly defined interfaces.

---

## 7. User-Centered Design

Software decisions should prioritize the needs of the majority of users.

Features useful only for rare situations should be implemented as optional modules rather than becoming part of the Core.

---

# Target Users

QuickNGS is designed for:

- Clinical laboratories
- Veterinary laboratories
- Research laboratories
- Universities
- Students
- Molecular biology laboratories
- Microbiology laboratories
- Virology laboratories

---

# Auto Mode

Auto Mode automatically performs analysis using the most appropriate workflow.

Examples include:

- File validation
- Sequence quality assessment
- Quality trimming
- Assembly
- Consensus generation
- Reference identification
- BLAST confirmation
- Serotype analysis (when applicable)
- Final report generation

---

# Expert Mode

Expert Mode provides access to all configurable parameters without changing the Core workflow.

---

# Project Goal

Build the most user-friendly open-source platform for Sanger sequence analysis.

The software should help users reach biological conclusions rather than simply generating sequence alignments.

---

# Core Values

- Simplicity
- Scientific accuracy
- Transparency
- Reliability
- Automation
- Extensibility
- Reproducibility

---

# Guiding Question

Before implementing any feature, ask:

"Will this improve the experience of most users?"

If the answer is No, the feature should probably become an optional plugin rather than part of the Core.

---

End of Document