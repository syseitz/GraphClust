# MMseqs2 Prefilter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the legacy BlastClust preprocessing backend in GraphClust with MMseqs2 while keeping the existing near-identical fragment filtering behavior as close as possible.

**Architecture:** Keep the current preprocessing control flow intact and swap only the clustering backend in `graphFasta.pl`. Introduce MMseqs2-named config keys and fall back from old BlastClust keys to preserve compatibility during migration. Verify behavior first on a small regression dataset, then with an end-to-end pipeline run that uses the local cmfinder checkout.

**Tech Stack:** Perl, shell regression tests, MMseqs2, existing GraphClust pipeline scripts

---
