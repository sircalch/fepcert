# FEPCert validation

```bash
pip install alchemlyb alchemtest networkx
git worktree add ../fep110 83ff97a            # FEPCert 1.1.0, for the comparison
python validation/benchmark_alchemtest.py --legacy ../fep110   # GROMACS legs vs alchemlyb/pymbar
python validation/benchmark_cycles.py --legacy ../fep110       # FEP+ maps vs networkx
python validation/validate_alchemtest.py                       # earlier 1.1.0 check (benzene)
```

Results: `results/benchmark.csv` (one row per leg and variant), `results/cycles.csv` (one row per map).
