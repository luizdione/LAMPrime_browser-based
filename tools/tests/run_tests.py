#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_tests.py - dependency-light runner for the LAMPrime reproducibility suite.
Discovers every test_*.py in this directory, runs each test_* function, and
prints PASS/FAIL. Works without pytest (pytest also collects these files).

Usage:  python tools/tests/run_tests.py
Needs:  Python 3.10+, Node.js 18+, biopython==1.86 (pip install -r requirements.txt)
"""
import os, sys, importlib, traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")


def main():
    mods = sorted(f[:-3] for f in os.listdir(HERE)
                  if f.startswith("test_") and f.endswith(".py"))
    total = passed = 0
    failures = []
    for m in mods:
        mod = importlib.import_module(m)
        for name in sorted(dir(mod)):
            if not name.startswith("test_"):
                continue
            fn = getattr(mod, name)
            if not callable(fn):
                continue
            total += 1
            try:
                fn()
                passed += 1
                print(f"  PASS  {m}.{name}")
            except Exception as e:
                failures.append((m, name, e))
                print(f"  FAIL  {m}.{name}: {e}")
                traceback.print_exc()
    print("\n" + "=" * 60)
    print(f"LAMPrime reproducibility suite: {passed}/{total} passed")
    print("=" * 60)
    if failures:
        for m, name, e in failures:
            print(f"  - {m}.{name}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
