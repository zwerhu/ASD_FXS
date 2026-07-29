# NumPy >=2.0 compatibility patches for MetaXcan / S-PrediXcan

The reference [MetaXcan](https://github.com/hakyimlab/MetaXcan) implementation of
S-PrediXcan predates NumPy 2.0 and fails in two places when run under a modern
NumPy (>=2.0) environment:

1. **Legacy `numpy.matrix` covariance objects.** The covariance-loading code
   returns `numpy.matrix` objects rather than plain `ndarray`s. Dot products
   of a vector against a `numpy.matrix` always return a 2-D (1x1) matrix
   rather than a scalar. NumPy < 2.0 silently coerced such single-element
   arrays to Python floats; NumPy >= 2.0 raises
   `TypeError: only 0-dimensional arrays can be converted to Python scalars`
   in `AssociationCalculation.association()`.

2. **Degenerate single-SNP covariance.** For genes whose prediction model
   contains only one SNP, the "covariance matrix" degenerates to a
   zero-dimensional array. An un-guarded diagnostic eigenvalue check
   (`numpy.linalg.eig(cov)`) that runs at verbose logging levels crashes on
   this shape, terminating the entire S-PrediXcan run for that
   tissue/GWAS pair with no output file produced and no partial results
   saved.


## How to apply

```bash
git clone https://github.com/hakyimlab/MetaXcan.git
cd MetaXcan/software
patch -p2 < /path/to/patches/numpy2_compatibility.patch
```

If the patch context does not apply cleanly against a newer MetaXcan
revision, the fix can be applied by hand: in
`metax/metaxcan/AssociationCalculation.py`, function `association()`:

- Wrap `cov` in `numpy.asarray(...)` before any linear-algebra operation
  (`numpy.linalg.eig`, `numpy.dot`), and guard the singularity-check block so
  it only runs when the resulting array is genuinely 2-D and square.
- Extract the final `sigma_g_2` scalar via
  `float(numpy.asarray(...).reshape(-1)[0])` instead of a bare `float(...)`
  call on the raw dot-product result.

Patched code produces numerically identical results to the reference
implementation for every case in which the reference implementation executed
successfully; the patch only prevents crashes, it does not change any
formula.
