The data is strutured as:

`data[i]=["MODEL",{params : params_val} , μ , np.array([ρ_avg]),np.array([ρ_Var]) ) `

The parameter part is a dictionary and could be any (noise, rewiring p, m for BA ecc.);
μ is the graph realization index.

Those are .pkl files, to load them correcly run

```
results_SWSF = []
filenames = []    

for pattern in [
    "SWFS_m3_p_0_N4000_chunk_*.pkl",
    "SWFS_m3_p_sweep_N4000_chunk_*.pkl",

]:
    for fn in sorted(glob.glob(pattern)):
        filenames.append(fn)
        with open(fn, "rb") as f:
            results_SWSF.extend(pickle.load(f))
```
if you want to merge multiple filenames into one object, or
```
chunk_files = sorted(glob.glob("SW_m3_p_sweep_N4000_chunk_*.pkl"))
results_SW = []
for fn in chunk_files:
    with open(fn, "rb") as f:
        results_SW.extend(pickle.load(f))
```
for just one filename.

The filenames are just examples, but they are present in the repository.
