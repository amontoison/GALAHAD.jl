# A [Julia](http://julialang.org) Interface to [GALAHAD](https://www.galahad.rl.ac.uk/)

## Custom Installation

To use your custom GALAHAD, set the environment variable `JULIA_GALAHAD_LIBRARY_PATH` to point to the shared library before `using GALAHAD`.

```bash
export JULIA_GALAHAD_LIBRARY_PATH=$(GALAHAD)/objects/pc64.lnx.gfo/double/shared/
```

The `JULIA_GALAHAD_LIBRARY_PATH` environment variable may be set permanently in the shell's startup file, or in `$HOME/.julia/config/startup.jl`.
