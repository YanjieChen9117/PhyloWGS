# Installing Python 2.7 for PhyloWGS

PhyloWGS requires **Python 2**. macOS no longer ships Python 2, and Homebrew has removed the `python@2` formula. Use **pyenv** to install and run Python 2.7.

## 1. Install pyenv (if needed)

```bash
brew install pyenv
```

Add to your shell config (e.g. `~/.zshrc`):

```bash
eval "$(pyenv init -)"
```

Then run `source ~/.zshrc` or open a new terminal.

## 2. Install Python 2.7.18

```bash
# Optional: install build dependencies (if install fails)
xcode-select --install
brew install openssl readline sqlite3 xz tcl-tk zlib pkg-config

# Install Python 2.7
pyenv install 2.7.18
```

## 3. Make `python2` available for PhyloWGS

**Option A – Use in this project only (recommended)**  
In the PhyloWGS repo root:

```bash
cd /Users/khanhngocdinh/Documents/YanjieChen/GitHub/PhyloWGS
pyenv local 2.7.18
```

Then ensure your R process sees this: run R from the same terminal (or from a shell where `pyenv` is in PATH), so that `Sys.which("python2.7")` finds the pyenv shim. With `pyenv local 2.7.18`, `python` and `python2.7` in that directory will point to 2.7.18.

Create a `python2` symlink so the script finds it:

```bash
# From PhyloWGS repo root, after pyenv local 2.7.18
ln -sf $(pyenv which python2.7) ./python2
export PATH="$PWD:$PATH"
```

Or run R with PATH set so `python2` is the pyenv interpreter:

```bash
export PATH="$(pyenv root)/shims:$PATH"
Rscript PhyloWGS_runtime_analysis/run_phylowgs_dream_runtime.R
```

**Option B – Global Python 2.7**  
Use Python 2.7 as the default for your user (affects all projects):

```bash
pyenv global 2.7.18
```

Then ensure `~/.pyenv/shims` is in your PATH (pyenv init does this). The shims provide `python2.7`; there is usually no `python2` by default. Either:

- Use the script’s fallback: it will use `python2.7` if `python2` is not found, or  
- Add a shim/symlink for `python2` that runs 2.7.18.

## 4. PhyloWGS Python 2 dependencies

With Python 2.7 active (via pyenv), install:

```bash
pip install --user numpy scipy ete2
```

**GSL** (GNU Scientific Library) is required to compile the C++ part:

```bash
brew install gsl
```

Then in the PhyloWGS repo root:

```bash
cd /Users/khanhngocdinh/Documents/YanjieChen/GitHub/PhyloWGS
g++ -o mh.o -O3 mh.cpp util.cpp $(gsl-config --cflags --libs)
```

## 5. scipy on Apple Silicon (arm64)

On **Apple Silicon Macs**, `pip install scipy` for Python 2.7 often **fails** to build (SuperLU/ARPACK C and Fortran issues with modern compilers). The following are already done on this machine where possible:

- **pyenv** and **Python 2.7.18** are installed.
- **numpy 1.16.6** is installed (from patched source to avoid `-faltivec` on arm64).
- **ete2** and **GSL** are installed; **mh.o** is compiled.
- **scipy** could not be built from source (C99/Fortran errors).

PhyloWGS **requires scipy** (e.g. `scipy.stats`, `scipy.special`). Options:

1. **Run PhyloWGS in a Linux environment** (e.g. Docker or a Linux VM) where `conda create -n py27 python=2.7 numpy scipy` works.
2. **Use an Intel Mac or x86_64 Linux** to generate the conda/pip environment, then copy the `site-packages` (or use a shared env).
3. If you obtain a pre-built **scipy wheel** for `macosx_arm64` + Python 2.7 (e.g. from a third-party or built on another setup), install it with `pip install that_wheel.whl`.

## 6. Check

```bash
cd /Users/khanhngocdinh/Documents/YanjieChen/GitHub/PhyloWGS
eval "$(pyenv init -)"
python2.7 --version   # Should print: Python 2.7.18
python -c "import numpy, ete2; print('numpy, ete2 OK')"
# python -c "import scipy; print('scipy OK')"   # only if scipy is installed
```

Run the R runtime script from a **terminal** where pyenv is in PATH (so `python2.7` is found):

```bash
export PATH="$(pyenv root)/shims:$PATH"
Rscript PhyloWGS_runtime_analysis/run_phylowgs_dream_runtime.R
```
