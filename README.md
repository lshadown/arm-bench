# arm-bench

## Installation

```bash
brew install llvm
brew install libomp
```


```bash
-DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang
-DCMAKE_CXX_COMPILER=/opt/homebrew/opt/llvm/bin/clang++
-DOpenMP_C_FLAGS="-fopenmp"
-DOpenMP_C_LIB_NAMES="omp"
-DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib
```



