# Bolt
![Build Status](https://github.com/sakkayaphab/bolt/workflows/Ubuntu/badge.svg?branch=master)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/sakkayaphab/bolt)](https://github.com/sakkayaphab/bolt/releases)




## Installation

### 1. Installing Bolt from source

Requirements
1. GCC >=v5.4
2. Intel tbb-devel >= 2019.9
3. Intel tbb-devel >= 2019.9
4. HTSlib >= 1.9


```sh
git clone https://github.com/sakkayaphab/bolt.git
mkdir build
cd build
cmake .. -DINSTALL_BIN_PREFIX=(directory for installation) -DINCLUDE_LIBRARY_PREFIX=(directory for include) -DLIBRARY_LINK_PREFIX=(directory for lib)
make
make install
```


### 2. Installing Bolt with conda


```sh
conda install -c conda-forge tbb
conda install -c bioconda bolt
```


## Usage
```sh
bolt call -b (aligment) -r (reference) -o (output)
```
for [user guide][UserGuide]


[UserGuide]:docs/README.md

## License
Bolt is distributed under the [MIT License][MITLicense]. However, Bolt includes several third-party open-source libraries, please see [Third Party Software Notices][LICENSETHIRDPARTY] for details.


[MITLicense]:LICENSE
[LICENSETHIRDPARTY]:THIRD-PARTY-LICENSE