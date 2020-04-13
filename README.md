# Bolt
![Ubuntu](https://github.com/sakkayaphab/bolt/workflows/Ubuntu/badge.svg?branch=master)
![MacOS](https://github.com/sakkayaphab/bolt/workflows/MacOS/badge.svg?branch=master)
![Docker](https://github.com/sakkayaphab/bolt/workflows/Docker/badge.svg?branch=master)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/sakkayaphab/bolt)](https://github.com/sakkayaphab/bolt/releases)
[![Conda](https://img.shields.io/conda/v/bioconda/bolt?color=blue&label=Anaconda%20Cloud)](https://anaconda.org/bioconda/bolt)



## Installation


### 1. Installing Bolt from source

Requirements
1. GCC >=v5.4
2. HTSlib >=v1.9
3. CMake >=v3.9
4. Intel TBB >=v2019.9

```sh
git clone https://github.com/sakkayaphab/bolt.git
mkdir build
cd build
cmake .. -DINSTALL_BIN_PREFIX=(directory path for the binary executable file of Bolt) -DINCLUDE_LIBRARY_PREFIX=(directory path for include of libraries) -DLIBRARY_LINK_PREFIX=(directory path for lib of libraries)
make
make install
```


### 2. Installing Bolt with conda

```sh
conda install -c intel tbb
conda install -c bioconda bolt
```

### 3. Run with docker
```
docker run -v /hostpath/:/opt/mount --name bolt sakkayaphab/bolt:latest /project/build/bolt call -b /opt/mount/sample.bam -r /opt/mount/ref.fa -t 100 -o /opt/mount/outputfolder
```

## Usage
```sh
bolt call -b (aligment) -r (reference) -t (number of threads) -o (output)
```


for [user guide][UserGuide]


[UserGuide]:docs/README.md

## License
Bolt is distributed under the [GPLv3][GPLv3]. However, Bolt includes several third-party open-source libraries, please see [Third Party Software Notices][LICENSETHIRDPARTY] for details.


[GPLv3]:LICENSE
[LICENSETHIRDPARTY]:THIRD-PARTY-LICENSE