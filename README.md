# shape-it-ob3


## Description

Code for shape-it with openbabel3


## INSTALL

- Following example is the basic way to install the tool.

```
git clone https://github.com/iwatobipen/shape-it-ob3.git
cd shape-it-ob3
mkdir biuild
cd build
cmake -DCMAKE_INSTALL_PREFIX=<where you want to insatall> ..
make
make install
```


- If you would like to use rdkit for shape-it please set BUILD_RDKIT_SUPPORT=ON. Also set BUID_PYTON_SUPPORT to ON, you can call shape-it as a library.


```
cd shape-it-ob3
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=<where you want to install> -DBUILD_RDKIT_SUPPORT=ON -DBUILD_PYTHON_SUPPORT=ON ..
make
make install
```

## Tips
- If make command failed, you shoud try to set RDKIT_INCLUDE_DIR, Boost_INCLUDE_DIR option.

```

cmake -DCMAKE_INSTALL_PREFIX=/home/iwatobipen/src/shape-it-ob3 -DRDKIT_INCLUDE_DIR=/home/iwatobipen/miniconda3/envs/chemoinfo/include/rdkit -DBUILD_RDKIT_SUPPORT=ON -DBUILD_PYTHON_SUPPORT=ON -DBoost_INCLUDE_DIR=/home/iwatobipen/miniconda3/envs/chemoinfo/pkgs/libboost-1.73.0-hf484d3e_11/include ..
make
make install
```


## Original code
- http://silicos-it.be.s3-website-eu-west-1.amazonaws.com/software/shape-it/1.0.1/shape-it.html



## History of Shape-it and how to cite

Shape-it is a rewrite of the original Pharao code that was developed in 2008 by Silicos (Jonatan Taminau, Gert Thijs and [Hans De Winter](https://github.com/hansdewinter)). If you use this code in your research, we would appreciate if you would include the following citation in your publication:

Taminau, J.; Thijs, G.; De Winter, H. (2008) ‘Pharao: Pharmacophore alignment and optimization’, J. Mol. Graph. Model. 27, 161-169

