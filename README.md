# dmprofile

**Current release:** v1.0.0

This `Python` package provides intuitive methods that perform some of the most common steps needed to study dark matter simulations using the `pynbody` package. 

#### Installing
To use `dmprofile` you will need the `pynbody` package. Please follow [these instructions](https://pynbody.github.io/pynbody/installation.html). 
To install `dmprofile` you have two options:

* Clone the repository:

```bash
git clone https://github.com/b-fontana/DarkMatter.git
cd DarkMatter/
```
and try the examples, which are located in `DarkMatter/dmprofile/examples`. A data sample composed of subhalos and corresponding snapshot is provided in the `data/` folder, so that one can run the examples and scripts and understand the basic workflow.


* Pip install (installs the latest release, which does not necessarily correspond to the code in this repository):

```bash
pip install dmprofile
```

**Classes**

* `Halos()`: performs the management of the halos stored in the simulations. The class was built such that `pynbody` never has to be explicitly invoked. The following code sample:

```python
h = Halos('halos_path', 'snapshot_path', min_size=300)  
```
will load all the halos with at least 300 particles and corresponding subhalos stored in `halos_path` and the snapshots (optional) stored in `snapshot_path`. One can use the argument `N` instead if a fixed number of halos is preferred. Mixing both `N` and `min_size` is also possible. 

The instantiation will create backup halos and subhalos in case the user later decides to filter them, making the recover of the original halos possible.

All operations can be performed directly using the instance of the `Halos()` class. For instance:

```python
h.concentration_200(idx=idx, sub_idx=sub_idx)
h.get_shape(idx=idx, sub_idx=sub_idx)
``` 
where `idx` and `sub_idx` are the indexes stored internally by pynbody. They are ordered according to the size of the halos. Inside a specific halo, the subhalo having `sub_idx=0` will be the largest subhalo, the subhalo with `sub_idx=1` will be the second largest, and so on.

For a more comprehensive example, please see `script1.sh`, or check the examples in the `/dmprofile/examples/` folder. To run them, just do 

```bash
python dmprofile/examples/<name of the example>` 
```
from a directory placed above the package directory.

* `Plot()`: Takes care of all the plotting steps. Please see `script1_plot.py`.

This class is currently able to plot:
* Shape parameters
* Concentration
* Profiles (density, mass, ...)
* Relaxation
* Simple particle distributions

[Shape]: https://github.com/b-fontana/DarkMatter/blob/master/Shape.png
Example of a plot produced with `script1_plot.py`, displaying shape parameters of several simulations:
![Shape][Shape]
 

* Other utilities:
1. Parser that can be extended, stored in `dmprofile/src/parser`. To use it and print the parsed arguments:

```python
FLAGS, _ = add_args(argparse.ArgumentParser())
print("Parsed arguments:")
for k,v in FLAGS.__dict__.items():
    print('{}: {}'.format(k,v))
```

2. General utilities like memory consumption check, I/O operations, and others. Stored in `dmprofile/src/utilities.py`.

3. NFW_fit (`dmprofile/src/fit.py`) and centering function for the whole simulation and individual halos and subhalos (`dmprofile/src/move.py`). The latter file includes an additional function that is currently not operational.

**Performing customized operations**

In case the user wants to perform some operation that is still not included in the `dmprofile` package, it can retrieve the wanted halo or subhalo in the following way:

```python
h.get_halo(idx=idx)
h.get_subhalo(idx=idx, sub_idx=sub_idx)

h.get_halos()
h.get_subhalos()
```
where the last two line retrieve all the available halos and subhalos, respectively.


**Acknowledgements** 

*This package was developed during an internship at Swinburne University of Technology, Melbourne. I thank my supervisor, Dr. Alan R. Duffy, for having given all the required technical and theoretical support.*

Contact: *bfontanasantosalves@swin.edu.au*
