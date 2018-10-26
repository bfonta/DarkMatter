# dmprofile

**Current release:** v1.0.0 [Package released!]

This `Python` package provides intuitive methods that perform some of the most common steps needed to study dark matter simulations using the `pynbody` package. In order to install `pynbody` please follow [these instructions](https://pynbody.github.io/pynbody/installation.html). A data sample is provided in the `data/` folder, so that one can run the examples and scripts and understand the basic workflow.

The package includes the following classes: 

*`Halos()`: performs the management of the halos stored in the simulations. The class was built such that `pynbody` never has to be explicitly invoked. The following code sample:

```python
h = Halos('halos_path', 'snapshot_path', min_size=300)  
```
will load all the halos with at least 300 particles and corresponding subhalos stored in `halos_path` and the snapshots (optional) stored in `snapshot_path`. One can use the argument `N` instead if a fixed number of halos is preferred. Mixing both `N` and `min_size` is also possible. The instantiation will also create backup halos and subhalos in case the user later decides to filter them. The previous halo can then be recovered.
All operations can be performed directly using the instance of the `Halos()` class. For instance:

```python
h.concentration_200(idx=idx, sub_idx=sub_idx)
h.get_shape(idx=idx, sub_idx=sub_idx)
``` 
where `idx` and `sub_idx` are the indexes stored internally by pynbody which are ordered based on halo size. Inside a specific halo, the subhalo with `sub_idx=0` will be the largest subhalo, the subhalo with `sub_idx=1` the second largest, and so on.
In case the user wants to perform some operation that is still not include in the `dmprofile` package, it can retrieve the wanted halo or subhalo in the following way:

```python
h.get_halo(idx=idx)
h.get_subhalo(idx=idx, sub_idx=sub_idx)

h.get_halos()
h.get_subhalos()
```
where the last two line retrieve all the available halos and subhalos, respectively. For a more comprehensive example, please see `script1.sh`, or check the examples in the `/dmprofile/examples/` folder. To run them, just do `python dmprofile/examples/<name of the example>` from a directory placed above the package directory.

* `Plot()`: Takes care of all the plotting steps. Please see `script1_plot.py`.

This class is currently able to plot:
*Shape parameters
*Concentration
*Profiles (density, mass, ...)
*Relaxation
*Simple particle distributions

[Shape]: https://github.com/b-fontana/DarkMatter/blob/master/Shape.png
![Shape][Shape]
 

*Other utilities:
**a
**b
**c

*This package is part of an internship at Swinburne University of Technology, Melbourne.*
