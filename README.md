# dmprofile
## Simplify your pynbody dark matter halos analysis

**Current release**

* v0.1.0: Not meant to be used. Further testing required.



This package provides intuitive methods that perform some of the most common steps needed to study `pynbody` dark matter simulations. In order to install `pynbody` please follow [these instructions](https://pynbody.github.io/pynbody/installation.html).

The package currently includes two major classes, `Halos()` and `Plot()` plus some additional convenience functions. The halos are managed with the first class while all plots are done with the second one. A simple example follows:

 
```python
import pynbody
import dmprofile
from dmprofile.move import centering_com

h = Halos("simulation_file", max_number_of_halos)

prof1 = h.get_profile(halo_idx, 'dm', bins=(10,45,BIN_NUMBER), bin_type='linear', normalize=False)
prof2 = h.get_profile(halo_idx, 'dm', bins=(10,45,BIN_NUMBER), bin_type='linear', normalize=False)

p = plot.Profile([prof1,prof2], h=9, w=8, name="Density.png")
p.set_all_properties(model='density_profile', xscale='log', yscale='log')
p.plot_all("radius", "density")
p.fit_and_plot_all('nfw')                                                                                     
p.savefig()
```

which produces the following:

[Density]: https://github.com/b-fontana/DarkMatter/blob/master/Density.png
![Density][Density]
 
*This package was created during an internship at Swinburne University of Technology, Melbourne.*
