executable('phonon_dispersion_relations',
           'activity.f90',
           'densityplots.f90',
           'densityplots_stuntscattering.f90',
           'main.f90',
           'options.f90',
           'type_fast_interpolation.f90',
           'unfold_phonons.f90',
           'velocitydos.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all)
