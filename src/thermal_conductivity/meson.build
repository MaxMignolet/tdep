executable('thermal_conductivity',
           'cumulative_kappa.f90',
           'kappa.f90',
           'main.f90',
           'options.f90',
           'scattering.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all)
