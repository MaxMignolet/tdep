executable('anharmonic_free_energy',
           'energy.f90',
           'epot.f90',
           'main.f90',
           'options.f90',
           link_with: [libolle, libflap],
           include_directories: ['../libolle', '../libflap'],
           dependencies: dep_all)
