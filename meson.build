project('Hyph-m', 'fortran',
        default_options: ['buildtype=debugoptimized', 'fortran_std=f2018'],
        version : '1.0',
	license : 'AGPL-3.0-or-later')

add_project_arguments('-ffree-line-length-none', language: 'fortran')
# add_project_arguments( '-Wall', '-Wextra', '-Wunused-parameter', '-Wtabs', '-pedantic', '-Waliasing', language: 'fortran')

# linear algebra
dep_blas = dependency('blas', required: false)
if dep_blas.found()
  dep_blas_vendor = 'blas'
else
  message('Did not find the blas library. Trying with openblas')
  dep_blas = dependency('openblas', required: true) # required as it is the last try
  dep_blas_vendor = 'openblas'
endif

if dep_blas_vendor == 'openblas'
  dep_linalg = [dep_blas] # openblas already provides a lapack implementation
else
  dep_lapack = dependency('lapack', required: true)
  dep_linalg = [dep_blas, dep_lapack]
endif


dep_stdlib = dependency('fortran_stdlib')

dep_all = [dep_linalg, dep_stdlib]

subdir('src')
