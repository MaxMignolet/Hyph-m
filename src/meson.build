# Copyright 2026 M. Mignolet
# SPDX-License-Identifier: AGPL-3.0-or-later

exec_name='hyph-m'

executable(exec_name,
           'main.F90',
           'commons.F90',
           'dynmat.F90',
           'lapack_wrapper.F90',
           'phonons.F90',
           'readInput.F90',
           'solve_saparov.F90',
           'solve_phmag.F90',
           dependencies: dep_all,
           install: true,
           install_dir: meson.project_build_root() / 'bin')

install_symlink(exec_name, pointing_to: meson.project_build_root() / 'src' / exec_name, install_dir: meson.project_source_root() / 'bin')
