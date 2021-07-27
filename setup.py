from skbuild import setup

setup(
    name="NumbaLSODA",
    packages=['NumbaLSODA'],
    version='0.1.0',
    license='MIT',
    install_requires=['numpy','numba'],
    author = 'Nicholas Wogan',
    author_email = 'nicholaswogan@gmail.com',
    description = 'Python wrapper of LSODA (solving ODEs)'+\
                  ' which can be called from within numba functions.',
    cmake_args=['-DSKBUILD=ON']
    )