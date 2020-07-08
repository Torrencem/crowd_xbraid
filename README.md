
This is the home REPO for the projects for RIPS 2020 (LLNL team). To get this code running:

    git clone https://github.com/Torrencem/crowd_xbraid.git --recursive
    
on Linux:

    sudo apt-get install liblapacke liblapacke-dev libopenblas-dev

or on Mac:

    brew install openblas
    brew install lapack

then simply:

    make

in the root directory to build. Please report any problems you have getting things to work.

If the build fails (probably on Mac) because of "-lgfortran", then try locating your `libgfortran.a` on your system (`sudo find /usr -iname 'libgfortran*.a'` might help), and symlinking it into the root directory of this project (`sudo ln -s /usr/<blahblah>/libgfortran.a ./libgfortran.a`).

Before committing code changes, please run:

    make fmt

To run the tests for the utility functions, run:

    make test
