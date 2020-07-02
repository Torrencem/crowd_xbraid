
This is the home REPO for the projects for RIPS 2020 (LLNL team). To get this code running:

    git clone https://Torrencem/crowd_xbraid.git --recursive
    
on Linux:

    sudo apt-get install liblapacke liblapacke-dev libopenblas-dev

or on Mac:

    brew install openblas
    brew install lapack

then simply:

    make

in the root directory to build. Please report any problems you have getting things to work.

Before committing code changes, please run:

    make fmt
