
This is the home REPO for the projects for RIPS 2020 (LLNL team). To get this code running:

    git clone https://github.com/Torrencem/crowd_xbraid.git --recursive

Then use ``make`` to build the example code you want to run (use ``make list`` to see the available targets) in the root directory.

Please report any problems you have getting things to work.

If the clone or build fails because of weirdness with the submodule (very annoying), do ``cd xbraid`` and ``git checkout trimgrit_cxx``.

Before committing code changes, please run:

    make fmt
