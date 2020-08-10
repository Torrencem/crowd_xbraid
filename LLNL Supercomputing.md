# Sign-in and Infrastructure

## General
You have two usernames.  One is your LLNL Official User Name (OUN).  The other is your LC username (LCU).  These are not the same, although you use the same password with both.  Your OUN will probably have a number stuck on the end.

If you accidentally use your OUN instead of your LCU or vice versa, you will get no helpful error message.  It will disallow the login as if you got your password wrong.

You will also have an RSA key.  Whenever you are asked for a password, you type in your fixed passcode and then you type in the six-digit number shown on the RSA key.  This is true pretty much everywhere, whether online or on the cluster.

## Web
LLNL's LC web dashboard is available at [https://lc.llnl.gov/lorenz/mylc/mylc.cgi](https://lc.llnl.gov/lorenz/mylc/mylc.cgi).  To access it, you will have to get through possibly two logins.  The first will have a blank background, and you log in with your OUN.  If you're already logged in with your OUN, you might skip this.  The second login will say "Welcome to LC's Collaboration Zone!" in the top left, and will have a black-and-white photograph of a data center as the background.  **You log in here with your LCU.**

The dashboard has a number of widgets for your entertainment and education, but the most interesting is _job management_.  Access this via the tab bar at the top of the screen.  It lists your active jobs, letting you check whether your calculation has finished without logging into the computer.

## Cluster
LLNL's newest cluster is called Sierra.  We are not on Sierra.  We are on Quartz.  The URL for all things quartz-related is `quartz.llnl.gov`.  This is where you go for shell logins and for file transfers.

Login with `ssh LCU@quartz.llnl.gov`.  Remember to put your RSA key number in after your password.

I would recommend that, for file transfer purposes, you install an SFTP client.  My favorite is [FileZilla](https://filezilla-project.org/).  You can connect to `quartz.llnl.gov` via SFTP on port `22`.  (Do not save your passwords, as you will have to retype them every time to include the correct RSA key numbers.)

# Building code
The default `mpicc` compiler uses Intel's `icc` on the backend.  This doesn't work very well for our code.  Before you build anything, whether by running a compiler directly or by running `make`, load the GNU compilers with `module load gcc`.  You have to do this every time you log in, but only once per login.
