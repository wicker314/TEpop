# TEpop
transposable element population analysis

+++++ WICKERsoft set up +++++++++++++++++ 

- all scripts depend on subroutines stored in the perl module Biomodule.pm. 

- place the Biomodule.pm and all WICKERsoft scripts into a dedicated directory which is in a common path (e.g. /usr/local/bin/) so that you can run them from anywhere in your system.

- change the path in the beginning of all WICKERsoft scripts to that directory

example, how the header should look:
#!/usr/bin/perl

use lib '/usr/local/bin/';
use Biomodule;

- make sure that all scripts are executable. if necessary, make them executable with chmod 755 <script_name>

- all WICKERsoft scripts have usage instructions to tell you what type of input files and parameters are needed. simply type the script name without arguments to see usage instructions



Comments:
- Wild cards (* or ?): WICKERsoft scripts do not recognize wild cards. So, if you are asked to provide a query for multiple files, provide a string common to all file names you want to process without wild cards. 

- some WICKERsoft scripts use the perl-tk module for graphics. Thes have the additional line 'use Tk;' in the header. The perl-tk module can be obtained e.g. from Ubuntu repositories. 
