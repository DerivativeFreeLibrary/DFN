-----------------------------------------------------------------------------------
 How to use the derivative-free optimizers CS-DFN and DFN_simple for Nonsmooth NLP
-----------------------------------------------------------------------------------

1- Gunzip and untar the archive in a folder on your computer.

2- Edit file problem.f90 to define your own MINLP problem.
   Be aware that an equality constraint h(x) = 0
   must be written as h(x) <= 0 and -h(x) <= 0.

3- At command prompt execute 

     $> make
 
   which will create the executables 'cs-dfn' and 'dfn_simple'

4- execute

     $> ./cs-dfn
     $> ./dfn_simple
