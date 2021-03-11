# parallel_pwmig
This package is a direct is a child of a more prototype package of research code called just "pwmig" that at this writing is 
found in same area.  i.e. pavlis/pwmig as compared this package pavlis/parallel_pwmig.  As the prepended "parallel" implies this 
version is a rewrite to work with the MsPASS (Massively Parallel Analysis System for Seismology) framework.  When finished this 
package will be an add on to MsPASS to accomplish a more specialized task than MsPASS.  That is, this package will allow a MsPASS 
user to assemble a data set through the stage of deconvolution.  The same data created by MsPASS deconvolution operators can then, 
after QC, be used as input to this package.  

The core algorithms planned for parallel_pwmig are the same as the original pwmig.   The difference is only that this version should 
run faster and be scalable to more processors and run much more efficiently on large clusters.   With luck it may one day even run on 
cloud systems. 
