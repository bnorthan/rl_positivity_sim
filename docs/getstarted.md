## rl_positivy_sim with new options and auto reports

This repo is forked from James Manton's [rl_positivity_sim](https://github.com/jdmanton/rl_positivity_sim)

Recently Weisong Zhao has started another fork of this repo to talk about the theory behind RL [see here}(https://weisongzhao.github.io/rl_positivity_sim/)

I realized I had also forked the same repo some time ago, because it was super useful to quickly examine a question I had about Weiner deconvolution.  

In the future it could be cool to produce a python version of this repo.  I mention that because if collaboration is the end goal, python may be better to facilitate that. 

For now I am staying in MATLAB though.  There is already nice code as a starting point, + MATLAB has several deconvolution algorithms with fancy options that would be interesting to test.  

So briefly here are the new features:

*  Reports are now automatically generated as markdown files [see](https://github.com/bnorthan/rl_positivity_sim/blob/master/reports/classic_rl_points_1000/report.md) for an example.
*  Can choose between method.
  1. classic RL (James Manton's bare bones implementaiton true to the orginal RL)
  2. deconlucy (MATLAB implementation that includes acceleration and options for noise supression)
  3. deconvblind (MATLAB implemtation similar to above that also updates the PSF (be careful with this option, the PSF almost never converges to the true PSF and results are ussually inferior)). 
*  Can choose between multiple simulation types
  1.  Lines same as in orginal code, 2 sets of lines forming a corss.
  2.  Points - Pairs of points starting from 3 pixels apart and increasing by 2 pixels apart each set.
  3.  circles varying size - a set of circles of increasing size
  4.  circles varying intensity - a set of circles of varying intensity. 
