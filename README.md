# ohanas_version_trap
This version of the ion trap code developed by Ohana Rodrigues and Walter Freeman to reproduce the stability bands plot or trapped ion animations.  

*************************************************************************************************************************************
*************************************************************************************************************************************
******************************************************* General instructions  *******************************************************
*************************************************************************************************************************************
*************************************************************************************************************************************

In order to use the following code, you need to have anim (from Walter Freeman) and root installed. You can find more about Anim here...
https://walterfreeman.github.io/phys307/notes/anim.html

ionTrap.h is where the class methods are. 
To reproduce the stability plots, use the BuildRootFile() to get the root output that you will later use as input to trapPlot.C. 
To reproduce the animations, use PrintForAnime() to get the output to be fed to anim. 
You will have to select one or the other in main.cpp before compiling and running it. 

vector.h is just a header that defines vector operations.

Edit the Makefile to point to your own root libs. 
